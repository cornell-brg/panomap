// SPDX-License-Identifier: MIT

#include "signal/event_pipelines/scrappie_event_pipeline.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace piru::signal {

namespace {

constexpr float kMadScale = 1.4826f;
constexpr float kMinMad = 1e-10f;

float median(std::vector<float>& vec) {
    const size_t n = vec.size();
    if (n == 0) {
        return 0.0f;
    }
    const size_t half = n / 2;
    std::nth_element(vec.begin(), vec.begin() + half, vec.end());
    float med = vec[half];
    if ((n % 2) == 0) {
        std::nth_element(vec.begin(), vec.begin() + half - 1, vec.end());
        med = (med + vec[half - 1]) * 0.5f;
    }
    return med;
}

float mad(std::vector<float>& vec, float med) {
    for (auto& v : vec) {
        v = std::abs(v - med);
    }
    return kMadScale * median(vec);
}

struct RawView {
    const std::vector<float>* raw{nullptr};
    int start{0};
    int end{0};
};

RawView trim_raw_by_mad(const std::vector<float>& raw, int chunk_size, float quantile) {
    const int n = static_cast<int>(raw.size());
    if (n == 0 || chunk_size <= 0) {
        return RawView{&raw, 0, 0};
    }
    if (n < chunk_size) {
        return RawView{&raw, 0, n};
    }

    const int nchunk = n / chunk_size;
    if (nchunk == 0) {
        return RawView{&raw, 0, n};
    }

    std::vector<float> mad_arr;
    mad_arr.reserve(nchunk);
    std::vector<float> chunk(chunk_size);

    for (int i = 0; i < nchunk; ++i) {
        std::copy(raw.begin() + i * chunk_size, raw.begin() + (i + 1) * chunk_size, chunk.begin());
        float m = median(chunk);
        mad_arr.push_back(mad(chunk, m));
    }

    if (mad_arr.empty()) {
        return RawView{&raw, 0, n};
    }

    std::vector<float> sorted_mad = mad_arr;
    std::sort(sorted_mad.begin(), sorted_mad.end());

    quantile = std::clamp(quantile, 0.0f, 1.0f);
    const int thresh_idx =
        std::clamp(static_cast<int>(quantile * nchunk), 0, std::max(nchunk - 1, 0));
    const float thresh = sorted_mad[thresh_idx];

    if (thresh < kMinMad) {
        return RawView{&raw, 0, n};
    }
    if (sorted_mad.front() == sorted_mad.back()) {
        return RawView{&raw, 0, n};
    }

    const bool use_gte = (thresh == sorted_mad.back());

    int start = 0;
    int end = n;
    for (int i = 0; i < nchunk; ++i) {
        if (use_gte ? (mad_arr[i] >= thresh) : (mad_arr[i] > thresh)) {
            break;
        }
        start += chunk_size;
    }
    for (int i = nchunk - 1; i >= 0; --i) {
        if (use_gte ? (mad_arr[i] >= thresh) : (mad_arr[i] > thresh)) {
            break;
        }
        end -= chunk_size;
    }

    return RawView{&raw, start, end};
}

void compute_sum_sumsq(const std::vector<float>& signal, std::vector<double>& sum,
                       std::vector<double>& sumsq) {
    const int n = static_cast<int>(signal.size());
    sum.assign(n + 1, 0.0);
    sumsq.assign(n + 1, 0.0);
    for (int i = 0; i < n; ++i) {
        sum[i + 1] = sum[i] + signal[i];
        sumsq[i + 1] = sumsq[i] + signal[i] * signal[i];
    }
}

std::vector<float> compute_tstat(const std::vector<double>& sum, const std::vector<double>& sumsq,
                                 int length, int wlen) {
    std::vector<float> tstat(length, 0.0f);
    const float eps = std::numeric_limits<float>::min();  // Match Scrappie's FLT_MIN

    if (length < 2 * wlen || wlen < 2) return tstat;

    for (int i = wlen; i < length - wlen; ++i) {
        const double sum1 = sum[i] - sum[i - wlen];
        const double sum2 = sum[i + wlen] - sum[i];
        const double sq1 = sumsq[i] - sumsq[i - wlen];
        const double sq2 = sumsq[i + wlen] - sumsq[i];

        const double mean1 = sum1 / wlen;
        const double mean2 = sum2 / wlen;
        const double var1 = sq1 / wlen - mean1 * mean1;
        const double var2 = sq2 / wlen - mean2 * mean2;

        const double combined_var = std::max(var1 + var2, static_cast<double>(eps));
        tstat[i] = static_cast<float>(std::abs(mean2 - mean1) / std::sqrt(combined_var / wlen));
    }

    return tstat;
}

std::vector<int> short_long_peak_detector(const std::vector<float>& t1,
                                          const std::vector<float>& t2, int window_length1,
                                          int window_length2, float threshold1, float threshold2,
                                          float peak_height) {
    std::vector<int> peak_pos;
    const int length = static_cast<int>(t1.size());

    float peak_val[2] = {std::numeric_limits<float>::infinity(),
                         std::numeric_limits<float>::infinity()};
    int peak_idx[2] = {-1, -1};
    int masked_to[2] = {0, 0};
    bool valid_peak[2] = {false, false};
    const float thresholds[2] = {threshold1, threshold2};
    const int wlens[2] = {window_length1, window_length2};
    const std::vector<float>* tstats[2] = {&t1, &t2};

    for (int i = 0; i < length; ++i) {
        for (int k = 0; k < 2; ++k) {
            if (masked_to[k] >= i) continue;
            const float curr = (*tstats[k])[i];

            if (peak_idx[k] == -1) {
                if (curr < peak_val[k]) {
                    peak_val[k] = curr;
                } else if (curr - peak_val[k] > peak_height) {
                    peak_val[k] = curr;
                    peak_idx[k] = i;
                }
            } else {
                if (curr > peak_val[k]) {
                    peak_val[k] = curr;
                    peak_idx[k] = i;
                }
                if (k == 0 && peak_val[0] > thresholds[0]) {
                    masked_to[1] = peak_idx[0] + wlens[0];
                    peak_idx[1] = -1;
                    peak_val[1] = std::numeric_limits<float>::infinity();
                    valid_peak[1] = false;
                }
                if (peak_val[k] - curr > peak_height && peak_val[k] > thresholds[k]) {
                    valid_peak[k] = true;
                }
                if (valid_peak[k] && (i - peak_idx[k]) > wlens[k] / 2) {
                    peak_pos.push_back(peak_idx[k]);
                    peak_idx[k] = -1;
                    peak_val[k] = curr;
                    valid_peak[k] = false;
                }
            }
        }
    }

    return peak_pos;
}

// Create event means from peak positions using prefix sums.
std::vector<float> create_event_means(const std::vector<int>& peaks, const std::vector<double>& sum,
                                      int nsample) {
    std::vector<float> means;
    if (peaks.empty()) {
        if (nsample > 0) {
            means.push_back(static_cast<float>(sum[nsample] / nsample));
        }
        return means;
    }

    // First event: [0, peaks[0])
    if (peaks.front() > 0) {
        means.push_back(static_cast<float>(sum[peaks.front()] / peaks.front()));
    }

    // Middle events: [peaks[i-1], peaks[i])
    for (size_t i = 1; i < peaks.size(); ++i) {
        int length = peaks[i] - peaks[i - 1];
        if (length > 0) {
            means.push_back(static_cast<float>((sum[peaks[i]] - sum[peaks[i - 1]]) / length));
        }
    }

    // Last event: [peaks.back(), nsample)
    int last_length = nsample - peaks.back();
    if (last_length > 0) {
        means.push_back(static_cast<float>((sum[nsample] - sum[peaks.back()]) / last_length));
    }

    return means;
}

// Apply z-score normalization to event means.
std::vector<float> zscore_normalize(const std::vector<float>& values, bool clip, float clip_min,
                                    float clip_max) {
    if (values.empty()) {
        return {};
    }

    const double sum = std::accumulate(values.begin(), values.end(), 0.0);
    const double mean = sum / static_cast<double>(values.size());

    double var = 0.0;
    for (const auto v : values) {
        const double diff = static_cast<double>(v) - mean;
        var += diff * diff;
    }
    var /= static_cast<double>(values.size());
    const double stddev = (var > 0.0) ? std::sqrt(var) : 1.0;

    std::vector<float> normalized;
    normalized.reserve(values.size());
    for (const auto v : values) {
        float z = static_cast<float>((static_cast<double>(v) - mean) / stddev);
        if (clip) {
            z = std::clamp(z, clip_min, clip_max);
        }
        normalized.push_back(z);
    }
    return normalized;
}

// Apply median-MAD normalization to event means.
std::vector<float> median_mad_normalize(const std::vector<float>& values, bool clip, float clip_min,
                                        float clip_max) {
    if (values.empty()) {
        return {};
    }

    std::vector<float> scratch = values;
    const float med = median(scratch);
    for (auto& v : scratch) {
        v = std::abs(v - med);
    }
    float mad_val = median(scratch) * kMadScale;
    if (mad_val == 0.0f) {
        mad_val = 1.0f;
    }

    std::vector<float> normalized;
    normalized.reserve(values.size());
    for (const auto v : values) {
        float norm = (v - med) / mad_val;
        if (clip) {
            norm = std::clamp(norm, clip_min, clip_max);
        }
        normalized.push_back(norm);
    }
    return normalized;
}

}  // namespace

ScrappieEventPipeline::ScrappieEventPipeline(EventPipelineConfig config)
    : config_(std::move(config)) {}

NormalizedSignal ScrappieEventPipeline::process(const io::RawRead& read) const {
    NormalizedSignal result;
    result.sampling_rate_hz = read.sampling_rate_hz;

    // Step 1: Convert raw ADC values to picoamps
    const float raw_unit = (read.digitisation == 0.0f) ? 1.0f : (read.range / read.digitisation);
    std::vector<float> raw;
    raw.reserve(read.raw_signal.size());
    for (auto value : read.raw_signal) {
        raw.push_back((static_cast<float>(value) + read.offset) * raw_unit);
    }

    // Step 2: Trim signal using MAD-based variance segmentation
    auto view = trim_raw_by_mad(raw, config_.varseg_chunk, config_.varseg_thresh);
    view.start = std::clamp(view.start + config_.trim_start, 0, static_cast<int>(raw.size()));
    view.end = std::clamp(view.end - config_.trim_end, 0, static_cast<int>(raw.size()));
    if (view.start >= view.end || view.raw == nullptr) {
        return result;
    }

    std::vector<float> clipped(view.raw->begin() + view.start, view.raw->begin() + view.end);
    const int nsample = static_cast<int>(clipped.size());

    // Step 3: Compute prefix sums for t-stat calculation
    std::vector<double> sum, sumsq;
    compute_sum_sumsq(clipped, sum, sumsq);

    // Step 4: Compute t-statistics with short and long windows
    auto t1 = compute_tstat(sum, sumsq, nsample, config_.window_length1);
    auto t2 = compute_tstat(sum, sumsq, nsample, config_.window_length2);

    // Step 5: Detect peaks using short/long peak detector
    auto peaks =
        short_long_peak_detector(t1, t2, config_.window_length1, config_.window_length2,
                                 config_.threshold1, config_.threshold2, config_.peak_height);

    // Step 6: Create event means from peaks
    auto event_means = create_event_means(peaks, sum, nsample);

    // Step 7: Normalize event means
    if (config_.norm_method == "median_mad") {
        result.samples = median_mad_normalize(event_means, config_.clip_outliers, config_.clip_min,
                                              config_.clip_max);
    } else {
        // Default to zscore
        result.samples = zscore_normalize(event_means, config_.clip_outliers, config_.clip_min,
                                          config_.clip_max);
    }

    return result;
}

const EventPipelineConfig& ScrappieEventPipeline::config() const { return config_; }

std::string ScrappieEventPipeline::name() const { return "scrappie"; }

}  // namespace piru::signal
