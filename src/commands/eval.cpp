#include "commands/eval.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "cli/parse.hpp"
#include "util/logging.hpp"
#include "version.hpp"

namespace {

// Ground truth info parsed from squigulator read names
// Format: read_id!reference!start!end!strand
struct GroundTruth {
    std::string read_id;
    std::string reference;
    int64_t start{0};
    int64_t end{0};
    char strand{'+'};
    bool valid{false};
};

// Mapping result from PAF
struct MappingResult {
    std::string query_name;
    int64_t query_length{0};
    int64_t query_start{0};
    int64_t query_end{0};
    char strand{'+'};
    std::string target_name;
    int64_t target_length{0};
    int64_t target_start{0};
    int64_t target_end{0};
    int64_t matches{0};
    int64_t block_length{0};
    int mapq{0};
};

// Evaluation result per read
struct EvalResult {
    std::string read_id;
    bool mapped{false};
    bool correct_target{false};
    bool correct_strand{false};
    bool correct_position{false};
    int64_t position_error{0};      // Distance from ground truth start
    double best_overlap_frac{0.0};  // Best overlap fraction with ground truth
    GroundTruth truth;
    MappingResult mapping;
};

// Compute overlap between two intervals [a_start, a_end) and [b_start, b_end)
int64_t computeOverlap(int64_t a_start, int64_t a_end, int64_t b_start, int64_t b_end) {
    int64_t overlap_start = std::max(a_start, b_start);
    int64_t overlap_end = std::min(a_end, b_end);
    return std::max(int64_t(0), overlap_end - overlap_start);
}

// Compute overlap fraction relative to mapping result length
// (what fraction of the mapping falls within the truth region)
double computeOverlapFraction(int64_t truth_start, int64_t truth_end, int64_t mapping_start,
                              int64_t mapping_end) {
    int64_t mapping_len = mapping_end - mapping_start;
    if (mapping_len <= 0) return 0.0;
    int64_t overlap = computeOverlap(truth_start, truth_end, mapping_start, mapping_end);
    return static_cast<double>(overlap) / static_cast<double>(mapping_len);
}

// Parse ground truth from squigulator read name
// Format: S1_1!gi|568815592:32578768-32589835!128!923!-
GroundTruth parseGroundTruth(const std::string& read_name) {
    GroundTruth gt;
    gt.read_id = read_name;

    // Split by '!'
    std::vector<std::string> parts;
    std::istringstream iss(read_name);
    std::string part;
    while (std::getline(iss, part, '!')) {
        parts.push_back(part);
    }

    if (parts.size() >= 5) {
        gt.reference = parts[1];
        try {
            gt.start = std::stoll(parts[2]);
            gt.end = std::stoll(parts[3]);
            gt.strand = parts[4].empty() ? '+' : parts[4][0];
            gt.valid = true;
        } catch (...) {
            gt.valid = false;
        }
    }

    return gt;
}

// Parse a PAF line
MappingResult parsePafLine(const std::string& line) {
    MappingResult result;
    std::istringstream iss(line);

    std::string strand_str;
    iss >> result.query_name >> result.query_length >> result.query_start >> result.query_end >>
        strand_str >> result.target_name >> result.target_length >> result.target_start >>
        result.target_end >> result.matches >> result.block_length >> result.mapq;

    result.strand = strand_str.empty() ? '+' : strand_str[0];
    return result;
}

// Load all mappings from a PAF file, keyed by query name
// Returns ALL mappings per read (not just primary) for haplotype-aware evaluation
std::unordered_map<std::string, std::vector<MappingResult>> loadPaf(const std::string& path) {
    std::unordered_map<std::string, std::vector<MappingResult>> mappings;
    std::ifstream in(path);
    if (!in) {
        LOG_ERROR("Failed to open PAF file: " + path);
        return mappings;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto result = parsePafLine(line);
        mappings[result.query_name].push_back(result);
    }
    return mappings;
}

// Load ground truth read names from truth PAF
std::vector<GroundTruth> loadGroundTruth(const std::string& path) {
    std::vector<GroundTruth> truths;
    std::ifstream in(path);
    if (!in) {
        LOG_ERROR("Failed to open ground truth PAF: " + path);
        return truths;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        std::string query_name;
        iss >> query_name;
        auto gt = parseGroundTruth(query_name);
        if (gt.valid) {
            truths.push_back(gt);
        }
    }
    return truths;
}

}  // namespace

int handle_eval(const std::vector<std::string>& args) {
    piru::cli::Parsed parsed;
    piru::cli::ParseConfig config;
    config.usage = "Usage: piru eval -t <truth.paf> -c <calls.paf> [options]";
    config.options = {
        {'h', "help", false, "Show help"},
        {'v', "verbose", false, "Show per-read results (pass/fail with details)"},
        {'t', "truth", true, "Ground truth PAF (from squigulator -c)"},
        {'c', "calls", true, "Mapping results PAF (from piru map)"},
        {'p', "tolerance", true, "Position tolerance in bp (default: 100)"},
        {'m', "min-overlap", true, "Min overlap fraction with truth (0.0-1.0, default: 0.5)"},
        {'l', "min-length", true, "Min read length in bp to include (default: 0, no filter)"},
        {'o', "output", true, "Output file (default: stdout)"},
        {'f', "format", true, "Output format: summary, tsv, json (default: summary)"},
    };
    config.on_error = [](const std::string&) { std::cerr << "eval: invalid option\n"; };

    if (!piru::cli::parse_args(args, config, parsed)) {
        piru::cli::print_help(config, std::cerr);
        return 1;
    }
    if (parsed.values.count("help")) {
        piru::cli::print_help(config, std::cout);
        return 0;
    }

    if (!parsed.values.count("truth") || !parsed.values.count("calls")) {
        LOG_ERROR("Both --truth and --calls are required");
        piru::cli::print_help(config, std::cerr);
        return 1;
    }

    const std::string truth_path = parsed.values.at("truth");
    const std::string calls_path = parsed.values.at("calls");
    const int64_t tolerance =
        parsed.values.count("tolerance") ? std::stoll(parsed.values.at("tolerance")) : 100;
    const double min_overlap =
        parsed.values.count("min-overlap") ? std::stod(parsed.values.at("min-overlap")) : 0.5;
    const int64_t min_length =
        parsed.values.count("min-length") ? std::stoll(parsed.values.at("min-length")) : 0;
    const std::string format =
        parsed.values.count("format") ? parsed.values.at("format") : "summary";

    // Load data
    auto truths = loadGroundTruth(truth_path);
    auto mappings = loadPaf(calls_path);

    if (truths.empty()) {
        LOG_ERROR("No ground truth entries found in " + truth_path);
        return 1;
    }

    // Filter by min_length if specified
    int64_t filtered_count = 0;
    if (min_length > 0) {
        std::vector<GroundTruth> filtered;
        for (const auto& gt : truths) {
            int64_t len = gt.end - gt.start;
            if (len >= min_length) {
                filtered.push_back(gt);
            } else {
                filtered_count++;
            }
        }
        truths = std::move(filtered);

        if (truths.empty()) {
            LOG_ERROR("No reads remaining after min-length filter (" + std::to_string(min_length) +
                      " bp)");
            return 1;
        }
    }

    // Compute read length stats from ground truth (after filtering)
    int64_t min_len = std::numeric_limits<int64_t>::max();
    int64_t max_len = 0;
    int64_t sum_len = 0;
    for (const auto& gt : truths) {
        int64_t len = gt.end - gt.start;
        min_len = std::min(min_len, len);
        max_len = std::max(max_len, len);
        sum_len += len;
    }
    double avg_len = truths.empty() ? 0.0 : static_cast<double>(sum_len) / truths.size();

    // Evaluate each read
    // For haplotype graphs: a read is correct if ANY of its candidate mappings matches ground truth
    std::vector<EvalResult> results;
    int total = 0;
    int mapped = 0;
    int correct_target = 0;
    int correct_strand = 0;
    int correct_position = 0;
    int correct_overlap = 0;

    for (const auto& gt : truths) {
        EvalResult eval;
        eval.read_id = gt.read_id;
        eval.truth = gt;
        total++;

        auto it = mappings.find(gt.read_id);
        if (it != mappings.end() && !it->second.empty()) {
            eval.mapped = true;
            mapped++;

            // Check ALL candidates - if ANY matches, consider it correct
            int64_t best_position_error = std::numeric_limits<int64_t>::max();
            double best_overlap = 0.0;
            MappingResult best_mapping;

            for (const auto& candidate : it->second) {
                // Check target
                bool target_match =
                    (candidate.target_name.find(gt.reference) != std::string::npos ||
                     gt.reference.find(candidate.target_name) != std::string::npos ||
                     candidate.target_name == gt.reference);

                // Check strand
                bool strand_match = (candidate.strand == gt.strand);

                // Check position (start within tolerance)
                int64_t pos_err = std::abs(candidate.target_start - gt.start);
                bool pos_match = (pos_err <= tolerance);

                // Check overlap fraction
                double overlap_frac = computeOverlapFraction(
                    gt.start, gt.end, candidate.target_start, candidate.target_end);
                bool overlap_match = (overlap_frac >= min_overlap);

                // Track best candidate (highest overlap, then lowest position error)
                if (overlap_frac > best_overlap ||
                    (overlap_frac == best_overlap && pos_err < best_position_error)) {
                    best_overlap = overlap_frac;
                    best_position_error = pos_err;
                    best_mapping = candidate;
                }

                // If ANY candidate matches, mark as correct
                if (target_match) eval.correct_target = true;
                if (strand_match) eval.correct_strand = true;
                if (pos_match) eval.correct_position = true;
                if (overlap_match)
                    eval.best_overlap_frac = std::max(eval.best_overlap_frac, overlap_frac);
            }

            eval.mapping = best_mapping;
            eval.position_error = best_position_error;
            eval.best_overlap_frac = best_overlap;

            if (eval.correct_target) correct_target++;
            if (eval.correct_strand) correct_strand++;
            if (eval.correct_position) correct_position++;
            if (eval.best_overlap_frac >= min_overlap) correct_overlap++;
        }

        results.push_back(eval);
    }

    // Output results
    std::ostream* out = &std::cout;
    std::ofstream file_out;
    if (parsed.values.count("output")) {
        file_out.open(parsed.values.at("output"));
        if (!file_out) {
            LOG_ERROR("Failed to open output file: " + parsed.values.at("output"));
            return 1;
        }
        out = &file_out;
    }

    if (format == "summary") {
        *out << "=== Evaluation Summary ===\n";
        *out << "Ground truth: " << truth_path << "\n";
        *out << "Calls: " << calls_path << "\n";
        *out << "Position tolerance: " << tolerance << " bp\n";
        *out << "Min overlap: " << std::fixed << std::setprecision(0) << (min_overlap * 100)
             << "%\n";
        if (min_length > 0) {
            *out << "Min length filter: " << min_length << " bp (" << filtered_count
                 << " reads excluded)\n";
        }
        *out << "\n";
        *out << "Read length (from truth): min=" << min_len << ", max=" << max_len
             << ", avg=" << std::fixed << std::setprecision(0) << avg_len << " bp\n";
        *out << "\n";
        *out << "Total reads:      " << total << "\n";
        *out << "Mapped:           " << mapped << " (" << std::fixed << std::setprecision(1)
             << (100.0 * mapped / total) << "%)\n";
        *out << "Correct target:   " << correct_target << " (" << std::fixed << std::setprecision(1)
             << (100.0 * correct_target / total) << "%)\n";
        *out << "Correct strand:   " << correct_strand << " (" << std::fixed << std::setprecision(1)
             << (100.0 * correct_strand / total) << "%)\n";
        *out << "Correct position: " << correct_position << " (" << std::fixed
             << std::setprecision(1) << (100.0 * correct_position / total) << "%, within "
             << tolerance << "bp)\n";
        *out << "Correct overlap:  " << correct_overlap << " (" << std::fixed
             << std::setprecision(1) << (100.0 * correct_overlap / total)
             << "%, >=" << std::setprecision(0) << (min_overlap * 100) << "% overlap)\n";

        // Verbose: per-read results
        if (parsed.values.count("verbose")) {
            *out << "\n=== Per-Read Results ===\n";
            for (const auto& r : results) {
                int64_t read_len = r.truth.end - r.truth.start;
                bool pass = r.mapped && r.best_overlap_frac >= min_overlap;
                *out << (pass ? "PASS" : "FAIL") << "  " << r.read_id << " (len=" << read_len
                     << "bp)";
                if (!r.mapped) {
                    *out << "  UNMAPPED";
                } else {
                    *out << "  mapped=" << r.mapping.target_name << ":" << r.mapping.target_start
                         << "-" << r.mapping.target_end << "(" << r.mapping.strand << ")"
                         << "  truth=" << r.truth.reference << ":" << r.truth.start << "-"
                         << r.truth.end << "(" << r.truth.strand << ")"
                         << "  overlap=" << std::fixed << std::setprecision(1)
                         << (r.best_overlap_frac * 100) << "%"
                         << "  pos_err=" << r.position_error << "bp";
                }
                *out << "\n";
            }
        } else {
            // Default: only show problem reads
            bool has_issues = false;
            for (const auto& r : results) {
                if (!r.mapped || r.best_overlap_frac < min_overlap) {
                    if (!has_issues) {
                        *out << "\n=== Problem Reads ===\n";
                        has_issues = true;
                    }
                    int64_t read_len = r.truth.end - r.truth.start;
                    *out << r.read_id << " (len=" << read_len << "bp): ";
                    if (!r.mapped) {
                        *out << "UNMAPPED";
                    } else {
                        *out << "overlap=" << std::fixed << std::setprecision(1)
                             << (r.best_overlap_frac * 100) << "%";
                        *out << ", pos_err=" << r.position_error << "bp";
                        if (!r.correct_strand) *out << ", wrong_strand";
                    }
                    *out << "\n";
                }
            }
        }
    } else if (format == "tsv") {
        *out << "read_id\tread_length\tmapped\tcorrect_target\tcorrect_strand\tcorrect_"
                "position\toverlap_frac\t"
             << "position_error\ttruth_start\ttruth_end\tmapping_start\tmapping_end\n";
        for (const auto& r : results) {
            int64_t read_len = r.truth.end - r.truth.start;
            *out << r.read_id << "\t" << read_len << "\t" << (r.mapped ? "1" : "0") << "\t"
                 << (r.correct_target ? "1" : "0") << "\t" << (r.correct_strand ? "1" : "0") << "\t"
                 << (r.correct_position ? "1" : "0") << "\t" << std::fixed << std::setprecision(3)
                 << r.best_overlap_frac << "\t" << r.position_error << "\t" << r.truth.start << "\t"
                 << r.truth.end << "\t"
                 << (r.mapped ? std::to_string(r.mapping.target_start) : "NA") << "\t"
                 << (r.mapped ? std::to_string(r.mapping.target_end) : "NA") << "\n";
        }
    } else if (format == "json") {
        *out << "{\n";
        *out << "  \"summary\": {\n";
        *out << "    \"total\": " << total << ",\n";
        *out << "    \"mapped\": " << mapped << ",\n";
        *out << "    \"correct_target\": " << correct_target << ",\n";
        *out << "    \"correct_strand\": " << correct_strand << ",\n";
        *out << "    \"correct_position\": " << correct_position << ",\n";
        *out << "    \"correct_overlap\": " << correct_overlap << ",\n";
        *out << "    \"tolerance_bp\": " << tolerance << ",\n";
        *out << "    \"min_overlap\": " << std::fixed << std::setprecision(2) << min_overlap
             << ",\n";
        *out << "    \"min_length_bp\": " << min_length << ",\n";
        *out << "    \"filtered_count\": " << filtered_count << ",\n";
        *out << "    \"mapping_rate\": " << std::fixed << std::setprecision(4)
             << (1.0 * mapped / total) << ",\n";
        *out << "    \"position_accuracy\": " << std::fixed << std::setprecision(4)
             << (1.0 * correct_position / total) << ",\n";
        *out << "    \"overlap_accuracy\": " << std::fixed << std::setprecision(4)
             << (1.0 * correct_overlap / total) << ",\n";
        *out << "    \"read_length_min_bp\": " << min_len << ",\n";
        *out << "    \"read_length_max_bp\": " << max_len << ",\n";
        *out << "    \"read_length_avg_bp\": " << std::fixed << std::setprecision(0) << avg_len
             << "\n";
        *out << "  }\n";
        *out << "}\n";
    } else {
        LOG_ERROR("Unknown format: " + format);
        return 1;
    }

    return 0;
}
