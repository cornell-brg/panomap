// SPDX-License-Identifier: MIT

#include "io/reads/slow5_provider.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <stdexcept>

#include "util/logging.hpp"

namespace piru::io {

namespace {

bool has_supported_extension(const std::filesystem::path& path) {
    auto ext = path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return ext == ".slow5" || ext == ".blow5";
}

}  // namespace

Slow5Provider::Slow5Provider(std::string path)
    : input_source_(std::move(path)), filenames_(collect_input_files(input_source_)) {
    if (filenames_.empty()) {
        LOG_ERROR("No slow5/blow5 files found at " + input_source_);
        throw std::runtime_error("No slow5/blow5 files found");
    }
    open();
}

Slow5Provider::~Slow5Provider() { close(); }

std::vector<std::string> Slow5Provider::collect_input_files(const std::string& path) {
    namespace fs = std::filesystem;

    fs::path input(path);
    std::error_code ec;

    if (!fs::exists(input, ec)) {
        LOG_ERROR("SLOW5/BLOW5 input path does not exist: " + path);
        throw std::runtime_error("Input path not found");
    }

    if (fs::is_regular_file(input, ec)) {
        return {input.string()};
    }

    if (fs::is_directory(input, ec)) {
        std::vector<std::string> files;
        fs::directory_options options = fs::directory_options::skip_permission_denied;
        std::error_code dir_ec;
        for (fs::recursive_directory_iterator it(input, options, dir_ec), end;
             it != end; it.increment(dir_ec)) {
            if (dir_ec) {
                LOG_WARN("Skipping entry while scanning " + path + ": " + dir_ec.message());
                dir_ec.clear();
                continue;
            }
            const auto& entry = *it;
            std::error_code status_ec;
            if (!entry.is_regular_file(status_ec)) {
                continue;
            }
            if (!has_supported_extension(entry.path())) {
                continue;
            }
            files.push_back(entry.path().string());
        }

        std::sort(files.begin(), files.end());
        return files;
    }

    LOG_ERROR("Input path is neither a regular file nor directory: " + path);
    throw std::runtime_error("Unsupported input path");
}

void Slow5Provider::open() {
    close();
    if (current_file_index_ >= filenames_.size()) {
        fp_ = nullptr;
        return;
    }

    const auto& filename = filenames_[current_file_index_];
    fp_ = slow5_open(filename.c_str(), "r");
    if (fp_ == nullptr) {
        LOG_ERROR("Failed to open slow5 file " + filename);
        throw std::runtime_error("slow5_open failed");
    }
    LOG_INFO("Opened slow5 file " + filename);
}

void Slow5Provider::close() {
    if (fp_ != nullptr) {
        slow5_close(fp_);
        fp_ = nullptr;
    }
}

bool Slow5Provider::advance_to_next_file() {
    if (current_file_index_ + 1 >= filenames_.size()) {
        close();
        current_file_index_ = filenames_.size();
        return false;
    }
    ++current_file_index_;
    open();
    return fp_ != nullptr;
}

bool Slow5Provider::get_next(RawRead& read) {
    if (fp_ == nullptr) {
        return false;
    }

    while (true) {
        slow5_rec_t* record = nullptr;
        int ret = slow5_get_next(&record, fp_);
        if (ret < 0) {
            if (ret == SLOW5_ERR_EOF) {
                if (!advance_to_next_file()) {
                    return false;
                }
                continue;
            }
            LOG_ERROR("Non-EOF error while loading slow5 data: " + std::to_string(ret));
            throw std::runtime_error("slow5_get_next failed");
        }

        read.read_id = record->read_id ? record->read_id : "";
        read.len_raw_signal = record->len_raw_signal;
        read.raw_signal.assign(record->raw_signal, record->raw_signal + record->len_raw_signal);
        read.range = static_cast<float>(record->range);
        read.digitisation = static_cast<float>(record->digitisation);
        read.offset = static_cast<float>(record->offset);
        read.sampling_rate_hz = static_cast<float>(record->sampling_rate);

        slow5_rec_free(record);
        return true;
    }
}

void Slow5Provider::reset() {
    close();
    current_file_index_ = 0;
    open();
}

}  // namespace piru::io
