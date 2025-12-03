#include "io/reads/slow5_provider.hpp"

#include <cerrno>
#include <cstring>
#include <stdexcept>

#ifdef PIRU_HAS_SLOW5
#include <slow5/slow5.h>
#endif

#include "util/logging.hpp"

namespace piru::io {

Slow5Provider::Slow5Provider(const std::string& path) : path_(path) {
#ifdef PIRU_HAS_SLOW5
    open();
#endif
}

Slow5Provider::~Slow5Provider() { close(); }

void Slow5Provider::open() {
#ifdef PIRU_HAS_SLOW5
    close();
    slow5_errno = 0;
    fp_ = slow5_open(path_.c_str(), "r");
    if (fp_ == nullptr) {
        LOG_ERROR("Failed to open slow5 file: " + path_ + " (" + std::string(strerror(errno)) +
                  ")");
    }
#endif
}

void Slow5Provider::close() {
#ifdef PIRU_HAS_SLOW5
    if (rec_ != nullptr) {
        slow5_rec_free(rec_);
        rec_ = nullptr;
    }
    if (fp_ != nullptr) {
        slow5_close(fp_);
        fp_ = nullptr;
    }
#endif
}

bool Slow5Provider::get_next(RawRead& read) {
#ifdef PIRU_HAS_SLOW5
    if (fp_ == nullptr) {
        if (!warned_) {
            LOG_ERROR("slow5lib not initialized for path '" + path_ + "'");
            warned_ = true;
        }
        return false;
    }

    slow5_errno = 0;
    int ret = slow5_get_next(&rec_, fp_);
    if (ret == 0 && rec_ != nullptr) {
        read.read_id = rec_->read_id ? rec_->read_id : "";
        read.len_raw_signal = rec_->len_raw_signal;
        read.raw_signal.assign(rec_->raw_signal, rec_->raw_signal + rec_->len_raw_signal);
        read.range = static_cast<float>(rec_->range);
        read.digitisation = static_cast<float>(rec_->digitisation);
        read.offset = static_cast<float>(rec_->offset);
        read.sampling_rate_hz = static_cast<float>(rec_->sampling_rate);
        return true;
    }
    if (slow5_errno == SLOW5_ERR_EOF) {
        return false;
    }
    LOG_ERROR("slow5 read error on '" + path_ + "' (slow5_errno=" + std::to_string(slow5_errno) +
              ")");
    return false;
#else
    if (!warned_) {
        LOG_ERROR("slow5lib not linked; slow5/blow5 reading unavailable for '" + path_ + "'");
        warned_ = true;
    }
    return false;
#endif
}

void Slow5Provider::reset() {
#ifdef PIRU_HAS_SLOW5
    open();
#endif
}

std::string Slow5Provider::get_format_name() const { return "slow5"; }

}  // namespace piru::io
