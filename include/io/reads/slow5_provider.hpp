// SPDX-License-Identifier: MIT
// Slow5/Blow5 read provider (stub unless slow5lib is linked).

#pragma once

#include <memory>
#include <string>

#include "io/reads/read_provider.hpp"

#ifdef PIRU_HAS_SLOW5
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wc++20-extensions"
#endif
#include <slow5/slow5.h>
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#endif

namespace piru::io {

class Slow5Provider : public ReadProvider {
public:
    explicit Slow5Provider(const std::string& path);
    ~Slow5Provider() override;

    bool get_next(RawRead& read) override;
    void reset() override;
    std::string get_format_name() const override;

private:
    std::string path_;
    void open();
    void close();
#ifdef PIRU_HAS_SLOW5
    slow5_file_t* fp_{nullptr};
    slow5_rec_t* rec_{nullptr};
#endif
    bool warned_{false};
};

using Slow5ProviderPtr = std::unique_ptr<Slow5Provider>;

}  // namespace piru::io
