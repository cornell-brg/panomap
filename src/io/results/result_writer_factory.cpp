#include "io/results/result_writer_factory.hpp"

#include <algorithm>

#include "io/results/gaf_writer.hpp"
#include "io/results/gam_writer.hpp"
#include "io/results/json_writer.hpp"
#include "util/logging.hpp"

namespace piru::io {

namespace {

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}

std::string extension_of(const std::string& path) {
    const auto pos = path.find_last_of('.');
    if (pos == std::string::npos) return "";
    return to_lower(path.substr(pos + 1));
}

}  // namespace

ResultWriterPtr make_result_writer(const std::string& path) {
    const std::string ext = extension_of(path);
    if (ext == "gaf") {
        return std::make_unique<GafWriter>(path);
    }
    if (ext == "gam") {
        return std::make_unique<GamWriter>(path);
    }
    if (ext == "json") {
        return std::make_unique<JsonWriter>(path);
    }
    LOG_ERROR("Unsupported result format for '" + path + "'");
    return nullptr;
}

}  // namespace piru::io
