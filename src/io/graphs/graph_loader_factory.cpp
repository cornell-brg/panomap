#include "io/graphs/graph_loader_factory.hpp"

#include <algorithm>
#include <cctype>

#include "io/graphs/gfa_loader.hpp"
#include "io/graphs/vg_loader.hpp"
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

GraphLoaderPtr make_graph_loader(const std::string& path) {
    const std::string ext = extension_of(path);
    if (ext == "gfa") {
        return std::make_unique<GfaLoader>(path);
    }
    if (ext == "vg" || ext == "json") {
#ifdef PIRU_HAS_LIBVGIO
        return std::make_unique<VgLoader>(path);
#else
        LOG_ERROR("libvgio not available; vg format unsupported for '" + path + "'");
        return nullptr;
#endif
    }
    LOG_ERROR("Unsupported graph format for '" + path + "'");
    return nullptr;
}

}  // namespace piru::io
