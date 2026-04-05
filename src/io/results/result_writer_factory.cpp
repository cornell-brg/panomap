#include "io/results/result_writer_factory.hpp"

#include <algorithm>

#include "io/results/gaf_writer.hpp"
#include "io/results/paf_writer.hpp"
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

ResultWriterPtr make_result_writer(const std::string& path, const index::FlatGraph& graph,
                                   bool primary_only) {
  return make_result_writer(path, extension_of(path), graph, primary_only);
}

ResultWriterPtr make_result_writer(const std::string& path, const std::string& format,
                                   const index::FlatGraph& graph, bool primary_only) {
  const std::string fmt = to_lower(format);
  if (fmt == "paf") {
    return std::make_unique<PafWriter>(path);
  }
  if (fmt == "gaf") {
    GafWriterConfig config;
    config.primary_only = primary_only;
    return std::make_unique<GafWriter>(path, graph, config);
  }
  LOG_ERROR("Unsupported result format '" + format + "' for '" + path + "'");
  return nullptr;
}

}  // namespace piru::io
