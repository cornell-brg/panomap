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
                                   bool primary_only, const std::vector<float>* node_1d_coords,
                                   const std::vector<std::uint32_t>* component_ids) {
  const std::string ext = extension_of(path);
  if (ext == "paf") {
    return std::make_unique<PafWriter>(path);
  }
  // Default to GAF for .gaf or any unrecognized extension
  GafWriterConfig config;
  config.primary_only = primary_only;
  config.node_1d_coords = node_1d_coords;
  config.component_ids = component_ids;
  return std::make_unique<GafWriter>(path, graph, config);
}

ResultWriterPtr make_result_writer_stdout(const index::FlatGraph& graph, bool primary_only,
                                          const std::vector<float>* node_1d_coords,
                                          const std::vector<std::uint32_t>* component_ids) {
  GafWriterConfig config;
  config.primary_only = primary_only;
  config.node_1d_coords = node_1d_coords;
  config.component_ids = component_ids;
  return std::make_unique<GafWriter>(graph, config);
}

}  // namespace piru::io
