// SPDX-License-Identifier: MIT
#include "index/vg_transform_factory.hpp"
#include "index/path_guided_transform.hpp"
#include "util/logging.hpp"

namespace piru::index {

std::unique_ptr<VGTransform> makeVGTransform(const std::string& backend,
                                               const TransformConfig& config) {
  // Currently only path-guided transform is implemented
  if (backend == "path_guided" || backend.empty()) {
    return std::make_unique<PathGuidedTransform>();
  }

  LOG_WARN("Unknown VG transform backend: '" + backend + "', falling back to path_guided");
  return std::make_unique<PathGuidedTransform>();
}

} // namespace piru::index
