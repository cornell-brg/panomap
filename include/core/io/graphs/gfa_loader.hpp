// SPDX-License-Identifier: MIT
// Loader for GFA-formatted graphs.

#pragma once

#include <memory>
#include <string>

#include "core/io/graphs/graph_loader.hpp"

namespace panomap::io {

class GfaLoader : public GraphLoader {
public:
  explicit GfaLoader(std::string path);
  ~GfaLoader() override = default;

  bool load(ImportedGraph& graph) override;
  std::string get_format_name() const override;

private:
  std::string path_;
};

using GfaLoaderPtr = std::unique_ptr<GfaLoader>;

}  // namespace panomap::io
