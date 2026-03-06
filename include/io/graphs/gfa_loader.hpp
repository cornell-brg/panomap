// SPDX-License-Identifier: MIT
// Loader for GFA-formatted graphs.

#pragma once

#include <memory>
#include <string>

#include "io/graphs/graph_loader.hpp"

namespace piru::io {

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

}  // namespace piru::io
