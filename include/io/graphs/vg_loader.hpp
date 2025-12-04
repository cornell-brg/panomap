// SPDX-License-Identifier: MIT
// Loader for vg graphs (expects JSON exported via `vg view -J`).

#pragma once

#include <memory>
#include <string>

#include "io/graphs/graph_loader.hpp"

namespace piru::io {

class VgLoader : public GraphLoader {
public:
    explicit VgLoader(std::string path);
    ~VgLoader() override = default;

    bool load(ImportedGraph& graph) override;
    std::string get_format_name() const override;

private:
    std::string path_;
    bool warned_{false};
};

using VgLoaderPtr = std::unique_ptr<VgLoader>;

}  // namespace piru::io
