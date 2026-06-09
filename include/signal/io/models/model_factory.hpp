// SPDX-License-Identifier: MIT
// Factory helpers for k-mer pore models.

#pragma once

#include <memory>
#include <string>

#include "signal/io/models/model.hpp"

namespace panomap::io {

using ModelPtr = std::unique_ptr<KmerModel>;

// Load a built-in model by name (e.g., "r10.4"). Returns nullptr if unknown.
ModelPtr load_builtin_model(const std::string& name);

// Load a model from file (to be implemented in later stages).
ModelPtr load_model_from_file(const std::string& path);

}  // namespace panomap::io
