// SPDX-License-Identifier: MIT
// Factory helpers to select a read provider based on format/path.

#pragma once

#include <memory>
#include <string>

#include "signal/io/reads/read_provider.hpp"

namespace piru::io {

// Attempt to choose a provider based on path/extension.
// Returns nullptr if no suitable provider is available.
ReadProviderPtr make_read_provider(const std::string& path);

}  // namespace piru::io
