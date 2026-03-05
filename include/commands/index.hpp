/**
 * index.hpp
 *
 * CLI entry point for `piru index`.
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <string>
#include <vector>

/** Handle the `piru index` subcommand. */
int handle_index(const std::vector<std::string>& args);
