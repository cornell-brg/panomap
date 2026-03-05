/**
 * map.hpp
 *
 * CLI entry point for `piru map`.
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <string>
#include <vector>

/** Handle the `piru map` subcommand. */
int handle_map(const std::vector<std::string>& args);
