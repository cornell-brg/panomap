/**
 * map.hpp
 *
 * CLI entry point for `panomap map`.
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <string>
#include <vector>

/** Handle the `panomap map` subcommand. */
int handle_map(const std::vector<std::string>& args);
