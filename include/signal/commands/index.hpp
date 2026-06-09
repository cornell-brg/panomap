/**
 * index.hpp
 *
 * CLI entry point for `panomap index`.
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <string>
#include <vector>

/** Handle the `panomap index` subcommand. */
int handle_index(const std::vector<std::string>& args);
