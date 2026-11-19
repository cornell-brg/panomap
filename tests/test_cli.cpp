// SPDX-License-Identifier: MIT
#include <cstdlib>
#include <doctest/doctest.h>
#include <filesystem>
#include <string>
#include <vector>

#include "cli/parse.hpp"

// This is a simple smoke test to check that the CLI commands run without crashing.
// It does not check for correctness of the output.

TEST_CASE("CLI smoke test: panomap index") {
  const std::string output_base = std::filesystem::temp_directory_path() / "panomap_cli_test_index";
  const std::string index_file = output_base + ".pirx";

  // Clean up any previous run
  std::filesystem::remove(index_file);

  const std::string graph_path = "../tests/data/graphs/sample.gfa";
  const std::string command = "./panomap index --output " + output_base + " " + graph_path;

  // The command needs to be run from the build directory.
  // The test runner (ctest) runs from panomap/build, so the relative path is correct.
  int ret = std::system(command.c_str());

  CHECK(ret == 0);

  // Check that the simple index file was created (default is now simple backend)
  CHECK(std::filesystem::exists(index_file));

  // Clean up
  std::filesystem::remove(index_file);
}

TEST_CASE("expand_preset puts preset arguments first so later flags override them") {
  const panomap::cli::PresetTable presets = {{"r10", {"--diff", "0", "--no-anchor-merge"}}};
  std::string error;

  std::vector<std::string> a = {"--threads", "4", "-x", "r10", "--diff", "0.1", "reads.blow5"};
  REQUIRE(panomap::cli::expand_preset(a, presets, &error));
  CHECK(a == std::vector<std::string>{"--diff", "0", "--no-anchor-merge", "--threads", "4",
                                      "--diff", "0.1", "reads.blow5"});

  std::vector<std::string> b = {"--preset=r10"};
  REQUIRE(panomap::cli::expand_preset(b, presets, &error));
  CHECK(b.size() == 3);

  std::vector<std::string> c = {"-xr10", "in.gfa"};
  REQUIRE(panomap::cli::expand_preset(c, presets, &error));
  CHECK(c.back() == "in.gfa");

  std::vector<std::string> d = {"-x", "r9"};
  CHECK_FALSE(panomap::cli::expand_preset(d, presets, &error));
  CHECK(error == "unknown preset 'r9'");

  std::vector<std::string> e = {"-x"};
  CHECK_FALSE(panomap::cli::expand_preset(e, presets, &error));
}
