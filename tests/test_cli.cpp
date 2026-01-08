// SPDX-License-Identifier: MIT
#include <doctest/doctest.h>

#include <cstdlib>
#include <filesystem>
#include <string>

// This is a simple smoke test to check that the CLI commands run without crashing.
// It does not check for correctness of the output.

TEST_CASE("CLI smoke test: piru index") {
    const std::string output_base = std::filesystem::temp_directory_path() / "piru_cli_test_index";
    const std::string index_file = output_base + ".pirx";

    // Clean up any previous run
    std::filesystem::remove(index_file);

    const std::string graph_path = "../tests/data/graphs/sample.gfa";
    const std::string command = "./piru index --output " + output_base + " " + graph_path;

    // The command needs to be run from the build directory.
    // The test runner (ctest) runs from piru/build, so the relative path is correct.
    int ret = std::system(command.c_str());

    CHECK(ret == 0);

    // Check that the simple index file was created (default is now simple backend)
    CHECK(std::filesystem::exists(index_file));

    // Clean up
    std::filesystem::remove(index_file);
}
