// SPDX-License-Identifier: MIT
#include <doctest/doctest.h>

#include <cstdlib>
#include <filesystem>
#include <string>

// This is a simple smoke test to check that the CLI commands run without crashing.
// It does not check for correctness of the output.

TEST_CASE("CLI smoke test: piru index") {
    const std::string temp_dir = std::filesystem::temp_directory_path() / "piru_cli_test_index";
    std::filesystem::create_directory(temp_dir);

    const std::string graph_path = "../tests/data/graphs/sample.gfa";
    const std::string command = "./piru index --graph-k 15 --output " + temp_dir + " " + graph_path;

    // The command needs to be run from the build directory.
    // The test runner (ctest) runs from piru/build, so the relative path is correct.
    int ret = std::system(command.c_str());

    CHECK(ret == 0);

    // Check that the three index files were created.
    CHECK(std::filesystem::exists(temp_dir + "/piru_cli_test_index.graph"));
    CHECK(std::filesystem::exists(temp_dir + "/piru_cli_test_index.signals"));
    CHECK(std::filesystem::exists(temp_dir + "/piru_cli_test_index.seeds"));

    // Clean up
    std::filesystem::remove_all(temp_dir);
}
