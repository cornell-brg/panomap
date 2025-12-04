#include <doctest/doctest.h>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "io/results/result.hpp"
#include "io/results/result_writer_factory.hpp"

TEST_CASE("GAF writer writes basic record") {
    const auto tmp_path =
        std::filesystem::temp_directory_path() / "piru_test_output.gaf";

    auto writer = piru::io::make_result_writer(tmp_path.string());
    REQUIRE(writer != nullptr);

    piru::io::AlignmentResult r;
    r.query_name = "read1";
    r.query_length = 100;
    r.query_start = 10;
    r.query_end = 90;
    r.strand = '+';
    r.target_path = "chr1";
    r.target_length = 1000;
    r.target_start = 100;
    r.target_end = 180;
    r.matches = 70;
    r.alignment_block_length = 80;
    r.mapq = 60;
    r.optional_fields = {"tp:A:P"};

    REQUIRE(writer->write(r));
    writer.reset();

    std::ifstream in(tmp_path);
    REQUIRE(in.good());
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        lines.push_back(line);
    }
    REQUIRE(lines.size() == 1);
    CHECK(lines.front() ==
          "read1\t100\t10\t90\t+\tchr1\t1000\t100\t180\t70\t80\t60\ttp:A:P");
}

TEST_CASE("Result writer factory produces json writer") {
    piru::io::AlignmentResult r;
    r.query_name = "q";
    r.query_length = 10;
    r.query_sequence = "ACGTACGTAA";
    r.target_path = "t";
    r.mappings.push_back({1, 0, false, {{4, 4, ""}}});

    auto json_writer = piru::io::make_result_writer("output.json");
    REQUIRE(json_writer != nullptr);
    CHECK(json_writer->write(r));
}

TEST_CASE("Result writer factory produces gam writer") {
    piru::io::AlignmentResult r;
    r.query_name = "q";
    r.query_length = 10;
    r.query_sequence = "ACGTACGTAA";
    r.target_path = "t";
    r.mappings.push_back({1, 0, false, {{4, 4, ""}}});

    auto gam_writer = piru::io::make_result_writer("output.gam");
    REQUIRE(gam_writer != nullptr);
    CHECK(gam_writer->write(r));
}
