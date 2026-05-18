#include <doctest/doctest.h>

#include "core/index/flat_graph.hpp"

using namespace piru::index;

TEST_CASE("FlatGraph 2-bit encoding round-trip") {
  // Build a simple graph with known sequences
  std::vector<char> seq_data = {'A', 'C', 'G', 'T', 'N', 'A', 'C', 'G'};
  std::vector<std::uint32_t> seq_offset = {0, 5};  // node 0: ACGTN, node 1: ACG
  std::vector<std::uint32_t> seq_len = {5, 3};
  std::vector<char> name_data = {'n', '1', 'n', '2'};
  std::vector<std::uint32_t> name_offset = {0, 2};
  std::vector<std::uint16_t> name_len = {2, 2};
  std::vector<std::uint32_t> edge_target;
  std::vector<std::uint32_t> out_edge_offset = {0, 0, 0};
  std::vector<std::uint32_t> step_data;
  std::vector<std::uint32_t> path_step_offset = {0};
  std::vector<std::uint32_t> path_name_offset;
  std::vector<std::uint16_t> path_name_len;
  std::vector<std::uint64_t> path_length;

  auto fg = FlatGraph::fromRawArrays(
      2, 0, std::move(seq_data), std::move(seq_offset), std::move(seq_len), std::move(name_data),
      std::move(name_offset), std::move(name_len), std::move(edge_target),
      std::move(out_edge_offset), std::move(step_data), std::move(path_step_offset),
      std::move(path_name_offset), std::move(path_name_len), std::move(path_length));

  CHECK(fg.nodeCount() == 2);
  CHECK(fg.seqLen(0) == 5);
  CHECK(fg.seqLen(1) == 3);

  // Verify 2-bit encoding
  CHECK(fg.base2bit(0, 0) == 0);  // A
  CHECK(fg.base2bit(0, 1) == 1);  // C
  CHECK(fg.base2bit(0, 2) == 2);  // G
  CHECK(fg.base2bit(0, 3) == 3);  // T
  CHECK(fg.isN(0, 4) == true);    // N
  CHECK(fg.isN(0, 0) == false);

  // Verify decode round-trip
  CHECK(fg.seqDecoded(0) == "ACGTN");
  CHECK(fg.seqDecoded(1) == "ACG");

  // Verify metadata
  CHECK(fg.name(0) == "n1");
  CHECK(fg.name(1) == "n2");
  CHECK(fg.isReverse(0) == false);
  CHECK(fg.isReverse(1) == true);
}

TEST_CASE("FlatGraph path accessors") {
  std::vector<char> seq_data = {'A', 'C', 'G', 'T'};
  std::vector<std::uint32_t> seq_offset = {0, 2};
  std::vector<std::uint32_t> seq_len = {2, 2};
  std::vector<char> name_data = {'a', 'b', 'p', '1'};
  std::vector<std::uint32_t> name_offset = {0, 1};
  std::vector<std::uint16_t> name_len = {1, 1};
  std::vector<std::uint32_t> edge_target = {1};
  std::vector<std::uint32_t> out_edge_offset = {0, 1, 1};
  std::vector<std::uint32_t> step_data = {0, 1};
  std::vector<std::uint32_t> path_step_offset = {0, 2};
  std::vector<std::uint32_t> path_name_offset = {2};
  std::vector<std::uint16_t> path_name_len = {2};
  std::vector<std::uint64_t> path_length = {4};

  auto fg = FlatGraph::fromRawArrays(
      2, 1, std::move(seq_data), std::move(seq_offset), std::move(seq_len), std::move(name_data),
      std::move(name_offset), std::move(name_len), std::move(edge_target),
      std::move(out_edge_offset), std::move(step_data), std::move(path_step_offset),
      std::move(path_name_offset), std::move(path_name_len), std::move(path_length));

  CHECK(fg.pathCount() == 1);
  CHECK(fg.pathName(0) == "p1");
  CHECK(fg.pathLength(0) == 4);
  CHECK(fg.pathStepCount(0) == 2);

  const auto* steps = fg.pathStepsBegin(0);
  CHECK(steps[0] == 0);
  CHECK(steps[1] == 1);

  // Edge check
  CHECK(fg.outDegree(0) == 1);
  CHECK(*fg.outBegin(0) == 1);
  CHECK(fg.outDegree(1) == 0);
}
