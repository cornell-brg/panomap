// SPDX-License-Identifier: MIT
// Index pipeline tests

#include <doctest/doctest.h>

#include "core/index/flat_graph.hpp"

using namespace panomap;

TEST_CASE("FlatGraph basic construction via fromRawArrays") {
  std::vector<char> seq = {'A', 'C', 'G', 'T'};
  std::vector<std::uint32_t> seq_off = {0, 2};
  std::vector<std::uint32_t> seq_len = {2, 2};
  std::vector<char> names = {'a', 'b'};
  std::vector<std::uint32_t> name_off = {0, 1};
  std::vector<std::uint16_t> name_len = {1, 1};
  std::vector<std::uint32_t> et = {1};
  std::vector<std::uint32_t> eo = {0, 1, 1};
  std::vector<std::uint32_t> sd;
  std::vector<std::uint32_t> pso = {0};
  std::vector<std::uint32_t> pno;
  std::vector<std::uint16_t> pnl;
  std::vector<std::uint64_t> pl;

  auto fg = index::FlatGraph::fromRawArrays(
      2, 0, std::move(seq), std::move(seq_off), std::move(seq_len), std::move(names),
      std::move(name_off), std::move(name_len), std::move(et), std::move(eo),
      std::move(sd), std::move(pso), std::move(pno), std::move(pnl), std::move(pl));

  CHECK(fg.nodeCount() == 2);
  CHECK(fg.edgeCount() == 1);
  CHECK(fg.seqDecoded(0) == "AC");
  CHECK(fg.seqDecoded(1) == "GT");
}
