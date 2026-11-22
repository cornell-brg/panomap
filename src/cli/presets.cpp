#include "cli/presets.hpp"

namespace panomap::cli {

namespace {

// Shared front ends. Both chemistries skip the raw diff filter, drop repeated
// tokens, and seed with minimizers (window 2).
const std::vector<std::string> kR9Index = {"-m", "r9.4", "--diff", "0", "--collapse-tokens",
                                           "--seed-type", "minimizer", "--minimizer-window", "2"};
const std::vector<std::string> kR10Index = {"-m", "r10.4", "--diff", "0", "--collapse-tokens",
                                            "--seed-type", "minimizer", "--minimizer-window", "2"};

std::vector<std::string> cat(std::vector<std::string> a, const std::vector<std::string>& b) {
  a.insert(a.end(), b.begin(), b.end());
  return a;
}

const std::vector<std::string> kR9Chain = {"--chain-pen-gap", "1.4", "--chain-pen-skip", "0.1"};
const std::vector<std::string> kR10Map = {"--chain-max-skip", "5", "--no-anchor-merge",
                                          "--chain-pen-gap", "1.2", "--map-end-mapq", "10"};

}  // namespace

const std::vector<Preset>& presets() {
  static const std::vector<Preset> kPresets = {
      {"r9", kR9Index, kR9Chain, "R9.4.1, genomes of bacterial size and up"},
      {"r9-fast", cat(kR9Index, {"--seed-k", "9"}), kR9Chain,
       "R9.4.1, about 3x faster; best with 4 or more chunks"},
      {"r9-small", kR9Index,
       {"--chain-bw", "100", "--chain-max-dist", "500", "--chain-min-score", "10",
        "--chain-pen-gap", "1.2", "--chain-pen-skip", "0.3", "--map-sc-min-anchors", "2",
        "--map-sc-ratio-lo", "0.5", "--map-sc-ratio-hi", "1.7"},
       "R9.4.1, small genomes such as viruses"},
      {"r10", kR10Index, kR10Map, "R10.4.1, genomes of bacterial size"},
      {"r10-fast", cat(kR10Index, {"--seed-k", "10"}), kR10Map,
       "R10.4.1, about 5x faster at lower recall"},
      {"r10-large", cat(kR10Index, {"--seed-k", "12"}), cat(kR10Map, {"--seed-freq-cap", "5000"}),
       "R10.4.1, large genomes such as human"},
  };
  return kPresets;
}

const Preset* find_preset(const std::string& name) {
  for (const auto& p : presets()) {
    if (p.name == name) return &p;
  }
  return nullptr;
}

PresetTable index_preset_table() {
  PresetTable t;
  for (const auto& p : presets()) t[p.name] = p.index_args;
  return t;
}

PresetTable map_preset_table() {
  PresetTable t;
  for (const auto& p : presets()) t[p.name] = p.map_args;
  return t;
}

std::string preset_help() {
  std::string s = "Preset; later flags override it. The index records it and map follows it:";
  for (const auto& p : presets()) s += "\n      " + p.name + ": " + p.help;
  return s;
}

}  // namespace panomap::cli
