#pragma once

#include <string>
#include <vector>

#include "cli/parse.hpp"

namespace panomap::cli {

// One preset: the index arguments and the map arguments it stands for.
// `index -x <name>` records the name in the index, and `map` applies the map
// arguments of that preset unless the command line names a preset itself.
struct Preset {
  std::string name;
  std::vector<std::string> index_args;
  std::vector<std::string> map_args;
  std::string help;
};

const std::vector<Preset>& presets();
const Preset* find_preset(const std::string& name);

PresetTable index_preset_table();
PresetTable map_preset_table();

// Help text for the -x option: one line per preset.
std::string preset_help();

}  // namespace panomap::cli
