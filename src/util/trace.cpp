// SPDX-License-Identifier: MIT

#include "util/trace.hpp"

#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <vector>

namespace piru::trace {

namespace {

bool is_tracing_active() {
  // Master switch: tracing only activates when PIRU_TRACE_DIR is set
  return std::getenv("PIRU_TRACE_DIR") != nullptr;
}

std::uint32_t parse_stages() {
  if (!is_tracing_active()) return 0;
  const char* env = std::getenv("PIRU_TRACE_STAGES");
  if (!env) return kAll;
  return static_cast<std::uint32_t>(std::strtoul(env, nullptr, 0));
}

std::vector<std::string> parse_read_filters() {
  const char* env = std::getenv("PIRU_TRACE_READS");
  std::vector<std::string> filters;
  if (!env || env[0] == '\0') return filters;
  std::istringstream ss(env);
  std::string token;
  while (std::getline(ss, token, ',')) {
    if (!token.empty()) filters.push_back(token);
  }
  return filters;
}

std::string get_trace_dir() {
  const char* env = std::getenv("PIRU_TRACE_DIR");
  if (!env) return {};
  std::string dir = env;
  std::filesystem::create_directories(dir);
  return dir;
}

}  // namespace

#ifdef PIRU_TRACE_ENABLED

std::uint32_t enabled_stages() {
  static std::uint32_t stages = parse_stages();
  return stages;
}

bool match_read(const std::string& read_id) {
  static auto filters = parse_read_filters();
  if (filters.empty()) return true;
  for (const auto& f : filters) {
    if (read_id.find(f) != std::string::npos) return true;
  }
  return false;
}

std::string trace_path(const std::string& tag, const std::string& read_id,
                        std::size_t chunk_idx) {
  static std::string dir = get_trace_dir();

  std::string gene = "unknown";
  std::string sid = "unknown";
  auto bang = read_id.find('!');
  if (bang != std::string::npos) {
    sid = read_id.substr(0, bang);
    auto colon = read_id.find(':', bang + 1);
    if (colon != std::string::npos) {
      gene = read_id.substr(bang + 1, colon - bang - 1);
    }
  }

  return dir + "/" + tag + "_" + gene + "_" + sid + "_c" + std::to_string(chunk_idx) + ".txt";
}

#else

// Stubs when tracing is disabled (never called at runtime due to if(false) in macro)
std::uint32_t enabled_stages() { return 0; }
bool match_read(const std::string&) { return false; }
std::string trace_path(const std::string&, const std::string&, std::size_t) { return {}; }

#endif

}  // namespace piru::trace
