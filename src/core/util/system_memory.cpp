// SPDX-License-Identifier: MIT

#include "core/util/system_memory.hpp"

#include <algorithm>
#include <cstdio>
#include <cstring>

namespace piru::util {

std::size_t getAvailableMemoryBytes() {
#ifdef __linux__
  FILE* f = fopen("/proc/meminfo", "r");
  if (!f) return 0;

  char line[256];
  while (fgets(line, sizeof(line), f)) {
    if (strncmp(line, "MemAvailable:", 13) == 0) {
      std::size_t kb = 0;
      if (sscanf(line + 13, " %zu", &kb) == 1) {
        fclose(f);
        return kb * 1024;
      }
    }
  }
  fclose(f);
#endif
  return 0;
}

std::size_t computeMemoryBudget(std::size_t max_memory, double fraction) {
  std::size_t available = getAvailableMemoryBytes();
  std::size_t budget = available > 0
                           ? static_cast<std::size_t>(static_cast<double>(available) * fraction)
                           : static_cast<std::size_t>(8ULL * 1024 * 1024 * 1024);  // 8 GB fallback

  if (max_memory > 0) {
    budget = std::min(budget, max_memory);
  }
  return budget;
}

std::string formatBytes(std::size_t bytes) {
  char buf[64];
  if (bytes >= 1024ULL * 1024 * 1024) {
    snprintf(buf, sizeof(buf), "%.1f GB", static_cast<double>(bytes) / (1024.0 * 1024 * 1024));
  } else if (bytes >= 1024ULL * 1024) {
    snprintf(buf, sizeof(buf), "%.1f MB", static_cast<double>(bytes) / (1024.0 * 1024));
  } else {
    snprintf(buf, sizeof(buf), "%.1f KB", static_cast<double>(bytes) / 1024.0);
  }
  return buf;
}

std::string formatCount(std::size_t count) {
  std::string s = std::to_string(count);
  int n = static_cast<int>(s.size());
  if (n <= 3) return s;
  std::string result;
  for (int i = 0; i < n; ++i) {
    if (i > 0 && (n - i) % 3 == 0) result += ',';
    result += s[static_cast<std::size_t>(i)];
  }
  return result;
}

std::string formatDuration(double seconds) {
  char buf[64];
  if (seconds >= 3600.0) {
    int h = static_cast<int>(seconds / 3600.0);
    int m = static_cast<int>((seconds - h * 3600.0) / 60.0);
    int s = static_cast<int>(seconds - h * 3600.0 - m * 60.0);
    snprintf(buf, sizeof(buf), "%dh %dm %ds", h, m, s);
  } else if (seconds >= 60.0) {
    int m = static_cast<int>(seconds / 60.0);
    int s = static_cast<int>(seconds - m * 60.0);
    snprintf(buf, sizeof(buf), "%dm %ds", m, s);
  } else if (seconds >= 1.0) {
    snprintf(buf, sizeof(buf), "%.1fs", seconds);
  } else {
    snprintf(buf, sizeof(buf), "%.3fs", seconds);
  }
  return buf;
}

void printProgress(std::size_t done, std::size_t total, const std::string& label) {
  if (total == 0) return;
  constexpr int kBarWidth = 40;
  double frac = static_cast<double>(done) / static_cast<double>(total);
  int filled = static_cast<int>(frac * kBarWidth);

  std::string bar(static_cast<std::size_t>(filled), '#');
  bar += std::string(static_cast<std::size_t>(kBarWidth - filled), '.');

  char buf[256];
  if (!label.empty()) {
    snprintf(buf, sizeof(buf), "\r       [%s] %zu/%zu %s", bar.c_str(), done, total, label.c_str());
  } else {
    snprintf(buf, sizeof(buf), "\r       [%s] %zu/%zu", bar.c_str(), done, total);
  }
  fprintf(stderr, "%s", buf);
  if (done >= total) fprintf(stderr, "\n");
  fflush(stderr);
}

}  // namespace piru::util
