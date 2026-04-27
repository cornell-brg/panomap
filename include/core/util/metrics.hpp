#pragma once

#include <chrono>
#include <ctime>

#ifdef _WIN32
#include <psapi.h>
#include <windows.h>
#else
#include <sys/resource.h>
#endif

namespace piru {

// Wall-clock time in seconds since an arbitrary steady-clock epoch.
inline double realtime() {
  using clock = std::chrono::steady_clock;
  const auto now = clock::now().time_since_epoch();
  return std::chrono::duration_cast<std::chrono::duration<double>>(now).count();
}

// CPU time in seconds consumed by the process.
inline double cputime() {
  return static_cast<double>(std::clock()) / static_cast<double>(CLOCKS_PER_SEC);
}

// Peak resident set size in bytes.
inline double peakrss() {
#ifdef _WIN32
  PROCESS_MEMORY_COUNTERS info;
  if (GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info))) {
    return static_cast<double>(info.PeakWorkingSetSize);
  }
  return 0.0;
#else
  struct rusage usage {};
  if (getrusage(RUSAGE_SELF, &usage) == 0) {
#ifdef __APPLE__
    // macOS reports in bytes.
    return static_cast<double>(usage.ru_maxrss);
#else
    // Linux reports in kilobytes.
    return static_cast<double>(usage.ru_maxrss) * 1024.0;
#endif
  }
  return 0.0;
#endif
}

}  // namespace piru
