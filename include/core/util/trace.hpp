/**
 * trace.hpp
 *
 * Compile-time optional pipeline tracing for debugging.
 * Build with -DPIRU_TRACE_ENABLED to enable. Zero overhead when compiled out.
 *
 * At runtime, tracing is OFF by default even when compiled in.
 * Set PIRU_TRACE_DIR to activate:
 *   PIRU_TRACE_DIR=/tmp/trace          Output directory (required to enable)
 *   PIRU_TRACE_STAGES=0x3f             Bitmask of stages to dump (default: all)
 *   PIRU_TRACE_READS=MSH2,NF2          Comma-separated read ID substrings (default: all)
 *
 * Usage in code:
 *   PIRU_TRACE_DUMP(kTokens, read_id, {
 *     std::ofstream ofs(trace_path("tokens", read_id, chunk_idx));
 *     for (auto t : tokens) ofs << t << "\n";
 *   });
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <cstdint>
#include <string>

namespace piru::trace {

enum Stage : std::uint32_t {
  kSignal = 1 << 0,  // Raw pA signal
  kNorm = 1 << 1,    // Normalized signal (post z-score)
  kEvents = 1 << 2,  // Detected events
  kTokens = 1 << 3,  // Quantized tokens (post diff filter)
  kSeeds = 1 << 4,   // Extracted seeds
  kHits = 1 << 5,    // Seed lookup hits
  kChains = 1 << 6,  // Chain scores / weighted decision
  kAll = 0x7f,
};

// Always declared (so code compiles even when tracing is off)
std::uint32_t enabled_stages();
bool match_read(const std::string& read_id);
std::string trace_path(const std::string& tag, const std::string& read_id,
                       std::size_t chunk_idx = 0);

#ifdef PIRU_TRACE_ENABLED

#define PIRU_TRACE_DUMP(stage, read_id, block)                                               \
  do {                                                                                       \
    if ((::piru::trace::enabled_stages() & (stage)) && ::piru::trace::match_read(read_id)) { \
      block                                                                                  \
    }                                                                                        \
  } while (0)

#else

// When tracing is disabled, compile the block but never execute it.
// This ensures the code stays valid (catches typos/refactoring breaks).
#define PIRU_TRACE_DUMP(stage, read_id, block) \
  do {                                         \
    if (false) {                               \
      (void)(read_id);                         \
      block                                    \
    }                                          \
  } while (0)

#endif

}  // namespace piru::trace
