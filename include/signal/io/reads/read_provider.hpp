// SPDX-License-Identifier: MIT
// Abstraction for streaming raw reads from multiple formats.

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace panomap::io {

struct RawRead {
  std::string read_id;
  std::string source_file;  // Input file this read came from (for per-file output routing)
  std::vector<int16_t> raw_signal;
  std::uint64_t len_raw_signal{0};
  float range{0.0f};
  float digitisation{0.0f};
  float offset{0.0f};
  float sampling_rate_hz{0.0f};
};

class ReadProvider {
public:
  virtual ~ReadProvider() = default;

  // Populate `read` with the next entry. Returns false at EOF.
  virtual bool get_next(RawRead& read) = 0;

  // Reset the provider (used primarily for testing/replays).
  virtual void reset() = 0;

  // Human-readable format name (e.g., slow5, pod5, fast5).
  virtual std::string get_format_name() const = 0;
};

using ReadProviderPtr = std::unique_ptr<ReadProvider>;

}  // namespace panomap::io
