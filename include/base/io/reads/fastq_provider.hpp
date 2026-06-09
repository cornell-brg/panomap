/**
 * fastq_provider.hpp
 *
 * Streaming FASTQ reader for base-mode input. Plain text and gzip
 * (.fastq / .fastq.gz). Single-line records only -- the basecaller
 * outputs we use (Dorado, Guppy) all emit one-line sequences. Quality
 * lines are read and discarded.
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace panomap::base::io {

struct FastqRead {
  std::string id;          // header up to first whitespace, no leading '@'
  std::string sequence;    // bases as-read (preserve case; seeder normalises)
  std::string source_file; // path the read came from (for per-file output routing)
};

class FastqProvider {
 public:
  // Open `path`. Auto-detects gzip via filename suffix or magic bytes
  // (gzopen handles both transparently).
  explicit FastqProvider(const std::string& path);
  ~FastqProvider();

  FastqProvider(const FastqProvider&) = delete;
  FastqProvider& operator=(const FastqProvider&) = delete;

  // Fill `out` with the next record. Returns false on EOF or error.
  // Skips empty lines and warns once on malformed records.
  bool next(FastqRead& out);

  // True if open() succeeded.
  bool is_open() const { return file_ != nullptr; }

  const std::string& path() const { return path_; }

 private:
  std::string path_;
  void* file_{nullptr};       // gzFile, opaque to keep zlib out of the header
  std::vector<char> line_buf_;
  bool warned_malformed_{false};
};

}  // namespace panomap::base::io
