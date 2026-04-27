/**
 * fastq_provider.cpp
 *
 * gzopen-based FASTQ reader. gzopen handles both plain and gzip files
 * transparently, so we don't need to branch on the filename suffix.
 *
 * Record layout (one record == 4 lines):
 *   @<id> [optional comment]
 *   <sequence>
 *   +
 *   <quality>
 *
 * SPDX-License-Identifier: MIT
 */

#include "base/io/reads/fastq_provider.hpp"

#include <cstring>
#include <zlib.h>

#include "core/util/logging.hpp"

namespace piru::base::io {

namespace {

constexpr std::size_t kInitialLineCap = 4096;

// Read one line from gzFile into `buf`, stripping trailing CR/LF. Returns
// false on EOF (with empty buf) or unrecoverable error.
bool read_line(gzFile gz, std::vector<char>& buf) {
  buf.clear();
  if (buf.capacity() < kInitialLineCap) buf.reserve(kInitialLineCap);

  while (true) {
    // Reserve room for at least 256 more chars; gzgets needs a writable
    // pointer with a known capacity.
    std::size_t old_size = buf.size();
    std::size_t want = std::max<std::size_t>(buf.capacity(), old_size + 256);
    if (buf.capacity() < want) buf.reserve(want);
    buf.resize(buf.capacity());

    char* dst = buf.data() + old_size;
    int chunk = static_cast<int>(buf.size() - old_size);
    if (chunk <= 1) {
      buf.reserve(buf.capacity() * 2);
      buf.resize(buf.capacity());
      dst = buf.data() + old_size;
      chunk = static_cast<int>(buf.size() - old_size);
    }

    if (gzgets(gz, dst, chunk) == nullptr) {
      // EOF or error.
      buf.resize(old_size);
      return old_size > 0;  // partial line at EOF is still data
    }

    std::size_t got = std::strlen(dst);
    buf.resize(old_size + got);

    if (!buf.empty() && buf.back() == '\n') {
      buf.pop_back();
      if (!buf.empty() && buf.back() == '\r') buf.pop_back();
      return true;
    }
    // No newline yet; loop to read more (long line).
  }
}

}  // namespace

FastqProvider::FastqProvider(const std::string& path) : path_(path) {
  file_ = gzopen(path.c_str(), "rb");
  if (!file_) {
    LOG_ERROR("fastq_provider: failed to open " + path);
  }
  line_buf_.reserve(kInitialLineCap);
}

FastqProvider::~FastqProvider() {
  if (file_) gzclose(static_cast<gzFile>(file_));
}

bool FastqProvider::next(FastqRead& out) {
  if (!file_) return false;
  gzFile gz = static_cast<gzFile>(file_);

  // Header line. Accept both '@' (FASTQ) and '>' (single-line FASTA);
  // the workspace eval data ships some "ground truth" sequences with
  // .fastq.gz extensions but FASTA contents (squigulator GT outputs).
  if (!read_line(gz, line_buf_)) return false;
  while (line_buf_.empty()) {
    if (!read_line(gz, line_buf_)) return false;
  }
  const char header_char = line_buf_.front();
  if (header_char != '@' && header_char != '>') {
    if (!warned_malformed_) {
      LOG_WARN("fastq_provider: expected '@' or '>' header in " + path_ + ", skipping rest");
      warned_malformed_ = true;
    }
    return false;
  }
  // Strip leading sigil, trim id at first whitespace.
  std::size_t end = 1;
  while (end < line_buf_.size() && line_buf_[end] != ' ' && line_buf_[end] != '\t') ++end;
  out.id.assign(line_buf_.data() + 1, end - 1);
  out.source_file = path_;

  // Sequence line
  if (!read_line(gz, line_buf_)) return false;
  out.sequence.assign(line_buf_.data(), line_buf_.size());

  if (header_char == '@') {
    // FASTQ: '+' separator then quality line (read and discard).
    if (!read_line(gz, line_buf_)) return false;
    if (line_buf_.empty() || line_buf_.front() != '+') {
      if (!warned_malformed_) {
        LOG_WARN("fastq_provider: expected '+' separator in " + path_ + ", skipping rest");
        warned_malformed_ = true;
      }
      return false;
    }
    if (!read_line(gz, line_buf_)) return false;
  }

  return true;
}

}  // namespace piru::base::io
