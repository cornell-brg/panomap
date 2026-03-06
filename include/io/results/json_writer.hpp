// SPDX-License-Identifier: MIT
// Writer for JSONL-formatted alignment results (debug-friendly).

#pragma once

#include <string>

#include "io/results/result_writer.hpp"

#ifdef PIRU_HAS_LIBVGIO
#include <vg/io/alignment_emitter.hpp>
#endif

namespace piru::io {

class JsonWriter : public ResultWriter {
public:
  explicit JsonWriter(const std::string& path);
  ~JsonWriter() override;

  bool write(const AlignmentResult& result) override;

private:
#ifdef PIRU_HAS_LIBVGIO
  std::unique_ptr<vg::io::AlignmentEmitter> emitter_;
#endif
};

using JsonWriterPtr = std::unique_ptr<JsonWriter>;

}  // namespace piru::io
