// SPDX-License-Identifier: MIT
// Stub GAM writer (placeholder until full GAM output is implemented).

#pragma once

#include "io/results/result_writer.hpp"

#ifdef PIRU_HAS_LIBVGIO
#include <vg/io/alignment_emitter.hpp>
#endif

namespace piru::io {

class GamWriter : public ResultWriter {
public:
    explicit GamWriter(const std::string& path);
    ~GamWriter() override = default;

    bool write(const AlignmentResult& result) override;

private:
    std::string path_;
#ifdef PIRU_HAS_LIBVGIO
    std::unique_ptr<vg::io::AlignmentEmitter> emitter_;
#endif
    bool warned_{false};
};

using GamWriterPtr = std::unique_ptr<GamWriter>;

}  // namespace piru::io
