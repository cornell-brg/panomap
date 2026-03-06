#include "io/results/json_writer.hpp"

#include "io/results/alignment_conversion.hpp"
#include "util/logging.hpp"

namespace piru::io {

JsonWriter::JsonWriter(const std::string& path) {
#ifdef PIRU_HAS_LIBVGIO
  emitter_ =
      vg::io::get_non_hts_alignment_emitter(path, "JSON", {}, /*threads=*/1, nullptr, nullptr);
#else
  (void)path;
#endif
}

JsonWriter::~JsonWriter() = default;

bool JsonWriter::write(const AlignmentResult& result) {
#ifdef PIRU_HAS_LIBVGIO
  if (!emitter_) {
    LOG_ERROR("JSON writer not initialized");
    return false;
  }
  emitter_->emit_single(to_vg_alignment(result));
  return true;
#else
  LOG_ERROR("libvgio not linked; JSON writing unavailable");
  return false;
#endif
}

}  // namespace piru::io
