#include "io/results/gam_writer.hpp"

#include "io/results/alignment_conversion.hpp"
#include "util/logging.hpp"

namespace piru::io {

GamWriter::GamWriter(const std::string& path) : path_(path) {
#ifdef PIRU_HAS_LIBVGIO
  emitter_ =
      vg::io::get_non_hts_alignment_emitter(path_, "GAM", {}, /*threads=*/1, nullptr, nullptr);
#else
  (void)path_;
#endif
}

bool GamWriter::write(const AlignmentResult& result) {
#ifdef PIRU_HAS_LIBVGIO
  if (!emitter_) {
    if (!warned_) {
      LOG_ERROR("GAM writer not initialized for '" + path_ + "'");
      warned_ = true;
    }
    return false;
  }
  vg::Alignment aln = to_vg_alignment(result);
  emitter_->emit_single(std::move(aln));
  return true;
#else
  if (!warned_) {
    LOG_ERROR("libvgio not linked; GAM writing unavailable for '" + path_ + "'");
    warned_ = true;
  }
  return false;
#endif
}

}  // namespace piru::io
