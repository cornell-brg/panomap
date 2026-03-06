#include "io/results/alignment_conversion.hpp"

#ifdef PIRU_HAS_LIBVGIO
#include <vg/vg.pb.h>
#endif

namespace piru::io {

#ifdef PIRU_HAS_LIBVGIO

vg::Alignment to_vg_alignment(const AlignmentResult& result) {
  vg::Alignment aln;
  aln.set_name(result.query_name);
  aln.set_sequence(result.query_sequence);
  if (!result.query_quality.empty()) {
    aln.set_quality(result.query_quality);
  }
  aln.set_mapping_quality(result.mapq);

  auto* path = aln.mutable_path();
  for (const auto& m : result.mappings) {
    auto* mapping = path->add_mapping();
    mapping->mutable_position()->set_node_id(static_cast<int64_t>(m.node_id));
    mapping->mutable_position()->set_offset(m.offset);
    mapping->mutable_position()->set_is_reverse(m.is_reverse);
    for (const auto& e : m.edits) {
      auto* edit = mapping->add_edit();
      edit->set_from_length(e.from_length);
      edit->set_to_length(e.to_length);
      if (!e.sequence.empty()) {
        edit->set_sequence(e.sequence);
      }
    }
  }
  return aln;
}

#endif

}  // namespace piru::io
