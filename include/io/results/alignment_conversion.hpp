// SPDX-License-Identifier: MIT
// Helpers to convert AlignmentResult to vg::Alignment for libvgio emitters.

#pragma once

#ifdef PIRU_HAS_LIBVGIO
#include <vg/vg.pb.h>
#endif

#include "io/results/result.hpp"

namespace piru::io {

#ifdef PIRU_HAS_LIBVGIO
vg::Alignment to_vg_alignment(const AlignmentResult& result);
#endif

}  // namespace piru::io
