# Roadmap (PIRU workspace)

A living list of follow-on tasks mentioned in DEV plans and code comments. This is informational only; not a commitment or schedule.

## Reads I/O
- Add POD5 backend (investigate lib/pipeline and conversion options).
- Keep slow5lib updated; verify C++17-friendly headers or upstream patches.

## Models
- None pending in DEV002 (built-ins + file loader are done).

## Graph I/O
- Add GBZ/GBWTGraph loader behind the graph factory.
- Expand GFA support (walks/containments, optional fields, stricter validation).
- Consider richer vg path metadata and optional edge/tag parsing.

## Results I/O
- Add richer tags/metadata to `AlignmentResult` and writers as needed by the mapper.
- Add optional debug/JSON schemas or pretty-print helpers.
- Integrate writers into mapping pipeline (select format from CLI).

## Alignment Graph (future)
- Define internal `AlnGraph` and conversion from `ImportedGraph`; clarify how overlaps are used.
- Provide optional GFA export for debugging internal graph state.

## Tooling/CI
- Consider adding tests around GAM output round-trips (vg view) when alignment data is richer.
- Keep third-party deps (libvgio, slow5lib) in sync and documented.
