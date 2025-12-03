# Read I/O

PIRU streams raw reads through a pluggable `ReadProvider` interface. Each backend implements `get_next`, `reset`, and `get_format_name` to hide format details from the mapper.

## Formats
- **slow5/blow5**: Primary backend via slow5lib (required dependency). Use `make_read_provider(path)` to get a provider; `map` subcommand can list reads from Blow5 test data.

## On the road-map
- **pod5**: Planned. Stub exists; will be implemented once a stable C++ reader is selected.
- **read until API**: no stubs yet, but planned to be added to support real-time mapping.

## Not supported
- **fast5**: Not supported. Users should convert fast5 to slow5/blow5 (e.g., with slow5tools) before running PIRU.

## Interface location
- Interface: `include/io/reads/read_provider.hpp`
- Factory: `include/io/reads/read_provider_factory.hpp`
- slow5 backend: `src/io/reads/slow5_provider.cpp`

## Notes
- slow5lib is required; CMake will fail if the submodule is missing.
- Extension-based selection is used in the factory; explicit selection can be added later if needed.
