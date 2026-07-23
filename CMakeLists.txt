cmake_minimum_required(VERSION 3.16)

project(panomap VERSION 0.0.2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

# Optimization/debug defaults: RelWithDebInfo = O3 + debug info + frame pointers; Release = O3 only.
if(MSVC)
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/O2 /Zi /Oy-" CACHE STRING "RelWithDebInfo flags" FORCE)
    set(CMAKE_CXX_FLAGS_RELEASE "/O2" CACHE STRING "Release flags" FORCE)
else()
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -fno-omit-frame-pointer" CACHE STRING "RelWithDebInfo flags" FORCE)
    set(CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE STRING "Release flags" FORCE)
endif()

option(PANOMAP_ENABLE_WARNINGS "Enable compiler warnings" ON)
option(PANOMAP_FETCH_TBB "Fetch oneTBB via FetchContent if not found" ON)
option(PANOMAP_FETCH_DOCTEST "Fetch doctest during configure (for tests)" ON)
option(PANOMAP_ENABLE_BASE_MODE "Build experimental base-mode binary (panomap-base)" OFF)
# zstd-compressed BLOW5 support. AUTO = enable if libzstd is found, otherwise
# build without it. Force with ON (errors if libzstd missing) or OFF.
set(PANOMAP_USE_ZSTD "AUTO" CACHE STRING "zstd-compressed BLOW5 support: ON, OFF, or AUTO")
find_package(ZLIB REQUIRED)

set(PANOMAP_ZSTD_ENABLED OFF)
if(NOT PANOMAP_USE_ZSTD STREQUAL "OFF")
    find_library(ZSTD_LIBRARY zstd)
    find_path(ZSTD_INCLUDE_DIR zstd.h)
    if(ZSTD_LIBRARY AND ZSTD_INCLUDE_DIR)
        set(PANOMAP_ZSTD_ENABLED ON)
        message(STATUS "zstd found (${ZSTD_LIBRARY}); enabling zstd-compressed BLOW5 support.")
    elseif(PANOMAP_USE_ZSTD STREQUAL "ON")
        message(FATAL_ERROR "PANOMAP_USE_ZSTD=ON but libzstd not found. Install libzstd-dev (or set PANOMAP_USE_ZSTD=OFF).")
    else()
        message(STATUS "libzstd not found; zstd-compressed BLOW5 files will not be readable. Install libzstd-dev to enable.")
    endif()
endif()

set(PANOMAP_GENERATED_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
file(MAKE_DIRECTORY ${PANOMAP_GENERATED_DIR})

set(PANOMAP_VERSION_STRING "${PROJECT_NAME} ${PROJECT_VERSION} (experimental scaffold)")
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/version.hpp.in
    ${PANOMAP_GENERATED_DIR}/version.hpp
    @ONLY
)

# oneTBB is required (parallel index/map backend). Use the system install if
# present, otherwise auto-fetch a pinned copy (PANOMAP_FETCH_TBB, on by default).
find_package(TBB QUIET)
if(NOT TBB_FOUND AND PANOMAP_FETCH_TBB)
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 4.0 AND NOT CMAKE_POLICY_VERSION_MINIMUM)
        message(FATAL_ERROR
            "Fetching pinned oneTBB v2021.10.0, whose CMakeLists needs a "
            "CMake <3.5 policy and errors under CMake ${CMAKE_VERSION}.\n"
            "Re-run:  cmake .. -DCMAKE_POLICY_VERSION_MINIMUM=3.5")
    endif()
    include(FetchContent)
    FetchContent_Declare(
        TBB
        GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
        GIT_TAG v2021.10.0
    )
    set(TBB_TEST OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(TBB)
endif()
if(NOT (TBB_FOUND OR TARGET TBB::tbb))
    message(FATAL_ERROR "oneTBB is required but was not found. Install it, or leave PANOMAP_FETCH_TBB=ON to fetch it automatically.")
endif()

# slow5lib (signal-mode read I/O). Pulled in early so panomap_signal can link against it.
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third_party/slow5lib/CMakeLists.txt)
    set(SLOW5_LINK_STATIC ON CACHE BOOL "" FORCE)
    if(PANOMAP_ZSTD_ENABLED)
        # Compiles slow5_press.c with -DSLOW5_USE_ZSTD (guarded code paths).
        set(SLOW5_USE_ZSTD ON CACHE BOOL "" FORCE)
        # slow5lib includes <zstd.h> and links a bare -lzstd, but adds no search
        # paths for either; supply them here for prefixes outside the default
        # (e.g. Homebrew /opt/homebrew, conda, /usr/local).
        include_directories(${ZSTD_INCLUDE_DIR})
        get_filename_component(PANOMAP_ZSTD_LIBDIR ${ZSTD_LIBRARY} DIRECTORY)
        link_directories(${PANOMAP_ZSTD_LIBDIR})
    endif()
    add_subdirectory(third_party/slow5lib)
else()
    message(FATAL_ERROR "slow5lib submodule not present; required for signal-mode read I/O.")
endif()

# ----- panomap_core: modality-agnostic static library -----
add_library(panomap_core STATIC
    src/core/mapping/path_chainer.cpp
    src/core/mapping/sort_chainer.cpp
    src/core/mapping/pan_chainer.cpp
    src/core/mapping/seed_merger.cpp
    src/core/index/flat_graph.cpp
    src/core/index/simple_expand.cpp
    src/core/index/sort_1d.cpp
    src/core/concurrency/executor.cpp
    src/core/util/signal_handlers.cpp
    src/core/util/system_memory.cpp
    src/core/util/timing.cpp
    src/core/util/trace.cpp
    src/core/io/graphs/gfa_exporter.cpp
    src/core/io/graphs/gfa_loader.cpp
    src/core/io/graphs/graph_loader_factory.cpp
    src/core/io/index/serialization.cpp
    src/core/io/results/gaf_writer.cpp
    src/core/io/results/paf_writer.cpp
    src/core/io/results/result_writer_factory.cpp
)
target_include_directories(panomap_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${PANOMAP_GENERATED_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/dirtyzipf
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/concurrentqueue
)
target_link_libraries(panomap_core PUBLIC TBB::tbb)
target_compile_definitions(panomap_core PUBLIC PANOMAP_USE_TBB)

# pthread for std::async / std::thread used in sort_1d.cpp parallel SGD
find_package(Threads REQUIRED)
target_link_libraries(panomap_core PUBLIC Threads::Threads)
if(PANOMAP_ENABLE_WARNINGS)
    if(MSVC)
        target_compile_options(panomap_core PRIVATE /W4 /permissive-)
    else()
        target_compile_options(panomap_core PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endif()

# ----- panomap_signal: signal-mode static library -----
add_library(panomap_signal STATIC
    src/signal/mapping/batch_mapper.cpp
    src/signal/index/bucket_indexer.cpp
    src/signal/index/index_pipeline.cpp
    src/signal/diff_filter.cpp
    src/signal/event_pipeline_factory.cpp
    src/signal/event_pipelines/standard_event_pipeline.cpp
    src/signal/tokenizers/tokenizer_factory.cpp
    src/signal/tokenizers/rh2_tokenizer.cpp
    src/signal/tokenizers/landmark_tokenizer.cpp
    src/signal/seed_extractors/seed_extractor_factory.cpp
    src/signal/seed_extractors/kmer_seed_extractor.cpp
    src/signal/seed_extractors/minimizer_seed_extractor.cpp
    src/signal/io/reads/pod5_provider.cpp
    src/signal/io/reads/read_provider_factory.cpp
    src/signal/io/reads/slow5_provider.cpp
    src/signal/io/models/builtin_models.cpp
    src/signal/io/models/model.cpp
    src/signal/commands/index.cpp
    src/signal/commands/inspect.cpp
    src/signal/commands/map.cpp
)
target_include_directories(panomap_signal
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
           ${CMAKE_CURRENT_SOURCE_DIR}/third_party/slow5lib/include
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(panomap_signal
    PUBLIC panomap_core slow5 ZLIB::ZLIB
)
target_compile_definitions(panomap_signal PUBLIC PANOMAP_HAS_SLOW5)
if(PANOMAP_ZSTD_ENABLED)
    # slow5 is static; -lzstd must appear on panomap_signal's final link line.
    target_link_libraries(panomap_signal PUBLIC ${ZSTD_LIBRARY})
endif()
if(PANOMAP_ENABLE_WARNINGS)
    if(MSVC)
        target_compile_options(panomap_signal PRIVATE /W4 /permissive-)
    else()
        target_compile_options(panomap_signal PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endif()

# ----- panomap executable -----
add_executable(panomap
    src/cli/signal_main.cpp
    src/cli/parse.cpp
)
target_link_libraries(panomap PRIVATE panomap_signal)
if(PANOMAP_ENABLE_WARNINGS)
    if(MSVC)
        target_compile_options(panomap PRIVATE /W4 /permissive-)
    else()
        target_compile_options(panomap PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endif()

# ----- panomap-base: experimental base-mode binary (opt-in) -----
# Off by default. Enable with -DPANOMAP_ENABLE_BASE_MODE=ON. Paper 1
# uses signal mode only; base mode is kept as future/dev material
# (see tickets/archive/dev-108-base-mode.md).
if(PANOMAP_ENABLE_BASE_MODE)
    add_library(panomap_base STATIC
        src/base/index/base_indexer.cpp
        src/base/io/reads/fastq_provider.cpp
        src/base/seeder/base_seeder.cpp
        src/base/mapping/base_mapper.cpp
        src/base/commands/index.cpp
        src/base/commands/inspect.cpp
        src/base/commands/map.cpp
    )
    target_include_directories(panomap_base
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
    )
    target_link_libraries(panomap_base PUBLIC panomap_core ZLIB::ZLIB)
    if(PANOMAP_ENABLE_WARNINGS)
        if(MSVC)
            target_compile_options(panomap_base PRIVATE /W4 /permissive-)
        else()
            target_compile_options(panomap_base PRIVATE -Wall -Wextra -Wpedantic)
        endif()
    endif()

    add_executable(panomap-base
        src/cli/base_main.cpp
        src/cli/parse.cpp
    )
    target_link_libraries(panomap-base PRIVATE panomap_base)
    if(PANOMAP_ENABLE_WARNINGS)
        if(MSVC)
            target_compile_options(panomap-base PRIVATE /W4 /permissive-)
        else()
            target_compile_options(panomap-base PRIVATE -Wall -Wextra -Wpedantic)
        endif()
    endif()
endif()

include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

add_subdirectory(tools/panolayout)
