cmake_minimum_required(VERSION 3.16)

project(piru VERSION 0.0.2 LANGUAGES CXX)

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

option(PIRU_ENABLE_WARNINGS "Enable compiler warnings" ON)
option(PIRU_USE_TBB "Enable oneTBB backend for concurrency" ON)
option(PIRU_FETCH_TBB "Fetch oneTBB via FetchContent if not found" ON)
option(PIRU_FETCH_DOCTEST "Fetch doctest during configure (for tests)" ON)
find_package(ZLIB REQUIRED)

set(PIRU_GENERATED_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
file(MAKE_DIRECTORY ${PIRU_GENERATED_DIR})

set(PIRU_VERSION_STRING "${PROJECT_NAME} ${PROJECT_VERSION} (experimental scaffold)")
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/version.hpp.in
    ${PIRU_GENERATED_DIR}/version.hpp
    @ONLY
)

set(PIRU_TBB_AVAILABLE OFF)
if(PIRU_USE_TBB)
    find_package(TBB QUIET)
    if(NOT TBB_FOUND AND PIRU_FETCH_TBB)
        include(FetchContent)
        FetchContent_Declare(
            TBB
            GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
            GIT_TAG v2021.10.0
        )
        set(TBB_TEST OFF CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(TBB)
    endif()
    if(TBB_FOUND OR TARGET TBB::tbb)
        set(PIRU_TBB_AVAILABLE ON)
    else()
        message(FATAL_ERROR "oneTBB requested but not found. Set PIRU_FETCH_TBB=ON to fetch automatically.")
    endif()
endif()

# slow5lib (signal-mode read I/O). Pulled in early so piru_signal can link against it.
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third_party/slow5lib/CMakeLists.txt)
    add_subdirectory(third_party/slow5lib)
else()
    message(FATAL_ERROR "slow5lib submodule not present; required for signal-mode read I/O.")
endif()

# ----- piru_core: modality-agnostic static library -----
add_library(piru_core STATIC
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
target_include_directories(piru_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${PIRU_GENERATED_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/dirtyzipf
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/concurrentqueue
)
if(PIRU_TBB_AVAILABLE)
    target_link_libraries(piru_core PUBLIC TBB::tbb)
    target_compile_definitions(piru_core PUBLIC PIRU_USE_TBB)
else()
    message(STATUS "oneTBB not found; using serial concurrency backend.")
endif()

# pthread for std::async / std::thread used in sort_1d.cpp parallel SGD
find_package(Threads REQUIRED)
target_link_libraries(piru_core PUBLIC Threads::Threads)
if(PIRU_ENABLE_WARNINGS)
    if(MSVC)
        target_compile_options(piru_core PRIVATE /W4 /permissive-)
    else()
        target_compile_options(piru_core PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endif()

# ----- piru_signal: signal-mode static library -----
add_library(piru_signal STATIC
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
target_include_directories(piru_signal
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
           ${CMAKE_CURRENT_SOURCE_DIR}/third_party/slow5lib/include
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(piru_signal
    PUBLIC piru_core slow5 ZLIB::ZLIB
)
target_compile_definitions(piru_signal PUBLIC PIRU_HAS_SLOW5)
if(PIRU_ENABLE_WARNINGS)
    if(MSVC)
        target_compile_options(piru_signal PRIVATE /W4 /permissive-)
    else()
        target_compile_options(piru_signal PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endif()

# ----- piru-signal executable -----
add_executable(piru-signal
    src/cli/signal_main.cpp
    src/cli/parse.cpp
)
target_link_libraries(piru-signal PRIVATE piru_signal)
if(PIRU_ENABLE_WARNINGS)
    if(MSVC)
        target_compile_options(piru-signal PRIVATE /W4 /permissive-)
    else()
        target_compile_options(piru-signal PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endif()

# ----- piru_base: base-mode static library -----
add_library(piru_base STATIC
    src/base/index/base_indexer.cpp
    src/base/io/reads/fastq_provider.cpp
    src/base/seeder/base_seeder.cpp
    src/base/mapping/base_mapper.cpp
    src/base/commands/index.cpp
    src/base/commands/inspect.cpp
    src/base/commands/map.cpp
)
target_include_directories(piru_base
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(piru_base PUBLIC piru_core ZLIB::ZLIB)
if(PIRU_ENABLE_WARNINGS)
    if(MSVC)
        target_compile_options(piru_base PRIVATE /W4 /permissive-)
    else()
        target_compile_options(piru_base PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endif()

# ----- piru-base executable -----
add_executable(piru-base
    src/cli/base_main.cpp
    src/cli/parse.cpp
)
target_link_libraries(piru-base PRIVATE piru_base)
if(PIRU_ENABLE_WARNINGS)
    if(MSVC)
        target_compile_options(piru-base PRIVATE /W4 /permissive-)
    else()
        target_compile_options(piru-base PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endif()

include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
