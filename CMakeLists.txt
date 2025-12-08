cmake_minimum_required(VERSION 3.16)

project(piru VERSION 0.0.1 LANGUAGES CXX)

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
option(PIRU_ENABLE_LIBVGIO "Enable vg loading via libvgio" ON)
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
        FetchContent_MakeAvailable(TBB)
    endif()
    if(TBB_FOUND OR TARGET TBB::tbb)
        set(PIRU_TBB_AVAILABLE ON)
    else()
        message(FATAL_ERROR "oneTBB requested but not found. Set PIRU_FETCH_TBB=ON to fetch automatically.")
    endif()
endif()

add_executable(piru
    src/main.cpp
    src/commands/index.cpp
    src/commands/map.cpp
    src/commands/eval.cpp
    src/commands/mt_test.cpp
    src/concurrency/executor.cpp
    src/io/graphs/graph_loader_factory.cpp
    src/io/graphs/gfa_loader.cpp
    src/io/graphs/vg_loader.cpp
    src/io/results/gaf_writer.cpp
    src/io/results/gam_writer.cpp
    src/io/results/json_writer.cpp
    src/io/results/alignment_conversion.cpp
    src/io/results/result_writer_factory.cpp
    src/io/reads/read_provider_factory.cpp
    src/io/reads/pod5_provider.cpp
    src/io/reads/slow5_provider.cpp
    src/io/models/builtin_models.cpp
    src/cli/parse.cpp
    src/signal/event_detectors/event_detector_factory.cpp
    src/signal/event_detectors/passthrough_event_detector.cpp
    src/signal/event_detectors/scrappie_event_detector.cpp
    src/signal/normalizers/normalizer_factory.cpp
    src/signal/normalizers/identity_normalizer.cpp
    src/signal/normalizers/zscore_normalizer.cpp
    src/signal/normalizers/median_mad_normalizer.cpp
    src/signal/fuzzy_quantizers/fuzzy_quantizer_factory.cpp
    src/signal/fuzzy_quantizers/passthrough_fuzzy_quantizer.cpp
    src/signal/fuzzy_quantizers/rh2_fuzzy_quantizer.cpp
    src/signal/alignment_quantizers/alignment_quantizer_factory.cpp
    src/signal/alignment_quantizers/passthrough_alignment_quantizer.cpp
    src/signal/alignment_quantizers/int_alignment_quantizer.cpp
    src/signal/seed_extractors/seed_extractor_factory.cpp
    src/signal/seed_extractors/kmer_seed_extractor.cpp
    src/index/aln_graph.cpp
    src/index/transform_dbg.cpp
    src/index/pseudo_linearize.cpp
    src/index/squigglize.cpp
    src/index/seed_builder.cpp
    src/mapping/batch_mapper.cpp
    src/util/signal_handlers.cpp
    src/util/timing.cpp
)

target_include_directories(piru PRIVATE include ${CMAKE_CURRENT_SOURCE_DIR}/src ${PIRU_GENERATED_DIR})

if(PIRU_ENABLE_WARNINGS)
    if(MSVC)
        target_compile_options(piru PRIVATE /W4 /permissive-)
    else()
        target_compile_options(piru PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endif()

if(PIRU_TBB_AVAILABLE)
    target_link_libraries(piru PRIVATE TBB::tbb)
    target_compile_definitions(piru PRIVATE PIRU_USE_TBB)
else()
    message(STATUS "oneTBB not found; using serial concurrency backend.")
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third_party/slow5lib/CMakeLists.txt)
    add_subdirectory(third_party/slow5lib)
    target_link_libraries(piru PRIVATE slow5 ZLIB::ZLIB)
    target_compile_definitions(piru PRIVATE PIRU_HAS_SLOW5)
    target_include_directories(piru PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/slow5lib/include)
else()
    message(FATAL_ERROR "slow5lib submodule not present; required for read I/O.")
endif()

set(PIRU_LIBVGIO_AVAILABLE OFF)
if(PIRU_ENABLE_LIBVGIO)
    set(PIRU_LIBVGIO_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party/libvgio")
    if(EXISTS ${PIRU_LIBVGIO_DIR}/CMakeLists.txt)
        set(PIRU_LIBVGIO_CAN_BUILD ON)
        find_package(PkgConfig QUIET)
        if(PKG_CONFIG_FOUND)
            pkg_check_modules(PIRU_JANSSON QUIET jansson)
            pkg_check_modules(PIRU_HTSLIB QUIET htslib)
            if(NOT PIRU_JANSSON_FOUND OR NOT PIRU_HTSLIB_FOUND)
                set(PIRU_LIBVGIO_CAN_BUILD OFF)
                message(WARNING "libvgio dependencies missing (jansson/htslib); set PIRU_ENABLE_LIBVGIO=OFF or install dependencies to enable vg loading.")
            endif()
        endif()

        if(PIRU_LIBVGIO_CAN_BUILD)
            add_subdirectory(${PIRU_LIBVGIO_DIR} libvgio)
            target_link_libraries(piru PRIVATE VGio::vgio)
            target_compile_definitions(piru PRIVATE PIRU_HAS_LIBVGIO)
            set(PIRU_LIBVGIO_AVAILABLE ON)
        endif()
    else()
        message(FATAL_ERROR "libvgio submodule not present; required for vg loading.")
    endif()
else()
    message(STATUS "PIRU_ENABLE_LIBVGIO=OFF; vg loading disabled.")
endif()

include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
