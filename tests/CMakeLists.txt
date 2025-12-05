# doctest setup and basic test target

include(FetchContent)
set(DOCTEST_AVAILABLE OFF)
set(DOCTEST_DISCOVERY_CMAKE "")

if(TARGET doctest::doctest)
    set(DOCTEST_AVAILABLE ON)
else()
    find_package(doctest CONFIG QUIET)
    if(doctest_FOUND)
        set(DOCTEST_AVAILABLE ON)
        if(EXISTS "${doctest_DIR}/doctest.cmake")
            set(DOCTEST_DISCOVERY_CMAKE "${doctest_DIR}/doctest.cmake")
        endif()
    endif()
endif()

if(NOT DOCTEST_AVAILABLE AND PIRU_FETCH_DOCTEST)
    FetchContent_Declare(
        doctest
        GIT_REPOSITORY https://github.com/doctest/doctest.git
        GIT_TAG v2.4.12
    )
    FetchContent_MakeAvailable(doctest)
    set(DOCTEST_AVAILABLE ON)
    set(DOCTEST_DISCOVERY_CMAKE "${doctest_SOURCE_DIR}/scripts/cmake/doctest.cmake")
    if(NOT TARGET doctest::doctest)
        add_library(doctest::doctest INTERFACE IMPORTED)
        target_include_directories(doctest::doctest INTERFACE ${doctest_SOURCE_DIR})
    endif()
endif()

if(DOCTEST_AVAILABLE)
    if(DOCTEST_DISCOVERY_CMAKE)
        include(${DOCTEST_DISCOVERY_CMAKE})
    endif()

    add_executable(piru_tests
        test_main.cpp
        test_metrics.cpp
        test_models.cpp
        test_graphs.cpp
        test_results.cpp
        test_reads.cpp
        test_signal_mocks.cpp
        test_signal_events.cpp
        test_signal_normalizers.cpp
        test_signal_fuzzy_quantize.cpp
        test_signal_alignment_quantize.cpp
        test_signal_seed_extractors.cpp
    )

    target_include_directories(piru_tests PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_BINARY_DIR}/generated
    )

    target_sources(piru_tests PRIVATE
        ${CMAKE_SOURCE_DIR}/src/io/graphs/graph_loader_factory.cpp
        ${CMAKE_SOURCE_DIR}/src/io/graphs/gfa_loader.cpp
        ${CMAKE_SOURCE_DIR}/src/io/graphs/vg_loader.cpp
        ${CMAKE_SOURCE_DIR}/src/io/results/gaf_writer.cpp
        ${CMAKE_SOURCE_DIR}/src/io/results/gam_writer.cpp
        ${CMAKE_SOURCE_DIR}/src/io/results/json_writer.cpp
        ${CMAKE_SOURCE_DIR}/src/io/results/alignment_conversion.cpp
        ${CMAKE_SOURCE_DIR}/src/io/results/result_writer_factory.cpp
        ${CMAKE_SOURCE_DIR}/src/io/models/builtin_models.cpp
        ${CMAKE_SOURCE_DIR}/src/signal/event_detectors/event_detector_factory.cpp
        ${CMAKE_SOURCE_DIR}/src/signal/event_detectors/passthrough_event_detector.cpp
        ${CMAKE_SOURCE_DIR}/src/signal/event_detectors/scrappie_event_detector.cpp
        ${CMAKE_SOURCE_DIR}/src/signal/normalizers/identity_normalizer.cpp
        ${CMAKE_SOURCE_DIR}/src/signal/normalizers/normalizer_factory.cpp
        ${CMAKE_SOURCE_DIR}/src/signal/normalizers/zscore_normalizer.cpp
        ${CMAKE_SOURCE_DIR}/src/signal/normalizers/median_mad_normalizer.cpp
        ${CMAKE_SOURCE_DIR}/src/signal/fuzzy_quantizers/fuzzy_quantizer_factory.cpp
        ${CMAKE_SOURCE_DIR}/src/signal/fuzzy_quantizers/passthrough_fuzzy_quantizer.cpp
        ${CMAKE_SOURCE_DIR}/src/signal/fuzzy_quantizers/rh2_fuzzy_quantizer.cpp
        ${CMAKE_SOURCE_DIR}/src/signal/alignment_quantizers/alignment_quantizer_factory.cpp
        ${CMAKE_SOURCE_DIR}/src/signal/alignment_quantizers/passthrough_alignment_quantizer.cpp
        ${CMAKE_SOURCE_DIR}/src/signal/alignment_quantizers/int_alignment_quantizer.cpp
        ${CMAKE_SOURCE_DIR}/src/signal/seed_extractors/seed_extractor_factory.cpp
        ${CMAKE_SOURCE_DIR}/src/signal/seed_extractors/kmer_seed_extractor.cpp
    )

    target_link_libraries(piru_tests PRIVATE doctest::doctest)
    if(PIRU_LIBVGIO_AVAILABLE)
        target_link_libraries(piru_tests PRIVATE VGio::vgio)
        target_compile_definitions(piru_tests PRIVATE PIRU_HAS_LIBVGIO)
    endif()

    if(DOCTEST_DISCOVERY_CMAKE)
        doctest_discover_tests(piru_tests
            TEST_PREFIX unit.
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        )
    else()
        add_test(NAME piru_tests COMMAND piru_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    endif()

    add_custom_target(check
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
        DEPENDS piru_tests
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
else()
    message(STATUS "doctest not found (set PIRU_FETCH_DOCTEST=ON to fetch); skipping tests.")
endif()
