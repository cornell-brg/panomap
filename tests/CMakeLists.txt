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
        test_aln_graph.cpp
        test_results.cpp
        test_reads.cpp
        test_signal_tokenize.cpp
        test_signal_seed_extractors.cpp
        test_mapping_batch.cpp
        test_index_pipeline.cpp
        test_cli.cpp
        test_logging.cpp
        test_gfa_exporter.cpp
        test_seed_merger.cpp
        test_path_chainer.cpp
        test_bucket_indexer.cpp
    )

    target_include_directories(piru_tests PRIVATE
        ${CMAKE_SOURCE_DIR}/src
    )

    target_link_libraries(piru_tests PRIVATE piru_signal doctest::doctest)

    if(DOCTEST_DISCOVERY_CMAKE)
        doctest_discover_tests(piru_tests
            TEST_PREFIX unit.
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        )
    else()
        add_test(NAME piru_tests COMMAND piru_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    endif()

    # Integration test: DRB1 pangenome mapping accuracy
    set(WORKSPACE_ROOT "${CMAKE_SOURCE_DIR}/..")
    add_test(
        NAME integration.drb1_accuracy
        COMMAND bash ${CMAKE_SOURCE_DIR}/tests/integration/test_drb1_accuracy.sh
                $<TARGET_FILE:piru-signal> ${WORKSPACE_ROOT}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    set_tests_properties(integration.drb1_accuracy PROPERTIES
        TIMEOUT 120
        LABELS "integration"
    )

    # Integration test: DRB1 with sort-chain backend
    add_test(
        NAME integration.drb1_sort_chain
        COMMAND bash ${CMAKE_SOURCE_DIR}/tests/integration/test_drb1_sort_chain.sh
                $<TARGET_FILE:piru-signal> ${WORKSPACE_ROOT}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    set_tests_properties(integration.drb1_sort_chain PROPERTIES
        TIMEOUT 120
        LABELS "integration"
    )

    # Integration test: DRB1 with pan-chain backend
    add_test(
        NAME integration.drb1_pan_chain
        COMMAND bash ${CMAKE_SOURCE_DIR}/tests/integration/test_drb1_pan_chain.sh
                $<TARGET_FILE:piru-signal> ${WORKSPACE_ROOT}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    set_tests_properties(integration.drb1_pan_chain PROPERTIES
        TIMEOUT 120
        LABELS "integration"
    )

    # Integration test: DRB1 with landmark tokenizer
    add_test(
        NAME integration.drb1_landmark
        COMMAND bash ${CMAKE_SOURCE_DIR}/tests/integration/test_drb1_landmark.sh
                $<TARGET_FILE:piru-signal> ${WORKSPACE_ROOT}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    set_tests_properties(integration.drb1_landmark PROPERTIES
        TIMEOUT 120
        LABELS "integration"
    )

    # Integration test: DRB1 with landmark tokenizer + sort-chain
    add_test(
        NAME integration.drb1_landmark_sort_chain
        COMMAND bash ${CMAKE_SOURCE_DIR}/tests/integration/test_drb1_landmark_sort_chain.sh
                $<TARGET_FILE:piru-signal> ${WORKSPACE_ROOT}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    set_tests_properties(integration.drb1_landmark_sort_chain PROPERTIES
        TIMEOUT 120
        LABELS "integration"
    )

    # Integration test: DRB1 with landmark tokenizer + pan-chain
    add_test(
        NAME integration.drb1_landmark_pan_chain
        COMMAND bash ${CMAKE_SOURCE_DIR}/tests/integration/test_drb1_landmark_pan_chain.sh
                $<TARGET_FILE:piru-signal> ${WORKSPACE_ROOT}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    set_tests_properties(integration.drb1_landmark_pan_chain PROPERTIES
        TIMEOUT 120
        LABELS "integration"
    )

    # Integration test: COVID linear reference (real reads, mm2 truth)
    add_test(
        NAME integration.covid_accuracy
        COMMAND bash ${CMAKE_SOURCE_DIR}/tests/integration/test_covid_accuracy.sh
                $<TARGET_FILE:piru-signal> ${WORKSPACE_ROOT}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    set_tests_properties(integration.covid_accuracy PROPERTIES
        TIMEOUT 120
        LABELS "integration"
    )

    # Integration test: COVID with RH2 viral preset
    add_test(
        NAME integration.covid_viral
        COMMAND bash ${CMAKE_SOURCE_DIR}/tests/integration/test_covid_viral.sh
                $<TARGET_FILE:piru-signal> ${WORKSPACE_ROOT}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    set_tests_properties(integration.covid_viral PROPERTIES
        TIMEOUT 120
        LABELS "integration"
    )

    # Integration test: COVID with landmark tokenizer (best config)
    add_test(
        NAME integration.covid_landmark_viral
        COMMAND bash ${CMAKE_SOURCE_DIR}/tests/integration/test_covid_landmark_viral.sh
                $<TARGET_FILE:piru-signal> ${WORKSPACE_ROOT}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    set_tests_properties(integration.covid_landmark_viral PROPERTIES
        TIMEOUT 120
        LABELS "integration"
    )

    add_custom_target(check
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
        DEPENDS piru_tests
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
else()
    message(STATUS "doctest not found (set PIRU_FETCH_DOCTEST=ON to fetch); skipping tests.")
endif()
