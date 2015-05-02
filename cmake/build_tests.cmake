add_executable(
    test_quantizer
    ${PROJECT_SOURCE_DIR}/test/test_quantizer.cpp
    ${PROJECT_SOURCE_DIR}/src/quantizer.cpp
    ${PROJECT_SOURCE_DIR}/src/sc_utilities.cpp)
if(MSVC)
    set_target_properties(test_quantizer PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
target_link_libraries(test_quantizer ${TEST_LIBS_FLAGS})
add_dependencies(test_quantizer gtest_main simplecluster_static openblas)

add_executable(
    test_encoder
    ${PROJECT_SOURCE_DIR}/test/test_encoder.cpp
    ${PROJECT_SOURCE_DIR}/src/encoder.cpp
    ${PROJECT_SOURCE_DIR}/src/sc_utilities.cpp)
if(MSVC)
    set_target_properties(test_encoder PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
target_link_libraries(test_encoder ${TEST_LIBS_FLAGS})
add_dependencies(test_encoder gtest_main simplecluster_static openblas)

add_executable(
    test_query
    ${PROJECT_SOURCE_DIR}/test/test_query.cpp
    ${PROJECT_SOURCE_DIR}/src/query.cpp
    ${PROJECT_SOURCE_DIR}/src/sc_utilities.cpp)
if(MSVC)
    set_target_properties(test_query PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
target_link_libraries(test_query ${TEST_LIBS_FLAGS})
add_dependencies(test_query gtest_main simplecluster_static openblas)

add_executable(
    test_multi_query
    ${PROJECT_SOURCE_DIR}/test/test_multi_query.cpp
    ${PROJECT_SOURCE_DIR}/src/query.cpp
    ${PROJECT_SOURCE_DIR}/src/multi_query.cpp
    ${PROJECT_SOURCE_DIR}/src/sc_utilities.cpp)
if(MSVC)
    set_target_properties(test_multi_query PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
target_link_libraries(test_multi_query ${TEST_LIBS_FLAGS})
add_dependencies(test_multi_query gtest_main simplecluster_static openblas)

add_executable(
    test_sc_quantizer
    ${PROJECT_SOURCE_DIR}/test/test_sc_quantizer.cpp
    ${PROJECT_SOURCE_DIR}/src/quantizer.cpp
    ${PROJECT_SOURCE_DIR}/src/sc_quantizer.cpp
    ${PROJECT_SOURCE_DIR}/src/sc_utilities.cpp)
if(MSVC)
    set_target_properties(test_sc_quantizer PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
target_link_libraries(test_sc_quantizer ${TEST_LIBS_FLAGS})
add_dependencies(test_sc_quantizer gtest_main simplecluster_static openblas)


add_executable(
    test_sc_encoder
    ${PROJECT_SOURCE_DIR}/test/test_sc_encoder.cpp
    ${PROJECT_SOURCE_DIR}/src/encoder.cpp
    ${PROJECT_SOURCE_DIR}/src/sc_encoder.cpp
    ${PROJECT_SOURCE_DIR}/src/sc_utilities.cpp)
if(MSVC)
    set_target_properties(test_sc_encoder PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
target_link_libraries(test_sc_encoder ${TEST_LIBS_FLAGS})
add_dependencies(test_sc_encoder gtest_main simplecluster_static openblas)

add_executable(
    test_sc_query
    ${PROJECT_SOURCE_DIR}/test/test_sc_query.cpp
    ${PROJECT_SOURCE_DIR}/src/query.cpp
    ${PROJECT_SOURCE_DIR}/src/sc_query.cpp
    ${PROJECT_SOURCE_DIR}/src/sc_utilities.cpp)
if(MSVC)
    set_target_properties(test_sc_query PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
target_link_libraries(test_sc_query ${TEST_LIBS_FLAGS})
add_dependencies(test_sc_query gtest_main simplecluster_static openblas)

add_executable(
test_evaluation
    ${PROJECT_SOURCE_DIR}/test/test_evaluation.cpp
    ${PROJECT_SOURCE_DIR}/src/evaluation.cpp
    ${PROJECT_SOURCE_DIR}/src/sc_utilities.cpp)
target_link_libraries(test_evaluation ${TEST_LIBS_FLAGS})
if(MSVC)
    set_target_properties(test_evaluation PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
add_dependencies(test_evaluation gtest_main simplecluster_static openblas)
