# You should have cmake >= 2.8 installed on your system.
cmake_minimum_required(VERSION 2.8)

# PROJECT'S NAME
project(SecondCandidate CXX)

# THE LOCATION OF OUTPUT BINARIES
set(CMAKE_LIBRARY_DIR ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_DIR})
 
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR})
link_directories(${CMAKE_LIBRARY_DIR})
    
# WE NEED OPENMP
find_package(OpenMP)
find_package(Boost REQUIRED)

# This project uses SIMPLE-CLUSTERS v1.0
# add_subdirectory(lib/simple-cluster)
add_custom_target(
    simplecluster
    COMMAND git submodule update --init &&
    CC=gcc CXX=g++ cmake -H. -Bbuild && 
    cmake --build build && 
    cp lib/libsimplecluster* ${CMAKE_LIBRARY_DIR}/
    WORKING_DIRECTORY ${CMAKE_LIBRARY_DIR}/simple-cluster)
include_directories(${PROJECT_SOURCE_DIR}/lib/simple-cluster/include)

# BUILD VLFEAT
add_custom_target(
    vlfeat
    COMMAND make && cp bin/*/libvl.* ${CMAKE_LIBRARY_DIR}/ && cp bin/*/libvl.* ${CMAKE_BINARY_DIR}/
    WORKING_DIRECTORY ${CMAKE_LIBRARY_DIR}/vlfeat)
include_directories(${CMAKE_LIBRARY_DIR}/vlfeat)

# BUILD OPENBLAS
add_custom_target(
    openblas
    COMMAND make && make PREFIX=${PROJECT_SOURCE_DIR} install
    WORKING_DIRECTORY ${CMAKE_LIBRARY_DIR}/OpenBLAS)
    
# This project uses Google Testing Framework and
# we added GoogleTest as a sub-project. The version
# of this distribution is 1.7.0.
add_subdirectory(lib/googletest)
include_directories(${PROJECT_SOURCE_DIR}/lib/googletest/include)

# For BOOST
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})


# SOURCE FILES
file(GLOB PROJECT_SRCS
    ${PROJECT_SOURCE_DIR}/src/*.cpp
    ${PROJECT_SOURCE_DIR}/lib/simple-cluster/src/*.cpp
)

# COMPILATION FLAGS
if(NOT MSVC)
	if (OPENMP_FOUND)
    	list (APPEND CMAKE_CXX_FLAGS 
    	"-fopenmp-simd -Wopenmp-simd -Wno-write-strings -O3 -std=c++0x ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	else()
    	list(APPEND CMAKE_CXX_FLAGS 
    	"-Wno-write-strings -O3 -std=c++0x ${CMAKE_CXX_FLAGS}")
	endif()
endif()

if(NOT MSVC)
    set(LIBS_FLAGS 
    m pthread simplecluster_static 
    ${CMAKE_LIBRARY_DIR}/libopenblas.a 
    gfortran ${Boost_LIBRARIES}
    vl)
else()
    set(LIBS_FLAGS m simplecluster_static openblas gfortran ${Boost_LIBRARIES} vl)
endif()

# Create a shared library file
add_library(${PROJECT_NAME} SHARED ${PROJECT_SRCS})
if(MSVC)
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
add_library(${PROJECT_NAME}_static STATIC ${PROJECT_SRCS})
if(MSVC)
    set_target_properties(${PROJECT_NAME}_static PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
target_link_libraries(${PROJECT_NAME} ${LIBS_FLAGS})
target_link_libraries(${PROJECT_NAME}_static ${LIBS_FLAGS})
add_dependencies(${PROJECT_NAME} vlfeat openblas)
add_dependencies(${PROJECT_NAME}_static vlfeat openblas)

if(NOT MSVC)
	set(TEST_LIBS_FLAGS 
	${CMAKE_LIBRARY_DIR}/libgtest_main.a 
	m gtest pthread simplecluster_static 
	${CMAKE_LIBRARY_DIR}/libopenblas.a 
	gfortran ${Boost_LIBRARIES}
	vl)
else()
	set(TEST_LIBS_FLAGS gtest simplecluster_static openblas gfortran ${Boost_LIBRARIES} vl)
endif()

## Tests
## set -DBUILD_TEST=ON to build tests
if(${BUILD_TEST})
    include(${PROJECT_SOURCE_DIR}/cmake/build_tests.cmake)
endif(${BUILD_TEST})