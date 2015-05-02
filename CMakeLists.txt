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

# This project uses SIMPLE-CLUSTERS v1.0
add_subdirectory(lib/simple-cluster)
include_directories(${PROJECT_SOURCE_DIR}/lib/simple-cluster/include)
    
# WE NEED OPENMP
find_package(OpenMP)
#find_package(Boost REQUIRED)

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
    	"-Wno-write-strings -O3 -std=c++0x ${CMAKE_CXX_FLAGS} -g -ftest-coverage -fprofile-arcs ${OpenMP_CXX_FLAGS}")
	else()
    	list(APPEND CMAKE_CXX_FLAGS 
    	"-Wno-write-strings -O3 -std=c++0x ${CMAKE_CXX_FLAGS} -g -ftest-coverage -fprofile-arcs")
	endif()
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

# BUILD OPENBLAS
add_custom_target(
    openblas
    COMMAND make && make PREFIX=${PROJECT_SOURCE_DIR} install
    WORKING_DIRECTORY ${CMAKE_LIBRARY_DIR}/OpenBLAS)
add_dependencies(${PROJECT_NAME} openblas)
add_dependencies(${PROJECT_NAME}_static openblas)

# OUTPUT EXECUTABLES
if(NOT MSVC)
	set(TEST_LIBS_FLAGS ${CMAKE_LIBRARY_DIR}/libgtest_main.a m gtest pthread simplecluster_static ${CMAKE_LIBRARY_DIR}/libopenblas.a gfortran ${Boost_LIBRARIES})
else()
	set(TEST_LIBS_FLAGS gtest simplecluster_static openblas gfortran ${Boost_LIBRARIES})
endif()

## Tests
## set -DBUILD_TEST=ON to build tests
if(${BUILD_TEST})
    include(${PROJECT_SOURCE_DIR}/cmake/build_tests.cmake)
endif(${BUILD_TEST})