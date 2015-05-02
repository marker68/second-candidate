# SIMPLE CLUSTERS: A simple library for clustering works.
# Copyright (C) 2014 Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# You should have cmake >= 2.8 installed on your system.
# For Ubuntu, just run 'sudo apt-get install cmake'
cmake_minimum_required(VERSION 2.8)

# PROJECT'S NAME
project(simplecluster CXX)

# DETECT THE SYSTEM
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
	message(WARNING "64-bit system is detected")
	set(EX_PLATFORM 64)
	set(EX_PLATFORM_NAME "x64")
elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
	message(WARNING "32-bit system is detected")
	set(EX_PLATFORM 32)
	set(EX_PLATFORM_NAME "x86")
else()
	message(WARNING "An ${CMAKE_SYSTEM_NAME} system is detected")
	set(EX_PLATFORM CMAKE_SYSTEM_NAME)
endif()

# BUILD TYPE
set(CMAKE_BUILD_TYPE Release)
message(WARNING "A ${CMAKE_BUILD_TYPE} build configuration is detected")

# THE LOCATION OF OUTPUT BINARIES
set(CMAKE_LIBRARY_DIR ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_DIR})
 
# INCLUDE DIRECTORIES
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
include_directories("${PROJECT_INCLUDE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}")

# CHECK FOR SUBMODULES AND LIBS
# Check for OpenCV, just for comparing with our methods
set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
#find_package(OpenCV)
find_package(OpenMP)

# SOURCE FILES
set(PROJECT_SRCS 
${PROJECT_SOURCE_DIR}/src/k-means.cpp 
${PROJECT_SOURCE_DIR}/src/utilities.cpp
${PROJECT_SOURCE_DIR}/src/kd-tree.cpp 
)

# COMPLATION FLAGS
if(NOT MSVC)
	if (OPENMP_FOUND)
    	list (APPEND CMAKE_CXX_FLAGS 
    	"-O3 -std=c++0x ${CMAKE_CXX_FLAGS} -g -ftest-coverage -fprofile-arcs ${OpenMP_CXX_FLAGS}")
	else()
    	list(APPEND CMAKE_CXX_FLAGS 
    	"-O3 -std=c++0x ${CMAKE_CXX_FLAGS} -g -ftest-coverage -fprofile-arcs")
	endif()
endif()

# OUTPUT LIBRARIES
add_library(${PROJECT_NAME} SHARED ${PROJECT_SRCS})
if(MSVC)
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
add_library(${PROJECT_NAME}_static STATIC ${PROJECT_SRCS})
if(MSVC)
    set_target_properties(${PROJECT_NAME}_static PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()