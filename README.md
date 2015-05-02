SecondCandidate
================================
[![Build Status](https://travis-ci.org/marker68/second-candidate.svg)](https://travis-ci.org/marker68/second-candidate)



## Introduction

This project is the source code of the paper "**Searching for nearest neighbors with a dense space partitioning**".
In this project, we achieved to create a library for nearest neighbors search with the following features:

* An implementation of **Product Quantization for Approximated Nearest Neighbor Search** 
* An implementation of **the inverted multi-index**
* The method that is described in **Searching for nearest neighbors with a dense space partitioning**

## References

## Installation

#### Prerequisites

* An Unix Operating System. Tested on Mac OS X 10.9 and Ubuntu 14.04.
* [CMake](http://www.cmake.org/) 2.8 or newer. For UNIX users, check your CMake version in terminal by `cmake -version`.
* An GNU C++ Compiler or LLVM compiler that support C++ 11 and OpenMP.

#### Build

We use CMake as the build system. On terminal,
```bash
$ git clone git@github.com:your/repo.git pq-learn
$ cd ./pq-learn
$ git submodule update --init
$ cmake -H. -Bbuild && cmake --build build -- -j4
# OR IF YOU WANT TO BUILD TESTS
$ cmake -DBUILD_TEST=ON -H. -Bbuild && cmake --build build -- -j4
```
This script will create binaries in your `bin/` and `lib/` directories. 

## Documentation

We are using Doxygen to generate the documentation. You can find that the style of comments for every methods that are implemented in these source codes are matching with Doxygen style.
To generate the documentation:

```bash
$ cd ./doc
$ doxygen config.doxygen
```

Doxygen will create the documentation files in HTML and LaTeX format. You can also find samples and documentation in the Wiki of this project.

## To developers

* This project was built under Eclipse CDT so we leave the configuration files `.cproject` and `.project` of the project here. If you are using Eclipse, just import the folder that contains the source code.

