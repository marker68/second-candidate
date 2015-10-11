SecondCandidate
================================

## Introduction

**[Our slide at ICIP 2015](http://tech-codes.com/wp-content/uploads/2015/09/icip-slide-v2.pdf)**

This project is the source code of the paper "**Searching for nearest neighbors with a dense space partitioning**".
In this project, we try to create a library for nearest neighbors search with the following features:

* An implementation of **Product Quantization for Approximated Nearest Neighbor Search** 
* An implementation of **the inverted multi-index**
* An implementation of the method that is described in **Searching for nearest neighbors with a dense space partitioning**

If you use this source code, please also cite the following reference:

```
@inproceedings{Nguyen15, 
	author={Nguyen, T. and Matsui, Y. and Yamasaki, T. and Aizawa, K.}, 
	booktitle={Proc. ICIP}, 
	title={Searching for nearest neighbors with a dense space partitioning}, 
	year={2015}, 
	pages={1--1}, 
	organization = {IEEE},
}
```

## References
```
[1] H. Jegou, M. Douze and C. Schmid, "Product quantization for nearest neighbor search," IEEE Trans. Pattern Anal. Mach. Intell., vol. 33, no. 1, pp. 117–128, 2011.

[2] A. Babenko and V. Lempitsky, "The inverted multi-index," IEEE Trans. Pattern Anal. Mach. Intell., vol. 37, no. 99, pp. 1247–1260, 2014.

[3] T. A. Nguyen, Y. Matsui, T. Yamasaki and K. Aizawa, "Searching for nearest neighbors with a dense space partitioning," in Proc. ICIP, IEEE, 2015.
```

## Installation

#### Prerequisites

* An Unix Operating System. Tested on Mac OS X 10.9 and Ubuntu 14.04.
* [CMake](http://www.cmake.org/) 2.8 or newer. For UNIX users, check your CMake version in terminal by `cmake -version`.
* An GNU C++ Compiler or LLVM compiler that support C++ 11 and OpenMP.

#### Build

We use CMake as the build system. On terminal,
```bash
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
* This project follows GPLv3 license.

## LICENSE STATEMENT

```
    SecondCandidate - A library for Nearest Neighbor Search
    Copyright (C) 2015  Tuan Anh Nguyen <tuan.nguyenanh@hotmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
```

