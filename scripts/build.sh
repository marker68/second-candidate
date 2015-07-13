#!/bin/bash
CMAKE=$2
CC=$3
CXX=$4

if [ "$1" == "all" ]; then
	CC=${CC} CXX=${CXX} ${CMAKE} -DBUILD_TEST=ON -Bbuild -H.
	${CMAKE} --build build -- -j4
elif [ "$1" == "clean" ]; then
	rm -rf bin
	rm -rf build
	rm -f lib/lib*
	cd lib/OpenBLAS
	make clean
	cd ../vlfeat
	make clean
	cd ../googletest
	make clean
	cd ../simple-cluster
	./clean.sh `pwd`
fi
	
