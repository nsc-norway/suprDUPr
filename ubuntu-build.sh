#!/bin/sh
g++ -fopenmp -std=c++11 -O3 -o duplicate-finder.subrange duplicate-finder.subrange.cpp -lz /usr/lib/x86_64-linux-gnu/libboost_iostreams.so
