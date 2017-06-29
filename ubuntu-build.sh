#!/bin/sh
#g++ -fopenmp -std=c++11 -O3 -o duplicate-finder.subrange duplicate-finder.subrange.cpp -lz /usr/lib/x86_64-linux-gnu/libboost_iostreams.so
#g++ -std=c++11 -O3 -o fastdup.read_id fastdup.read_id.cpp -lz -lboost_program_options -lboost_iostreams
g++ -std=c++11 -O3 -o fastdup fastdup.cpp -lz -lboost_program_options -lboost_iostreams
