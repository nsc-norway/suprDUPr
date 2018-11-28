#!/bin/bash

# Compile the files statically and create a tar ball
make clean
export CFLAGS=-static
make
FILES="filter.sh suprDUPr suprDUPr.read_id filterfq"
tar cfz suprDUPr.tar.gz $FILES

