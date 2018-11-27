#!/bin/bash

# Compile the files statically and create a tar ball
make clean
export CFLAGS=-static
make
tar cf suprDUPr.tar suprDUPr suprDUPr.read_id filterfq

