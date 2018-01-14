#!/bin/bash

# Standard benchmark to make sure we don't make it much slower
# or introduce bugs.

time ./fastdup ../../nsc/data/Undetermined_S0_L001_R1_001.fastq.gz 
### NUM_READS   READS_WITH_DUP  DUP_RATIO
### 12192602    356128  0.0292085
### 
### real    0m50,594s
### user    0m49,819s
### sys 0m0,576s


