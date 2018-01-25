# fastdup

Tool for identifying duplicates reads in close physical proximity on patterned Illumina flow cells.

The tool is designed to analyse FASTQ files with un-trimmed reads, in the format output by BCL
conversion software (e.g. bcl2fastq). 

## Installation

A compiler compatible with C++11 is required to build the program. Furthermore, it relies on the
following two Boost libraries:

  * program-options
  * iostreams (with gzip support)

The program can be built by entering the source directory and issuing the command:

    make

This will produce binaries `fastdup` and `fastdup.read_id` in the current directory. Run the binaries to get
a list of options.

