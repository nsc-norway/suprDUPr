# suprDUPr

Tool for identifying duplicates reads in close physical proximity on patterned Illumina flow cells.

The tool is designed to analyse FASTQ files with un-trimmed reads, in the format output by BCL
conversion software (e.g. bcl2fastq).

## News

  - Version 1.3: Fixed some compilation errors (!) and provided downloadable
                 binaries.
  - Version 1.2: Added multithreaded input. See Multithreading section
  - Version 1.1: Introduces native support for analysing both reads of paired-
                 end data. Because two command-line arguments are needed to 
                 specify the input files, it is no longer possible to specify
                 an output file on the command line. Instead it is recommended
                 to redirect STDOUT to a file using other methods.


## Programs and pipelines

### Description

  - The main program is called `suprDUPr`. It examines sequence reads in a
    fastq file (optionally gzip compressed), and computes the fraction
    of reads which are "local" duplicates.
  - A seconary program is `suprDUPr.read_id`. It outputs part of the FASTQ
    headers for pairs of reads identified as duplicates.
  - See the scripts below for more functionality.


### Main program usage

The options for `suprDUPr` and `suprDUPr.read_id` are the same.

    usage: ./suprDUPr [options] input_file_r1 [input_file_r2]

    Allowed options:
      -x [ --winx ] arg (=2500)  x coordinate window, +/- pixels
      -y [ --winy ] arg (=2500)  y coordinate window, +/- pixels
      -s [ --start ] arg (=10)   First position in reads to consider
      -e [ --end ] arg (=60)     Last position in reads to consider
      -r [ --region-sorted ]     Assume the input file is sorted by region (tile),
                                 but not by (y, x) coordinate within the region.
      -u [ --unsorted ]          Process unsorted file (a large hash-size is
                                 recommended, see --hash-size). This mode requires
                                 all data to be stored in memory, and it is not
                                 well optimised.
      -1 [ --single ]            Disable multithreading
      --hash-size arg (=4194304) Hash table size (bytes), must be a power of 2.
                                 (increase if winy>2500).
      -h [ --help ]              Show this help message
    
    Specify - for input_file_r1 to read from stdin.

If two files are provided on the command line, they are assumed to be from paired-end
sequencing. Then the same substring is always used in both reads, and `-s` and `-e`
specify the positions to use in both reads.

There was previously an option for output filename instead of `input_file_r2`.
**If you are going to analyse PE data, make sure that the suprDUPr version is
greater or equal to 1.1, or it may instantly overwrite your second FASTQ file,
as an output file.**


### Basic examples

Run with default options:

    $ suprDUPr data.fastq

Search for duplicates in a region of +/- 5000 pixels in x coordinate and +/- 3000 pixels in y coordinate:

    $ suprDUPr -x 5000 -y 3000 data.fastq

Use the first 60 bases in each read for comparison. Note that no mismatches are allowed
in this region:

    $ suprDUPr -s 0 -e 59 data.fastq

Process a file which is not sorted:

    $ suprDUPr --unsorted data.fastq

### Example output

    $ ./suprDUPr TEST_R1.fq.gz
    -- suprDUPr v1.1 --
    Using positions from 10 to 60
    Started reading FASTQ file...
    Completed. Analysed 10000 records.
    NUM_READS	READS_WITH_DUP	DUP_RATIO
    10000	233	0.0233

The two bottom lines are written to STDOUT, the rest is status information written
to STDERR. The table contains three columns:

 - `NUM_READS`: Number of reads analysed. For PE, it's the number of read pairs.
 - `READS_WITH_DUP`: Number of reads identified as duplicate. See below.
 - `DUP_RATIO`: Fraction of the reads which are duplicate. This is usually the most
   interesting. This is `READS_WITH_DUP` divided by `NUM_READS`.

#### Definition of duplicate reads

The simplest is if duplicates only occur as pairs of identical sequences. Then one
read out of each pair is counted as duplicate. When there are more than two, all
reads except one are duplicates. The column is named `READS_WITH_DUP` because of an
equivalent method of counting: the program counts all reads for which we already know
about a read with identical sequence (within the search area), a "duplicate". Each
examined pair is required to be within the distance threshold.

#### Read-identifier output

The alternative program `suprDUPr.read_id` can be used for further analysis of
duplicates. It outputs a line for each pair of sequences which are identical in the 
defined substring match, if they are within the distance threshold. Each line contains
two tab-separated values, with a substring of the FASTQ headers of the identical
sequences.


    ST-E00159:10:H05PHALXX:1:1101:27126:1326	ST-E00159:10:H05PHALXX:1:1101:27105:1326
    ST-E00159:10:H05PHALXX:1:1101:29227:1344	ST-E00159:10:H05PHALXX:1:1101:29196:1326
    ST-E00159:10:H05PHALXX:1:1101:30303:1449	ST-E00159:10:H05PHALXX:1:1101:30282:1449
    ST-E00159:10:H05PHALXX:1:1101:31683:1555	ST-E00159:10:H05PHALXX:1:1101:31663:1555
    ST-E00159:10:H05PHALXX:1:1101:27288:1748	ST-E00159:10:H05PHALXX:1:1101:27298:1731
    ST-E00159:10:H05PHALXX:1:1101:28912:1748	ST-E00159:10:H05PHALXX:1:1101:27796:1432
    ST-E00159:10:H05PHALXX:1:1101:18833:1766	ST-E00159:10:H05PHALXX:1:1101:19107:1748
    ST-E00159:10:H05PHALXX:1:1101:22203:1836	ST-E00159:10:H05PHALXX:1:1101:22142:1801

This output will usually contain more lines than `READS_WITH_DUP`, if there are any
n-plicates in the file for n>2. The output is only for pairwise comparisons, but a
single read may produce multiple lines if it is equal to several other reads.  The
number of unique items in either of the columns is equal to `READS_WITH_DUP`. The
first column always contains the read which occurs later in the file. `suprDUPr.read_id`
will output the normal statistics table to STDERR, as STDOUT is reserved for the read
identifiers.


### Multithreading

For gzip'd input files, the decompression runs in separate threads by default.  If
the input is not compressed, only a single thread is used. gzip decompression requires
more CPU than the analysis, so the core for analysis is never fully utilised. For
single-read data, the CPU usage is around 140 % of a core, for paired end it is
240 %. If only one core is available, it may be slightly more efficient to constrain
it to single-threaded operation using the -1 option.


### Filtering pipeline (duplicate removal)

The `filter.sh` script can be used to remove duplicates in single-read or paired-end data.
Note that only the first file (Read 1) is considered for the purpose of identifying 
duplicates. It uses the `suprDUPr.read_id` program and a C++ program `filterfq`
to perform the duplicate removal. `filter.sh` itself is mainly a wrapper program,
which uses a named pipe to apply the same "read-ID" filtering to both data files
in case of paired-end data.
Example:

    $ ./filter.sh sample_R1.fastq sample_R2.fastq filtered_R1.fastq filtered_R2.fastq

Produces files:

    filtered_R1.fastq    filtered_R2.fastq

The input of `filter.sh` must be real files, not a pipe. The read-1 file is read twice in
parallel, by `suprDUPr.read_id` and by `filterfq`.

If the inputs have extension ".gz", the input will be treated as compressed, and the
output will also be compressed, regardless of its extension.


### Advanced filtering script

The package includes a program to remove duplicates in a file based on the output of
`suprDUPr.read_id`, called `filterfq`.

This program is used by the above wrapper script `filter.sh`, but it can also be used directly.

The following command outputs reads which are not "sequencinge duplicates" in data.fastq,
to a file called filtered.fastq:

    $ ./suprDUPr.read_id data.fastq | ./filterfq data.fastq > filtered.fastq

Note that the file must be specified as an input to both `suprDUPr.read_id` and
`filterfq`.


### Paired-end analysis

Usually it is sufficient to analyse just one of the reads in paired-end data.
Analysis of both reads may be useful if the library has been sequenced very
deeply. By analysing both ends, one can avoid counting fragments which have the
same sequence in read one, but a different insert size. These fragments are
not caused by "sequencing-related" effects, such as ExAmp -- instead they are
genuinely different reads.

Previously, PE analysis was implemented using `suprDUPr.read_id` and standard
UNIX shell commands. As suprDUPr now has native PE support, this approach is no
longer needed. Keep in mind that a shorter matching string length per read
(start/end position) may be in required, as the effective matching string is
twice as long when running on PE data.


## How to compile the program

make, and a compiler compatible with C++11 are required to build the program. Furthermore, it relies
on the following two Boost libraries, at version 1.48 or newer:

  * program-options
  * iostreams (with gzip support)

### General instructions

You can obtain the source directory by downloading this git repository as a zip file, using the
"Clone or download" button in Github. Alternatively clone it using:

    $ git clone https://github.com/nsc-norway/suprDUPr

This produces a source directory `suprDUPr/` which you should then enter. The program can be built using

    $ make

once the dependencies are satisfied.
This will produce binaries `suprDUPr` and `suprDUPr.read_id` in the current directory. Run the binaries to get
a list of options, and see the "Programs and Pipelines" section below for descriptions.
There is no installation script -- you can copy the executable to `/usr/bin` or some other directory on
the system PATH.


### CentOS 7 / RHEL 7 / Fedora

Run yum to install required packages (it will skip them if they are already installed):

    # yum install make gcc-c++ zlib-devel boost-{devel,iostreams,program-options}

The boost-devel package is just needed for build, not to run the program. If you already have a C++ compiler
you like, you can omit gcc-c++. Then run `make` in the source directory.

Tip: if in a mixed environment, the software won't run on different versions, due to boost version conflicts.
The boost148 packages can be used, as they are available on both versions of RHEL/CentOS. For example, build
it on Centos 6 as described below.


### Ubuntu / Debian

This works at least as far back as Ubuntu 16.04, Debian stretch (9), maybe also earlier versions.

    $ sudo apt install make g++ libboost-{iostreams,program-options}-dev libz-dev
    $ make

This pulls in the boost libraries in addition to the development headers, as dependencies. You can omit the -dev
suffix in the command above to just install libraries, to run but not build suprDUPr.


### macOS

Install boost using one of the package systems for Mac, or build directly from source. See this question on
Stackoverflow for how to install boost:

https://stackoverflow.com/questions/104322/how-do-you-install-boost-on-macos

One of these commands might do it, if you have a package system such as MacPorts/Homebrew installed:

    $ sudo ports install boost
    $ brew install boost
    $ fink install boost # ?

You may also need to install the package libz. Specify the package system's include and library directory, or if
you installed boost directly, you may need to specify its directories. For MacPorts it's:

    CFLAGS="-I/opt/local/include -L/opt/local/lib" make

I believe this may work for Homebrew:

    CFLAGS="-I/usr/local/include -L/usr/local/lib" make

In case you get the following error, the boost library files may be installed with a suffix:

    ld: library not found for -lboost_program_options

Try this (use a CFLAGS that worked before):

    CFLAGS="-I/opt/local/include -L/opt/local/lib" BOOST_LIB_SUFF="-mt" make


### CentOS 6 / RHEL 6

It is possible to use scl (Software Collections) to get access to a recent enough compiler. A recent enough
boost can be found in EPEL, so EPEL is also required. First, install SCL and EPEL:

    # yum install epel-release centos-release-scl		# CentOS
    # yum install epel-release scl-utils			# RHEL (may require others)

Then install the devtoolset scl and other dependencies:

    # yum install devtoolset-4 boost148-{devel,iostreams,program-options} zlib-devel

After installing dependencies, enter a scl subshell and compile it:

    $ scl enable devtoolset-4 bash
    $ CFLAGS="-I/usr/include/boost148 -L/usr/lib64/boost148" make
    $ exit



## Model

Some statistical models were developed to estimate the effect of global duplicates
on the estimate for local duplicates. For example, how many randomly positioned 
duplicates are necessary in order to get an appreciable result from suprDUPr?

These scripts are located in the model/ subdirectory. They are not useful for 
analysing data.


## Bugs / comments

If you find a bug, we would greatly appreciate if you submit an issue on Github.

## Publications

None so far.
