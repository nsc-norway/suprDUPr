# suprDUPr

Tool for identifying duplicates reads in close physical proximity on patterned Illumina flow cells.

The tool is designed to analyse FASTQ files with un-trimmed reads, in the format output by BCL
conversion software (e.g. bcl2fastq).

## Programs and pipelines

### Examples

Run with default options:

    $ suprDUPr data.fastq

Search for duplicates in a region of +/- 5000 pixels in x coordinate and +/- 3000 pixels in y coordinate:

    $ suprDUPr -x 5000 -y 3000 data.fastq

Process a file which is not sorted:

    $ suprDUPr --unsorted data.fastq


### Description

  - The main program is called `suprDUPr`. It examines sequence reads in a
    fastq file (optionally gzip compressed), and computes the fraction
    of reads which are "local" duplicates.
  - A seconary program is `suprDUPr.read_id`. It outputs part of the FASTQ
    header for all reads identified as duplicates.


### Filtering pipeline

The package also includes a PERL script to filter a file based on the output of
`suprDUPr.read_id`. For usage information, run:

    ./filter-dups.pl

This command outputs reads which are not "sequencinge duplicates" in data.fastq,
to a file called filtered.fastq:

    $ ./suprDUPr.read_id data.fastq | ./filter-dups.pl data.fastq > filtered.fastq

The filter-pe script can be used to filter  paired-end reads, but note that only
the first file is considered for the purpose of identifying duplicates. This
script will produce two new fastq files, with the prefix `filtered_`.

    $ ./filter-pe.sh sample_R1.fastq sample_R2.fastq


### Paired-end analysis

Usually it is sufficient to analyse just one of the reads in paired-end data.
Analysis of both reads may be useful if the library has been sequenced very
deeply. By analysing both ends, one can avoid counting fragments which have the
same sequence in read one, but a different insert size. These fragments are
not caused by "sequencing-related" effects, such as ExAmp -- instead they are
genuinely different reads.

It is possible to use common Unix commands, and the suprDUPr.read_id program,
to count only reads which are duplicates in both read 1 and read 2.
Assuming we have two files, R1.fastq and R2.fastq, the following command
counts paired-end duplicates:

    comm -12 <( ./suprDUPr.read_id R1.fastq | cut -f1 | uniq | sort ) <( ./suprDUPr.read_id R2.fastq | cut -f1 | uniq | sort ) | wc -l

It will produce the normal suprDUPr output for both files, and then the last
line of output contains the number of paired-end duplicates. If working with
unsorted fastq files, you may have to swap the `sort` and `uniq` commands
around.


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



### Model

Some statistical models were developed to estimate the effect of global duplicates
on the estimate for local duplicates. For example, how many randomly positioned 
duplicates are necessary in order to get an appreciable result from suprDUPr?

These scripts are located in the model/ subdirectory. They are not useful for 
analysing data.


## Bugs / comments

If you find a bug, we would greatly appreciate if you submit an issue on Github.

## Publications

None so far.
