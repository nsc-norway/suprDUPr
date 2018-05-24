# suprDUPr

This program was previously known as fastdup, that's why the repo name and the binary are still
called fastdup. We will try to change the name of the binary in the future.

Tool for identifying duplicates reads in close physical proximity on patterned Illumina flow cells.

The tool is designed to analyse FASTQ files with un-trimmed reads, in the format output by BCL
conversion software (e.g. bcl2fastq).

## How to compile the program

make, and a compiler compatible with C++11 are required to build the program. Furthermore, it relies
on the following two Boost libraries, at version 1.48 or newer:

  * program-options
  * iostreams (with gzip support)

### General instructions

You can obtain the source directory by downloading this git repository as a zip file, using the
"Clone or download" button in Github. Alternatively clone it using:

    $ git clone https://github.com/nsc-norway/fastdup

This produces a source directory `fastdup/` which you should then enter. The program can be built using

    $ make

once the dependencies are satisfied.
This will produce binaries `fastdup` and `fastdup.read_id` in the current directory. Run the binaries to get
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


## Programs and pipelines

### Examples

Run with default options:

    $ fastdup data.fastq

Search for duplicates in a region of +/- 5000 pixels in x coordinate and +/- 3000 pixels in y coordinate:

    $ fastdup -x 5000 -y 3000 data.fastq

Process a file which is not sorted:

    $ fastdup --unsorted data.fastq


### Description

  - The main program is called `fastdup`. It examines sequence reads in a
    fastq file (optionally gzip compressed), and computes the fraction
    of reads which are "local" duplicates.
  - A seconary program is `fastdup.read_id`. It outputs part of the FASTQ
    header for all reads identified as duplicates.


### Filtering pipeline

The package also includes a PERL script to filter a file based on the output of
`fastdup.read_id`. For usage information, run:

    ./filter-dups.pl

This command outputs reads which are not "sequencinge duplicates" in data.fastq,
to a file called filtered.fastq:

    $ ./fastdup.read_id data.fastq | ./filter-dups.pl data.fastq > filtered.fastq

The filter-pe script can be used to filter  paired-end reads, but note that only
the first file is considered for the purpose of identifying duplicates. This
script will produce two new fastq files, with the prefix `filtered_`.

    $ ./filter-pe.sh sample_R1.fastq sample_R2.fastq


### Python scripts

The python scripts are included as a kind of documentation. They are only used to simulate how
many duplicates are detected in certain well defined scenarios.


## Bugs / comments

If you find a bug, we would greatly appreciate if you submit an issue on Github.

## Publications

None so far.
