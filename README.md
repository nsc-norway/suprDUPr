## fastdup

Tool for identifying duplicates reads in close physical proximity on patterned Illumina flow cells.

The tool is designed to analyse FASTQ files with un-trimmed reads, in the format output by BCL conversion software (e.g. bcl2fastq). 

# Installation

The tool requires a compiler compatible C++11 


# Alternative programs

The programs ''duplicate-finder.subrange.cpp'' and ''fastdup.read_id.cpp'' are used to investigate the subtleties of
identified duplicates in fastq files. They are not intended for production use. Both programs print the "read-id" 
of each duplicate, which is is the first part of the fastq header. For any group of equal sequences found, all of 
the sequences are printed (this is different from the MarkDuplicates tool in Picard Tools, which keeps one "primary"
read as non-duplicate, and marks others as duplicate).

  * fastdup.read_id: Works like the standard fastdup tools, but outputs the read-ID.
  * duplicate-finder.subrange: Works like fastdup.read_id, but uses the Levenshtein distance to compare sequences,
			instead of requiring an exact match. Takes most options at compile time (again, it's not
			production ready). At least an order of magnitude slower, but SMP-enabled.

