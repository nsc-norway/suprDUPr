#!/usr/bin/env perl

use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
my $fastq_file_path = shift;

if (!$fastq_file_path) {
    print(<<END
    
    == FASTQ Filter script ==

  This script may be used in conjunction with fastdup.read_id
  to filter FASTQ files for local duplicates.
  Quick start: To output a FASTQ file with duplicates removed to
  standard output, use:
 
  \$ ./fastdup.read_id INPUT_FILE | ./filter-dups.pl INPUT_FILE

  (note that the input file has to be specified twice)

   For paired-end reads, see the wrapper script filter-pe.sh.

 Out of any group of duplicate reads, the read occurring at the 
 first position in the file is chosen as the 'original', and is 
 written to the output. The subsequent duplicates are discarded.

END
);
 die('Missing required argument INPUT_FILE.');
}

my $fastq_file = new IO::Uncompress::Gunzip($fastq_file_path,
                                        Transparent => 1,
                                        MultiStream => 1);
open(my $fastdup_input, '-');

my $prev_fqheader = "";
while (<$fastdup_input>) {
    my @parts = split;
    my $fqheader = '@' . @parts[0];
    if ($fqheader eq $prev_fqheader) {
        next; # Skip repeated (duplicate) header lines
    }
    $prev_fqheader = $fqheader;
    my $found = 0;
    while(defined(my $fqdata = <$fastq_file>)) {
        @fqparts = split / /, $fqdata;
        if ($found) {
            if (++$found == 3) {
                last;
            }
        }
        elsif (@fqparts[0] eq $fqheader) {
            $found = 1;
        }
        else {
            print $fqdata;
        }
    }
}

