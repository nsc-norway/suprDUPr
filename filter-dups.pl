#if [ "$1" == "-h" ]
#then
#    echo " -- FASTQ Filter script --"
#    echo ""
#    echo " This script may be used in conjunction with fastdup.read_id"
#    echo " to filter FASTQ files for local duplicates."
#    echo " To output a FASTQ file with duplicates removed to standard"
#    echo " output, use:"
#    echo ""
#    echo " $ ./fastdup.read_id INPUT_FILE | ./filter-dups.sh INPUT_FILE"
#    echo ""
#    echo " (note that the input file has to be specified twice)"
#    
#fi
open(my $fastq_file, $1);
# TODO: this can probably be written better
open(my $stdin, '-');
while(<$stdin>) {
    $_ =~ s/ .*//; # ugh
    my $fqheader = '@' . $_;
    my @lines =  TODO.... 
    while(<$fastq_file>) {
        if TODO field 1 == $fqheader {
            # discard 3 more lines
            break;
        }
        else {
            print;
        }
    }
}

