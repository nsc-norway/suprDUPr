./fastdup.read_id -s 10 -e60 ../Undetermined_S0_L001_R1_001.fastq.gz | sort > fd.txt
./duplicate-finder.subrange ../Undetermined_S0_L001_R1_001.fastq.gz |sort > df.txt
