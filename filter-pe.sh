#!/usr/bin/env bash

if [[ "$2" == "" ]]
then
	echo "usage: filter-pe.sh R1_FILE R2_FILE"
	exit 1
fi

fifo=`mktemp -u`
mkfifo $fifo
./suprDUPr.read_id $1 > $fifo &
./filter-dups.pl $1 > filtered_`basename ${1%.gz}` < $fifo &
./filter-dups.pl $2 > filtered_`basename ${2%.gz}` < $fifo &
wait

