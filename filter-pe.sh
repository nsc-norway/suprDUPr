#!/bin/env bash

if [[ "$2" == "" ]]
then
	echo "usage: filter-pe.sh R1_FILE R2_FILE"
	exit 1
fi

fifo=`mktemp`
mkfifo $fifo
./filter-dups.pl $1 > filtered_$1 < $fifo &
./filter-dups.pl $2 > filtered_$2 < $fifo &
./fastdup.read_id $1 > $fifo &
wait

