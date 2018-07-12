#!/usr/bin/env bash

single_read=false
supr_args=()
while [[ ! -z "$1" ]];
do
	if [[ $1 == "-1" ]]
	then
		single_read=true
	elif [[ $1 == -* ]]
	then
		supr_args+=($1)
		shift
	else
		break
	fi
done

if $single_read
then
	r1file="$1"
	r1output="$2"
else
	r1file="$1"
	r1output="$3"
	r2file="$2"
	r2output="$4"
fi

if [[ ! -f "$r1file" ]]
then
	echo "Error: '$r1file' is not a file."
	exit 1
fi

# Re-compress if outputs are compress, pipe through gzip
if [[ "$r1output" == *.gz ]]
then
	outputfifo1=`mktemp -u`_outputfifo1
	mkfifo $outputfifo1
	gzip -c > "$r1output" < $outputfifo1 &
	r1output=$outputfifo1
fi
if [[ "$r2output" == *.gz ]]
then
	outputfifo2=`mktemp -u`_outputfifo2
	mkfifo $outputfifo2
	gzip -c > "$r2output" < $outputfifo2 &
	r2output=$outputfifo2
fi

if [[ -z "$1" ]]
then
	echo ""
	echo "usage: filter-pe.sh [-1] [SUPR_OPTIONS] R1_FILE [R2_FILE] R1_OUTPUT [R2_OUTPUT]"
	echo ""
	echo "       Use option -1 to process single-read data."
	echo ""
	exit 1
fi


if $single_read;
then
	./suprDUPr.read_id ${supr_args[*]} "$r1file" | ./filter-dups.pl "$r1file" > "$r1output"
else
	read_id_fifo=`mktemp -u`_readidfifo
	mkfifo $read_id_fifo
	./suprDUPr.read_id ${supr_args[*]} "$r1file" | tee $read_id_fifo | ./filter-dups.pl "$r1file" > "$r1output" &
	./filter-dups.pl "$r2file" > "$r2output" < $read_id_fifo &
	rm $read_id_fifo $outputfifo1 $outputfifo2
	wait
fi

