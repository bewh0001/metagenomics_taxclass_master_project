#!/bin/bash

## Usage: subsample_reads.sh <input fastq files>

set -ueo

BASES="5000MB"
CPUS=16

while test $# -gt 0
do
	case "$1" in
		# input fastq files
		*)
			if [[ "$1" =~ (.fastq|.fq|.fastq.gz|.fq.gz)$ ]]
			then
				FILES+=($(abs_path $1))
			else
				echo "Please provide only fastq input files"
				exit 1
			fi
		;;
	esac
	shift
done

echo "Subsampling input reads to $BASES bases"
mkdir -p rasusa
printf "%s\n" ${FILES[@]} \
	| xargs -P $CPUS -n 1 -I{} sh -c '
		rasusa reads \
			-o rasusa/$(basename $1).subsampled.fq.gz \
			-b $2 \
			$1
		' -- {} $BASES
echo "Done"
