#!/bin/bash

## Usage: nanoq.sh --min_length <minimum read length> --min_qual <minimum read quality> <input fastq files>

set -ueo

source script_includes.sh

MINLEN=0
MINQUAL=0
FILES=()

while test $# -gt 0
do
	case "$1" in
		# minimum read length (bp)
		--min_length)
			shift
			if [ $# -gt 0 ]; then
				MINLEN=$1
			fi
		;;
		# minimum read quality (bp)
		--min_qual)
			shift
			if [ $# -gt 0 ]; then
				MINQUAL=$1
			fi
		;;
		# input fastq files
		*)
			if [[ "$1" =~ (.fastq|.fq|.fastq.gz|.fq.gz)$ ]]
			then
				FILES+=$1
			else
				echo "Please provide only fastq input files"
				exit 1
			fi
		;;
	esac
	shift
done

echo "Filtering fastq input files"
mkdir -p nanoq
printf "%s\n" ${FILES[@]} \
	| xargs -P $CPUS -n 1 -I{} sh -c '
		shopt -s extglob
		nanoq \
			--min-len $2 \
			--min-qual $3 \
			--output nanoq/$(basename $1).nanoq.fq.gz \
			--input $1
		' -- {} $MINLEN $MINQUAL
echo "Done"
