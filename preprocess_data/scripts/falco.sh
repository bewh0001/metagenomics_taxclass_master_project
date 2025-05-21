#!/bin/bash

## Usage: falco.sh <input fastq files>

set -ueo

HERE="$(dirname $0)"
CPUS=16
FILES=()

while test $# -gt 0
do
	case "$1" in
		*)
			if [[ "$1" =~ (.fastq|.fq|.fastq.gz|.fq.gz)$ ]]
			then
				FILES+="$1"
			fi
		;;
	esac
	shift
done

echo "Compiling FastQC reports from input files"
mkdir -p fastqc
cd fastqc
printf "%s\n" ${FILES[@]} \
	| xargs -P $CPUS -n 1 -I{} sh -c '
		falco -o $(basename $1) $1
		' -- {}
echo "Done"
