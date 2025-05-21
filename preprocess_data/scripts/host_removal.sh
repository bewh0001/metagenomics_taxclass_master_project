#!/bin/bash

## Usage: host_removal.sh <input fastq files>

HERE=$(dirname $0)
REF_GENOME="${HERE}/../inputs/T2T-CHM13v2.0_reference_genome.fna"
CPUS=16

set -ueo

while test $# -gt 0
do
	case "$1" in
		# input fastq files
		*.fastq|*.fq|*.fastq.gz|*.fq.gz)
			FILES+="$1"
		;;
	esac
	shift
done

# Check if Minimap2 index exists. If not, generate it
INDEX=${REF_GENOME}.mmi

if [ ! -s $INDEX ]
then
	echo "Building reference genome index"
	minimap2 -x map-ont -d $INDEX $REF_GENOME
	echo "Done"
fi

echo "Minimap2: mapping reads to host reference genome"
mkdir -p minimap2
printf "%s\n" ${FILES[@]} \
	| xargs -n 1 -I{} sh -c '
		minimap2 \
			-a \
			-x map-ont \
			-t $1 \
			-o minimap2/$(basename $3).sam \
			$2 $3
		' -- $CPUS $INDEX {}
echo "Done"

echo "Extracting reads"
mkdir -p unmapped
mkdir -p mapped
printf "%s\n" ${FILES[@]} \
	| xargs -n 1 -P $CPUS -I{} bash -c '
		samtools fastq \
			-f 0x4 \
			minimap2/$(basename $1).sam \
			| gzip -c \
			> unmapped/$(basename $1).hostrem.fastq.gz
		samtools fasta \
			-F 0x4 \
			minimap2/$(basename $1).sam \
			| grep "^>" \
			| cut -c 2- \
			> mapped/$(basename $1).host_mapped_reads.txt
		' -- {}
echo "Done"
