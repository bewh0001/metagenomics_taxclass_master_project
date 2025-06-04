#!/bin/bash

set -ueo

HERE="$(dirname $0)"
TAXONOMY="${HERE}/../../taxonomy"
SAMPLESHEET="${HERE}/../inputs/samplesheet.tsv"
DB="${HERE}/../../results/diamond_db/diamond_db.dmnd"

SAMPLES=($(cut -d "," -f 2 $SAMPLESHEET | tail -n +2))

FILES=()
for sample in $SAMPLES
do
	FILES+="$sample"
done

# block-size (b) and index-chunks (c) affect performance and memory use
# peak memory consumption should be roughly 18b/c + 2b (in Gb)

# DIAMOND errors out on some fastq files. Use seqtk seq to convert to fasta
# and pipe to blastx

cd ${HERE}/../../

echo "Aligning reads with DIAMOND blastx"
mkdir -p results/diamond

printf "%s\n" ${FILES[@]} \
	| xargs -n 1 -I{} sh -c '
		seqtk seq -a $1 | diamond blastx \
			--db $2 \
			--out results/diamond/$(basename $1).diamond.tsv \
			--threads 16 \
			--outfmt 102 \
			--long-reads \
			--block-size 6 \
			--index-chunks 1 \
			--query -
		' -- {} $DB
echo "Done"
