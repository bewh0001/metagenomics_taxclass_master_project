#!/bin/bash

set -ueo

HERE="$(dirname $0)"
TAXONOMY="${HERE}/../../taxonomy"
SAMPLESHEET="${HERE}/../inputs/samplesheet.tsv"
DB="${HERE}/../../results/metabuli_db"

SAMPLES=($(cut -d "," -f 2 $SAMPLESHEET | tail -n +2))

FILES=()
for sample in $SAMPLES
do
	FILES+="$sample"
done

cd ${HERE}/../../

echo "Classifying reads with Metabuli"
mkdir -p results/metabuli
printf "%s\n" ${FILES[@]} \
	| xargs -n 1 -I{} sh -c '
		metabuli classify \
			--seq-mode 3 \
			$1 \
			$2 \
			results/metabuli \
			$(basename $1).metabuli \
			--threads 16 \
			--taxonomy-path $3
		' -- {} $DB $TAXONOMY
echo "Done"
