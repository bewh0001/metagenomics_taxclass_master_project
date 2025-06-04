#!/bin/bash

set -ueo

HERE="$(dirname $0)"
SAMPLESHEET="${HERE}/../inputs/samplesheet.tsv"
DB="${HERE}/../../results/metacache_db"

SAMPLES=($(cut -d "," -f 2 $SAMPLESHEET | tail -n +2))

FILES=()
for sample in $SAMPLES
do
	FILES+="$sample"
done

cd ${HERE}/../../

echo "Querying reads against Metacache database"
mkdir -p results/metacache

queries=""
for file in ${FILES[@]}
do
	queries="${queries}
		-out results/metacache/$(basename $file).metacache.txt
		-abundances results/metacache/$(basename $file).metacache_abundances.txt
		${file}\n"
done

echo -e ${queries} \
	| metacache query \
		${DB}/metacache_db \
		-threads 16 \
		-taxids \
		-separate-cols \
		-lowest species
echo "Done"

