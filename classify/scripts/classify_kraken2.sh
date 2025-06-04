#!/bin/bash

set -ueo

HERE="$(dirname $0)"
SAMPLESHEET="${HERE}/../inputs/samplesheet.tsv"
DB="${HERE}/../../results/kraken2_db"

SAMPLES=($(cut -d "," -f 2 $SAMPLESHEET | tail -n +2))

FILES=()
for sample in $SAMPLES
do
	FILES+="$sample"
done

# Copy database to shared memory for memory mapping"
SHM_DB=/dev/shm/${SLURM_JOB_ID}/k2-db
echo "Mapping database to shared memory (${SHM_DB})"
cp -r $DB $SHM_DB
echo "Done"

cd ${HERE}/../../

echo "Classifying reads with Kraken2"
mkdir -p results/kraken2
printf "%s\n" ${FILES[@]} \
	| xargs -n 1 -I{} sh -c '
		k2 classify \
			--db $2 \
			--threads 16 \
			--memory-mapping \
			--output results/kraken2/$(basename $1).k2_output.tsv \
			--report results/kraken2/$(basename $1).k2_report.tsv \
			$1
		' -- {} $SHM_DB $CPUS $CONFIDENCE
echo "Done"

echo "Unmapping database from shared memory"
rm -r $SHM_DB
echo "Done"
