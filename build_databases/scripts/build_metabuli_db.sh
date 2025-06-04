#!/bin/bash

HERE="$(dirname $0)"
TAXONOMY="${HERE}/../../taxonomy"
SAMPLESHEET="${HERE}/../inputs/samplesheet.csv"

build_db () {
	mkdir metabuli_tmp

	# Extract list of input FASTA file paths
	cut -d "," -f 3 ${SAMPLESHEET} \
		| tail -n +2 \
		> metabuli_tmp/fastas.txt

	echo "Merging taxonomic accession maps"
	# Merge taxmap files
	cat \
		${TAXONOMY}/nucl_wgs.accession2taxid \
		${TAXONOMY}/nucl_gb.accession2taxid \
	> metabuli_tmp/merged.accession2taxid
	echo "Done"

	mkdir -p results/metabuli_db

	echo "Building metabuli database"
	metabuli build \
		results/metabuli_db \
		metabuli_tmp/fastas.txt \
		metabuli_tmp/merged.accession2taxid \
		--taxonomy-path $TAXONOMY \
		--threads 16
	echo "Done"

	rm -r metabuli_tmp
}

set -ueo
build_db
