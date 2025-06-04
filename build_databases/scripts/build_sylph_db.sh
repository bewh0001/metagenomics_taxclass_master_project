#!/bin/bash

HERE="$(dirname $0)"
TAXONOMY="${HERE}/../../taxonomy"
SAMPLESHEET="${HERE}/../inputs/samplesheet.csv"
VIRAL_GENOMES="${HERE}/../inputs/sylph_viral_genomes.txt"

build_db () {
	cd ${HERE}/../../

	mkdir results/sylph_db
	cd results/sylph_db

	echo "Building Sylph database"
	cut -d "," -f 3 ${SAMPLESHEET} \
		| tail -n +2 \
		> fastas.txt

	echo "Building separate viral database with c=1 subsampling rate"
	(fgrep -f $VIRAL_GENOMES fastas.txt || true;) \
		> viral_fastas.txt
	(fgrep -v -f viral_fastas.txt fastas.txt || true;) \
		> other_fastas.txt

	sylph sketch \
		-c 1 \
		-l viral_fastas.txt \
		-t 16 \
		-o viral_database
	echo "Done"

	echo "Building non-viral database with c=100"
	sylph sketch \
		-c 100 \
		-l other_fastas.txt \
		-t 16 \
		-o non_viral_database
	echo "Done"

	echo "Generating Sylph taxonomy"
	python3 ${HERE}/../../taxtools/build_sylph_taxonomy.py \
		--samplesheet $SAMPLESHEET \
		--names_dmp ${TAXONOMY}/names.dmp \
		--nodes_dmp ${TAXONOMY}/nodes.dmp \
		--use-taxids \
		--output ${PWD}/sylph_tax_ids.tsv 
	echo "Done"
}

set -ueo
build_db
