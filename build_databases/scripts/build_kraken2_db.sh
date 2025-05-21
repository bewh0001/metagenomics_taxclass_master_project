#!/bin/bash

HERE="$(dirname $0)"
TAXONOMY="${HERE}/../../taxonomy"
SAMPLESHEET="${HERE}/../inputs/samplesheet.csv"

build_db () {
	# NCBI taxonomy files needed to build database
	TAX_FILES=(
		"nucl_gb.accession2taxid"
		"nucl_wgs.accession2taxid"
		"names.dmp"
		"nodes.dmp"
	)

	# Link taxonomy files
	mkdir -p k2_db/taxonomy

	for f in "${TAX_FILES[@]}"; do
		ln -s ${TAXONOMY}/${f} k2_db/taxonomy/${f}
	done

	echo "Adding sequences to Kraken2 library..."
	cut -d "," -f 3 ${SAMPLESHEET} \
		| tail -n +2 \
		| xargs -P 16 -I{} -n 1 \
			k2 add-to-library --db "k2_db" --threads 1 --masker-threads 1 --file {}
	echo "Done"

	echo "Building Kraken2 database"
	kraken2-build \
		--build \
		--db k2_db \
		--threads 16
	echo "Done"

	echo "Cleaning up library build files"
	rm -r k2_db/library
	echo "Done"
}

set -ueo
build_db

