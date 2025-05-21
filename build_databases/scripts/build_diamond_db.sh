#!/bin/bash

HERE="$(dirname $0)"
TAXONOMY="${HERE}/../../taxonomy"
SAMPLESHEET="${HERE}/../inputs/samplesheet.csv"

build_db () {
echo "Building DIAMOND protein library from samplesheet $SAMPLESHEET"
	LIBRARY="$(pwd)/library.faa"

	# Merge input FASTAs
	cut -d "," -f 4 $SAMPLESHEET \
		| tail -n +2 \
		| xargs cat \
			> ${LIBRARY}.duplicated.faa
	echo "Done"

	# Remove duplicate sequences with SeqKit
	echo "Deduplicating library"
	seqkit rmdup ${LIBRARY}.duplicated.faa \
		> "$LIBRARY"
	rm ${LIBRARY}.duplicated.faa
	echo "Done"

	echo "Building DIAMOND database"
	mkdir -p diamond_db
	cd diamond_db

	diamond makedb \
		--in ${LIBRARY} \
		--db "diamond-db" \
		--taxonmap ${TAXONOMY}/prot.accession2taxid.FULL \
		--taxonnodes ${TAXONOMY}/nodes.dmp \
		--taxonnames ${TAXONOMY}/names.dmp \
		--threads 16
	cd ..

	echo "Cleaning up temporary files"
		rm ${LIBRARY}
	echo "Done"
}

set -ueo

build_db
