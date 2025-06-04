#!/bin/bash

HERE="$(dirname $0)"
TAXONOMY="${HERE}/../../taxonomy"
SAMPLESHEET="${HERE}/../inputs/samplesheet.csv"

cd ${HERE}/../../
mkdir tmp

build_library () {
	echo "Building metacache genomic library from samplesheet $SAMPLESHEET"
	LIBRARY="${PWD}/tmp/library"
	mkdir $LIBRARY

	cut -d "," -f 3 ${SAMPLESHEET} \
		| tail -n +2 \
		| xargs -I{} cp {} $LIBRARY/
	echo "Done"
}

build_db () {
	echo "Linking taxonomy files"
	mkdir tmp/taxonomy
	find $TAXONOMY \
		-maxdepth 1 \
		\( \
			-name "nucl_*.accession2taxid" -o \
			-name "nodes.dmp" -o \
			-name "names.dmp" \
		\) \
		-exec ln -s {} tmp/taxonomy/ \;
	echo "Done"

	echo "Building metacache database"
	mkdir -p results/metacache_db
	cd results/metacache_db

	metacache build \
		metacache_db \
		${LIBRARY} \
		-taxonomy ../../taxonomy/ \
		-taxpostmap ../../taxonomy/*.accession2taxid
	echo "Done"
	cd ..
}

set -ueo
build_library
build_db

echo "Cleaning up temporary files"
rm -r $LIBRARY
echo "Done"
