#!/bin/bash

# BLAST sample fastq reads against a viral database

HERE="$(dirname $0)"
TAXONOMY="${HERE}/../../taxonomy"
SAMPLESHEET="${HERE}/../inputs/samplesheet.tsv"
DB="nt_viruses"
TAXIDLIST="${HERE}/../inputs/blast_taxids.txt"

BLASTDIR="/sw/data/blast_databases"

cd ${HERE}/../../
mkdir results/blast
cd results/blast

echo "Copying BLAST taxonomy files to working directory. This seems to be necessary for some reason"
cp -n --no-preserve=mode ${BLASTDIR}/taxonomy4blast.sqlite3 .
cp -n --no-preserve=mode ${BLASTDIR}/taxdb.btd .
cp -n --no-preserve=mode ${BLASTDIR}/taxdb.bti .
echo "Done"

mkdir queries

echo "Extracting FASTA queries from input fastq files"
cut -f 2 $SAMPLESHEET \
	| tail -n +2 \
	| xargs -n 1 -P 16 -I{} sh -c '
		seqkit fq2fa -j 1 $1 \
		> queries/$(basename $1).fa
	' -- {}
echo "Done"

# For some reason using more than one thread made BLAST extremely slow,
# so just run queries single-threaded in parallel instead
echo "Blasting reads"
find queries -type f \
	| xargs -P 16 -n 1 -I{} sh -c '
		query=$1
		db=$2
		taxidlist=$3
		out=./$(basename -s ".fa" $query).blast
		blastn -db "$db" -query "$query" \
			-max_target_seqs 10 -word_size 28 \
			-task megablast \
			-perc_identity 90 -evalue 1e-5 \
			-dust yes -qcov_hsp_perc 50 \
			-outfmt "10 qseqid sseqid evalue pident qlen qstart qend sstart send" \
			-out $out -num_threads 1 \
			-taxidlist $taxidlist
	' -- {} $DB $TAXIDLIST
echo "Done"

echo "Cleaning up temporary work files"
rm -f taxonomy4blast.sqlite3 taxdb.btd taxdb.bti
echo "Done"
