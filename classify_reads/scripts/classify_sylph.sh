#!/bin/bash

set -ueo

HERE="$(dirname $0)"
TAXONOMY="${HERE}/../../taxonomy"
SYLPH_TAX="${HERE}/../../results/sylph_db/sylph_tax_ids.tsv"
SAMPLESHEET="${HERE}/../inputs/samplesheet.tsv"
DB="${HERE}/../../results/sylph_db"

ASM_SUMMARY=${TAXONOMY}/assembly_summary_refseq.txt
MIN_COVERAGE_DEPTH=1

SAMPLES=($(cut -d "," -f 2 $SAMPLESHEET | tail -n +2))

FILES=()
for sample in $SAMPLES
do
	FILES+="$sample"
done

cd ${HERE}/../../

echo "Sketching sample reads"
mkdir -p results/sylph/reads
cd sylph/reads

sylph sketch \
	-t 16 \
	-c 100 \
	-r ${FILES[@]} 
cd ..
echo "Done"

echo "Profiling reads with Sylph"
sylph profile \
	$DB/*.syldb \
	reads/*.sylsp \
	-u \
	--min-number-kmers 3 \
	--minimum-ani 95 \
	-t 16 \
	-o sylph.results.tsv
echo "Done"

echo "Generating taxonomic Sylph profile"
sylph-tax taxprof \
	sylph.results.tsv \
	-t $SYLPH_TAX
echo "Done"

echo "Estimating read counts"
for SAMPLE in "${FILES[@]}"
do
	SAMPLE_BASE=$(basename $SAMPLE)
	if [ ! -f ${SAMPLE_BASE}.sylphmpa ]; then
		continue
	fi

	# get total number of sample reads
	N_READS=$(seqtk seq -A $SAMPLE | grep "^>" -c)

	# get estimated read counts from sequence abundances
	echo "read_count" > ${SAMPLE_BASE}.read_counts.tmp
	tail -n +3 ${SAMPLE_BASE}.sylphmpa \
		| awk -v n_reads="$N_READS" \
			'{printf "%1.0f\n", $3 * n_reads / 100}' \
		>> ${SAMPLE_BASE}.read_counts.tmp

	# add read counts to taxonomic output
	head -n 1 ${SAMPLE_BASE}.sylphmpa \
		> ${SAMPLE_BASE}.with_read_counts.sylphmpa
	paste \
		<(tail -n +2 ${SAMPLE_BASE}.sylphmpa | cut -f 1-3) \
		${SAMPLE_BASE}.read_counts.tmp \
		<(tail -n +2 ${SAMPLE_BASE}.sylphmpa | cut -f 4-) \
		>> ${SAMPLE_BASE}.with_read_counts.sylphmpa
	rm ${SAMPLE_BASE}.read_counts.tmp
done
echo "Done"

mkdir minimap2

echo "Mapping sample reads to positive genomes to extract positive reads"
mkdir -p minimap2/idx
mkdir -p queries

while read -r fastq genome
do
	sample=$(basename $fastq | cut -d "." -f 1)

	# get species taxid of genome assembly
	taxid=$({ grep -w $(basename -s "_genomic.fna" $genome) \
		$ASM_SUMMARY || true; } | cut -f 7)

	if [ -z "$taxid" ]; then
		echo "Genome $genome not found in assembly summary file. Skipping"
		continue
	fi

	# index genome once and reuse
	genome_idx="minimap2/idx/$(basename -s .fna $genome).mmi"
	if [ ! -f "$genome_idx" ]; then
		echo "Indexing genome"
		minimap2 -x map-ont -d $genome_idx $genome
		echo "Done"
	fi

	bam="minimap2/bam/${sample}.${taxid}/${sample}.${taxid}.$(basename -s "_genomic.fna" $genome).bam"
	mkdir -p $(dirname $bam)

	echo "Mapping sample reads to genome with minimap2"
	minimap2 -a -x map-ont -t 16 \
		$genome_idx $fastq \
		| samtools view -b - \
		| samtools sort -m 7G -o $bam -
	echo "Done"

	echo "Building samtools index"
	samtools index $bam
	echo "Done"

	echo "Getting genome length"
	genome_len=$(samtools view -H $bam \
		| grep '^@SQ' \
		| cut -d ':' -f 3 \
		| awk '{sum+=$1} END {print sum}')
	echo "Done"

	echo "Getting taxon name and rank"
	read -r taxon_name taxon_rank \
		<<< "$(taxonkit lineage -L -n -r <<< $taxid \
			| cut -f 2-2 | tr -s ' ' '_')"
	echo "Done"

	echo "Extracting reads mapped to genome"
	FASTA="${PWD}/queries/${sample}.taxid_${taxid}.fa"
	samtools fasta \
		-F 0x4 \
		$bam \
		>> $FASTA
	echo "Done"

	echo "Extracting read mapping scores"
	samtools view -F 0x4 $bam \
		| cut -f 1,5 \
		> ${bam%.bam}.mapq.tsv
	echo "Done"

	echo "$sample" "$taxid" "$FASTA" \
		| column -t -o $'\t' \
		>> query_list.txt.tmp

	mv ${bam%.bam}.coverage.tsv minimap2/
	mv ${bam%.bam}.mapq.tsv minimap2/
done < <(tail -n +2 sylph.results.tsv | cut -f 1-2)
echo "Done"

sort -u query_list.txt.tmp > query_list.txt
rm query_list.txt.tmp

echo "Extracting mapped reads per sample and taxid"
mkdir mapped_reads

for sample in $(cut -f 1 query_list.txt | sort -u)
do
	echo "read_id" "taxid" \
		| column -t -o $'\t' \
		> ${sample}.mapped_reads.txt
done

while read -r sample taxid fasta
do
	grep "^>" $fasta \
		| cut -c 2- \
		> mapped_reads/${sample}.taxid_${taxid}.reads.txt
	
	awk -v OFS="\t" -v taxid="$taxid" \
		'{print $1, taxid}' \
		mapped_reads/${sample}.taxid_${taxid}.reads.txt \
		>> ${sample}.mapped_reads.txt
done < query_list.txt
echo "Done"
