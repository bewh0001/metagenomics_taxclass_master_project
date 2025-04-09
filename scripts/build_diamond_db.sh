#!/bin/bash
#SBATCH -A sens2022554
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 10-00:00:00
#SBATCH -C mem256GB

set -ueo
source script_includes.sh

echo "Loading modules"
module load bioinfo-tools
module load SeqKit
echo "Done"

echo "Loading conda environment"
source ./load_conda.sh
echo "Done"

while test $# -gt 0
do
        case "$1" in
                # database name
                --dbname)
                        shift
                        if [ $# -gt 0 ]; then
                                DBNAME=$1
                        fi
                ;;
                #  database output directory
                --dbdir)
                        shift
                        if [ $# -gt 0 ]; then
                                DBDIR=$(abs_path $1)
                        fi
                ;;
                # work directory (optional)
                --work)
                        shift
                        if [ $# -gt 0 ]; then
                                WORK=$(abs_path $1)
                        fi
                ;;
                # directory with NCBI taxonomy files
                --taxonomy)
                        shift
                        if [ $# -gt 0 ]; then
                                TAXONOMY=$(abs_path $1)
                        fi
                ;;
                # Accession to taxid map file
                --accession2taxid)
                        shift
                        if [ $# -gt 0 ]; then
                                TAXMAP=$(abs_path $1)
                        fi
                ;;
                # Taxprofiler style csv samplesheet with fields
                # id,taxid,fasta_dna,fasta_aa
                --samplesheet)
                        shift
                        if [ $# -gt 0 ]; then
                                SAMPLESHEET=$(abs_path $1)
                        fi
                ;;
                *) echo "unknown argument $1"
                ;;
        esac
        shift
done

mkdir -p $WORK
cd $WORK

echo "Building DIAMOND protein library from samplesheet $SAMPLESHEET"
        LIBRARY="${PWD}/library/protein/protein.library.faa"
        mkdir -p library
        cut -d "," -f 4 $SAMPLESHEET \
                | tail -n +2 \
                | xargs cat \
                        > ${LIBRARY}.tmp
        echo "Done"

        # Remove duplicate sequences using SeqKit
        echo "Deduplicating DIAMOND library"
        seqkit rmdup ${LIBRARY}.tmp \
                > "$LIBRARY"
        rm ${LIBRARY}.tmp
echo "Done"

echo "Building diamond database"
mkdir -p db
cd db

diamond makedb \
        --in ${LIBRARY} \
        --db ${DBNAME} \
        --taxonmap ${TAXMAP} \
        --taxonnodes ${TAXONOMY}/nodes.dmp \
        --taxonnames ${TAXONOMY}/names.dmp \
        --threads ${CPUS}
cd ..

echo "Writing database to $DBDIR..."
if [ -d $DBDIR ]; then
        echo "Output database directory exists and will be overwritten"
        rm -r $DBDIR
fi

mkdir -p $DBDIR
mv db/* $DBDIR/
echo "Done"

echo "Cleaning up temporary files"
        rm -r library
echo "Done"
echo "Finished"
