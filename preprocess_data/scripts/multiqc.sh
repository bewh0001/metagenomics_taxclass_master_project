#!/bin/bash

set -ueo

HERE="$(dirname $0)"
INDIR=${HERE}/fastqc

echo "Compiling MultiQC report from directory $INDIR"
mkdir -p multiqc
multiqc -o multiqc $INDIR
echo "Done"
