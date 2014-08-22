#!/bin/bash

INDIR=04_fasta
OUTDIR=05_align
mkdir -p $OUTDIR
for FASTA in $(ls -1 04_fasta | grep fasta); do
    if [[ ! -f $OUTDIR/$FASTA ]]; then
        mafft --localpair --maxiterate 1000 $INDIR/$FASTA > $OUTDIR/$FASTA
    else
        echo "Outfile $OUTDIR/$FASTA already exists"
    fi
done
