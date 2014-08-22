#!/bin/bash

THREADS=6
INDIR=05_align
OUTDIR=06_tree

for FASTA in $(ls $INDIR | grep fasta); do
  ID=$(echo $FASTA | cut -d '.' -f 1)
  if [[ ! -f $OUTDIR/$ID.nwk ]]; then
    raxmlHPC-PTHREADS-AVX -T $THREADS -m GTRGAMMA -n $ID -s $INDIR/$FASTA -p 1 -w $(pwd)/$OUTDIR
    cp $OUTDIR/RAxML_bestTree.$ID $OUTDIR/$ID.nwk
    #FastTreeMP -nt -gtr -nosupport < $FASTA > $(pwd)/5_tree/$ID.nwk
  else
    echo "Outfile $OUTDIR/$ID.nwk already exists"
  fi
done
