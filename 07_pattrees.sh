#!/bin/bash

OUTDIR=07_pattrees
/usr/bin/env Rscript 07_pattrees.r
for PDF in $(ls $OUTDIR); do
    pdfjam --quiet --nup 2x1 $OUTDIR/$PDF --outfile $OUTDIR/$PDF
    pdfcrop $OUTDIR/$PDF $OUTDIR/$PDF
done
pdfunite $OUTDIR/*.pdf 07_pattrees.pdf
