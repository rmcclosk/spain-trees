#!/bin/sh

03_figures/optimize.py | tail -n 4 > 03_figures/params.r
python 03_figures/parse_spain.py > 03_figures/parse_spain.csv
03_figures/frequencyByDay.R
03_figures/fitnessCurve.R
pdfcrop 03_figures/legend.pdf 03_figures/legend.pdf
pdfcrop 03_figures/legend_detailed.pdf 03_figures/legend_detailed.pdf
