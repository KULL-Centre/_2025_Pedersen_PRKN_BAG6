#!/bin/bash
module load gcc/11.2.0
module load R

Rscript make_mageck_input.r

module load miniconda
# conda create -c bioconda -n mageck python=3.9 mageck
conda activate mageck

mageck test -k counts.txt -t T12.high -c T12 -n mageck_park_R42P_T12high
mageck test -k counts.txt -t T12.low -c T12 -n mageck_park_R42P_T12low
mageck test -k counts.txt -t T12 -c T12 -n mageck_park_R42P_T12
