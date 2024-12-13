module load eautils

ls ./fastq/*R1_001.fastq.gz  |  parallel -P 16 ./call_zerotol_unpaired.sh {};

grep LABEL run.out | head -n1 > paring_stat.txt
grep RESULT run.out >> paring_stat.txt

module load gcc/11.2.0
module load R
Rscript merge_and_map.r samples.csv *_counts.txt > merge_and_map.out

