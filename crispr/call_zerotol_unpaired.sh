#!/bin/bash

usage="usage: call_zerotol_paired.sh  file1_R1_001.fastq.gz  [file2_R1_001.fastq.gz ...]"
[[ -f $1 ]] || { echo $usage; exit 0; }

for gzfile in $@; do
    if [ ! -f $gzfile ]; then
	echo "ERROR: Cannot find $gzfile - skip"
	continue
    fi
    if [ ! "${gzfile:(-9)}" == ".fastq.gz" ]; then
	echo "ERROR: File $gzfile does not have extension .fastq.gz - skip"
	continue
    fi
    if [ ! "${gzfile:(-16)}" == "_R1_001.fastq.gz" ]; then
	echo "ERROR: File $gzfile does not seem to be forward read '*_R1_001.fastq.gz' - skip"
	continue
    fi
    
    echo "" >&2
    gzfile_dir=$(dirname $gzfile)
    gzfile_bn=$(basename $gzfile)
    fileid=${gzfile_bn%_R1_001.fastq.gz}
    r1_gzfile="${gzfile_dir}/${fileid}_R1_001.fastq.gz"
    # r2_gzfile="${gzfile_dir}/${fileid}_R2_001.fastq.gz"
    # if [ ! -f $r2_gzfile ]; then
    # 	echo "ERROR: Cannot find reverse read fastq '$r2_gzfile' - skip"
    # 	continue
    # fi

    outfile="${fileid}.out"
    date > $outfile
    
    r1_reads=$(gunzip -c $r1_gzfile | wc -l | awk '{print $1/4}')
    
    # Align FASTQ to reference
    if [ -f ${fileid}.txt ]; then
	echo "Using excisting file ${fileid}.txt" >> $outfile
    else
	echo "=== Remove adapters $fileid"  >> $outfile
	# 5' adaptet: TTGTGGAAAGGACGAAACACCG  --  20 baser --  3' adapter: GTTTTAGAGCTAGAAATAGCAAGTAGATCGGA
        # Jeg fjerner de 3 første og 9 sidste baser i ovenstående fordi der typisk er flere fejl her
        # --match-read-wildcards : N matcher alt, reads ser ud til at have ret mange N'er? Hvis ellers der er et unikt match på insert synes jeg det er et fint match
        # -e 0.2 : Max fejlrate 20% dvs med 18 og 23 baser i adapter er op til 3 og 4 mismatch tilladt
        # -O 10 : Mindst 10 baseres overlap en adapter er identificeret 
        # -m 20 -M 20 : Overvej at fjerne reads hvor insertet ikke er 20 baser. Det har nok mest at gøre med fil størrelser og frekvenser (måske)
	#  --discard-untrimmed
	cutadapt -e 0.2 -O 10 --match-read-wildcards -g TGGAAAGGACGAAACACCG...GTTTTAGAGCTAGAAATAGCAAG -m 20 -M 20 -o ${fileid}.fastq  $r1_gzfile  >> $outfile

	echo "=== Re-format FATSQ to one-read-per-line format $fileid"  >> $outfile
	awk 'BEGIN { id="NA"; seq="NA"; qual="NA"}; 
                   { if(NR%4==1) {id=$1} 
                     else if(NR%4==2) {seq=$1} 
                     else if (NR%4==3) { if($1!="+") {printf("ERROR: Bad format at line %d\n",NF)>"/dev/stderr"; exit(2)} } 
                     else { qual=$1; printf("%-45s  %18s  %18s\n",id,seq,qual); id="NA"; seq="NA"; qual="NA" }
                   }' ${fileid}.fastq > ${fileid}.txt
	rm ${fileid}.fastq
    fi

    # Count aligned peptides
    n_filtered=$(cat ${fileid}.txt | wc -l)
    echo "=== Counting variants from ${fileid}.txt"  >> $outfile
    awk '{print $2}' ${fileid}.txt | sort | uniq -c \
	| awk -v tot_counts="$n_filtered" '{ print sprintf("%-10s %4d %7.2f", $2, $1, $1*1.0e6/tot_counts) }' > ${fileid}_counts.txt
    n_unq=$(cat ${fileid}_counts.txt | wc -l)

    # # Only keep peptides with a significant number of observations
    # min_obs=2
    # awk -v min_obs=$min_obs '{if ($2 >= min_obs) {print $0}}' ${fileid}_counts.txt > ${fileid}_counts_sig.txt
    # n_minobs=$(cat ${fileid}_counts_sig.txt | wc -l)
    
    echo "" >> $outfile
    echo "=== Done $fileid" >> $outfile
    echo "    Found $n_unq unique variants from $n_filtered length-filtered reads ($r1_reads R1 reads)" >> $outfile
    echo "" >> $outfile
    
    # echo "$fileid $r1_reads $n_filtered $n_unq $n_minobs"
    counted_pct=$(((n_filtered*100)/r1_reads))
    printf "LABELS: %30s  %9s  %9s  %7s  %3s\n" "file"  "reads1"  "filtered"  "unique" "pct"
    printf "RESULT: %30s  %9d  %9d  %7d  %3d\n" $fileid  $r1_reads  $n_filtered  $n_unq $counted_pct

    # [[ $r1_reads -eq $r2_reads ]] || { echo "WARNING: Different number of reads in R1 and R2"; }
    [[ $counted_pct -ge 80 ]] || { echo "WARNING: Final counts contain less than 80% of reads"; }

    gzip ${fileid}.txt
    # tar -zcf ${fileid}.tgz ${fileid}.txt ${fileid}.un1 ${fileid}.un2
    # rm ${fileid}.txt ${fileid}.un1 ${fileid}.un2
    
    # # Convergence profile
    # echo "=== Building convergence profile from ${fileid}_aligned.sam" >&2
    # echo "0  0" > ${fileid}_profile.txt
    # # echo "0  0  0" > ${fileid}_profile.txt
    # n_points=30
    # incr=$((n_aligned / (n_points*100) * 100)) # increment rounded to whole 100
    # for n in $(seq -f %.0f $incr $incr $n_aligned); do
    # 	echo "   count  $n  of  $n_aligned"  >&2
    # 	n_unq_h=$(head -n $n ${fileid}_aligned.sam  |  awk '{print $3}'  |  sort  |  uniq -c  |  wc -l)
    # 	# Could also use sample to do random selections but 1 test shows that it is not worth doing this
    # 	# n_unq_t=$(tail -n $n ${fileid}_aligned.sam  |  awk '{print $3}'  |  sort  |  uniq -c  |  wc -l)
    # 	echo "$n  $n_unq_h  $n_unq_t" >> ${fileid}_profile.txt
    # done
    # echo "$n_aligned  $n_unq" >> ${fileid}_profile.txt
    # # echo "$n_aligned  $n_unq  $n_unq" >> ${fileid}_profile.txt
done
