options(width=160, digits=4, stringsAsFactors=F)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    samples_file = "samples.csv"
    files = c("T6-5Lower_S4_L001_counts.txt","T6-5Lower_S4_L002_counts.txt","T6-5Lower_S4_L003_counts.txt")
} else if (length(args) < 2) {
    print("")
    print("usage: Rscript merge_and_map.r  <samples.csv>  <counts1.txt>  [counts2.txt  ...]")
    quit(save="no")
} else {
    samples_file = args[1]
    files = args[2:length(args)]
}
print(sprintf("Args: %s",paste0(args, collapse=" ")))


##
## Settings
##
settings = list()

# number of read counts required to consider DNA covered in summary table
settings$coverage_min_counts = 20

# store unmapped reads for later analysis
settings$store_unmapped = TRUE
settings$unmapped_min_reads = 5
settings$unmapped_good_lengths = c(20)


###
###  Libraries
###
# Read library files
# These are made from which(lib$even1), which(lib$even2), etc and should be non-redundant and including controls
raw = read.csv(gzfile("tkov3_guide_sequence.csv.gz"))
stopifnot(all(colnames(raw) == c("GENE","SEQUENCE","GUIDE_ID","TARGET.EXON")))
name_cols = c("gene","dna","GUIDE_ID","TARGET.EXON")
colnames(raw) = name_cols


###
###  Control sequences
###
i_ctrl = which(raw$TARGET.EXON=="CTRL")
ctrl = data.frame(gene = raw[i_ctrl,"gene"], dna = raw[i_ctrl,"dna"])

samples = read.csv(samples_file)

###
### Read count files
###
raw_read_counts = list()
if (settings$store_unmapped) { unmapped_raw = list() }

# for each file, map counts to library and put them in a data frame
for (file in files) {
    # check filename and extract file_id
    fl = strsplit(file, "/")[[1]]
    stopifnot(substr(fl[length(fl)], nchar(fl[length(fl)])-10, nchar(fl[length(fl)])) == "_counts.txt")
    file_id = substr(fl[length(fl)], 1, nchar(fl[length(fl)])-11)

    # determine sample index
    if ( grepl("_L00",file_id) ) {
        si = which(samples$file == substr(file_id, 1, nchar(file_id)-5))
    } else {
	si = which(samples$file == file_id)
    }
    stopifnot(length(si) == 1)

    # read counts file: col 1 should be the DNA sequence, col 1 the counts
    cf = read.table(file)
    colnames(cf) = c("dna","counts","rpm")

    # map file read counts to library
    i_mapped = match(raw[,"dna"],cf$dna)
    raw[,file_id] = cf[i_mapped,"counts"]
    
    # set un-observed variants to zero counts
    i_na = is.na( raw[,file_id] )
    raw[i_na,file_id] = 0

    # store number of unique variants and total read counts
    n_unq_dna = nrow(cf)
    n_raw_counts = sum(cf$counts)
    raw_read_counts[[file_id]] = c(n_unq_dna, n_raw_counts)

    # Report
    n_mapped_dna = sum(raw[,file_id] > 0)
    n_mapped_counts = sum(raw[,file_id])
    # n_filtered_counts = sum(cf$counts)
    # print(sprintf("Mapped %d dna seq covering %.2f%% of library and %d of %d length-filtered read counts (%.2f%%)",
    #               n_mapped_dna, n_mapped_dna/nrow(raw)*100,
    # 		  n_mapped_counts, n_filtered_counts, n_mapped_counts/n_filtered_counts*100))

    print(sprintf("Mapped %d out of %d raw counts (%.2f%%) from %s",n_mapped_counts,n_raw_counts,n_mapped_counts/n_raw_counts*100,file_id))

    # store unmapped redas for later analysis
    if (settings$store_unmapped) {
        i_above_threshold = which(cf$counts >= settings$unmapped_min_reads)
        unmapped_raw[[file_id]] = cf[setdiff(i_above_threshold,i_mapped),c("dna","counts")]
    }
    
    # print(table(unlist(strsplit(cf[,"dna"],""))))
}


# merge lane counts
# counts = list()
if (settings$store_unmapped) { unmapped_counts = list() }
lane_names = list()
samples$obs = FALSE
    
# init counts data frame
counts = raw[,name_cols]
    
# column names of columns with reads
cn = colnames(raw)[(length(name_cols)+1):ncol(raw)]
    
# column names with lane info removed
cn_nolane = unique( sapply(cn, function(s){ if (grepl("_L00",s)) {substr(s,1,nchar(s)-5)} else {s} }) )
    
for (sample_name in cn_nolane) {
    stopifnot(sample_name %in% samples$file)
    si = which(samples$file==sample_name)
    stopifnot(! samples[si,"obs"])
    samples[si,"obs"] = TRUE
	
    sample_lane_names = cn[which(grepl(sample_name, cn))]
    counts[,sample_name] = apply(raw[,sample_lane_names], MARGIN=1, sum)

    # consider merging raw_read_counts[[file_id]] = c(n_unq_dna, n_raw_counts)

    # report correlations of raw[,sample_lane_names]
    lane_names[[sample_name]] = sample_lane_names
    print(sprintf("Merging %d lanes into %s", length(sample_lane_names), sample_name))

    n = length(sample_lane_names)
    if (n > 1) {
        rps = c()
        for (i in seq(0,n-1)) for (j in seq(i,n-1)) if (i!=j) {
            rp = cor(raw[,sample_lane_names[i+1]], raw[,sample_lane_names[j+1]], method="pearson", use="complete.obs")
	    rps = c(rps,rp)
        }
	print(sprintf("  average pearson %.3f, range %.3f-%.3f", mean(rps), min(rps), max(rps)))
    }

    # merge lanes for unmapped read counts
    if (settings$store_unmapped) {
        sample_unmapped_dna = unique(unlist(sapply(unmapped_raw[sample_lane_names], "[", "dna")))

        unmapped_counts[[sample_name]] = data.frame(dna = sample_unmapped_dna)
        for (sn in sample_lane_names) {
        unmapped_counts[[sample_name]][,sn] = unmapped_raw[[sn]][match(sample_unmapped_dna,unmapped_raw[[sn]][,"dna"]),"counts"]
	}
	# print(sprintf("Merged lanes for %s resulting in %d unmapped reads from %d unique DNA sequences",
	#               sample_name, sum(unmapped_counts[[sample_name]][,"sum"]), length(sample_unmapped_dna)))
    }
}
# print(sprintf("Done lane merging, distribution of lanes per sample:"))
# print(table(unlist(n_lanes)))

# counts per sample, info on lanes and unmapped reads not present
save(counts, samples, name_cols, ctrl, settings, file="counts.rda")

# data on unmapped reads and lanes
if (settings$store_unmapped) {
    save(unmapped_counts, raw_read_counts, lane_names, file="counts_unmapped.rda")
}

# table of read counts
# all_lib_dna = unique(unlist(sapply(raw, "[", "dna")))
counts_summary = data.frame(sample=NULL, coverage = NULL, ctrl_cover=NULL,
                            all=NULL, mapped=NULL, controls=NULL, unmapped=NULL, unmapped_bad_len=NULL)

for (sn in sort(names(unmapped_counts))) {
    si = which(sn == samples$file)
    # lib = samples[si,"lib"]
    sample_lane_names = lane_names[[sn]]

    # i_ctrl = which(counts[[lib]][,"dna"] %in% ctrl$dna) 
    # i_ctrl = which(counts[[lib]][,"name"] %in% ctrl$name)
    stopifnot(length(i_ctrl) == nrow(ctrl))

    mapped_other_lib = 0
    unmapped = 0
    unmapped_bad_len = 0

    if (settings$store_unmapped) {
        unmapped_counts[[sn]][,"sum"] = apply(unmapped_counts[[sn]][,sample_lane_names], MARGIN=1, sum, na.rm=T)
	
        unmapped = sum(unmapped_counts[[sn]][,"sum"])

        # unmapped DNA with bad length, i.e. length not in
        unmapped_counts[[sn]][,"length"] = sapply(unmapped_counts[[sn]][,"dna"], nchar)
        ibl = which(! unmapped_counts[[sn]][,"length"] %in% settings$unmapped_good_lengths)
        unmapped_bad_len = sum(unmapped_counts[[sn]][ibl,"sum"])
    }
    
    df = data.frame(sample = sn,
                    coverage = sum(counts[,sn] >= settings$coverage_min_counts)/nrow(counts)*100,
		    ctrl_cover = sum(counts[i_ctrl,sn] >= settings$coverage_min_counts)/nrow(ctrl)*100,
                    all = sum(sapply(raw_read_counts[sample_lane_names], "[", 2)),
                    mapped = sum(counts[,sn]),
	            controls = sum(counts[i_ctrl,sn]),
#	            mapped_other_lib = mapped_other_lib,
		    unmapped = unmapped,
		    unmapped_bad_len = unmapped_bad_len)
    counts_summary = rbind(counts_summary, df)
}
print("")
print("Read counts summary")
print(counts_summary)
write.csv(counts_summary, file="counts_summary.csv")

pct_summary = counts_summary
cns = colnames(counts_summary)[5:ncol(counts_summary)]
for (ri in seq(nrow(pct_summary))) {
    pct_summary[ri,cns] = pct_summary[ri,cns]/pct_summary[ri,"all"]*100
}
print("")
print("Read percent summary")
print(pct_summary)
