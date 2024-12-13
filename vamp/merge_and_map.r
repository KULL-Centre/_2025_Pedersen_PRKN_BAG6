options(width=160)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    input_file = "samples.csv"
} else if (length(args) < 1) {
    print("")
    print("usage: Rscript merge_and_map.r  <samples.csv>")
    quit(save="no")
} else {
    input_file = args[1]
}
print(sprintf("Args: %s",paste0(args, collapse=" ")))


##
## Settings (default)
##
settings = list()

# number of read counts required to consider DNA covered in summary table
settings$coverage_min_counts = 20

# store unmapped reads for later analysis
settings$store_unmapped = TRUE
settings$include_unmapped = FALSE
settings$unmapped_min_reads = 1

settings$save_csv = TRUE
settings$auto_detect_controls = FALSE

settings$library_is_dna = FALSE
settings$library_is_protein = FALSE
settings$library_is_barcode = FALSE

settings$allow_truncated_translation = FALSE

# map the class of settings for casting when string have been read from input file
bool_settings = c()
num_settings = c()
for (setting in names(settings)) {
    if (class(settings[[setting]]) == "logical") { bool_settings = c(bool_settings,setting) }
    if (class(settings[[setting]]) == "numeric") { num_settings = c(num_settings,setting) }
}
    

##
## Read input file, may contain settings
##

# read input file
max_col = max(count.fields(input_file, sep=",", comment.char="")) 
if (max_col < 5) {
    print(sprintf("ERROR: Input file '%s' does ont seem to contain any files (Too few columns %d < 5)",input_file,max_col))
    stopifnot(max_col>=5)
}
input = read.csv(input_file, sep=",", header=F, fill=T, strip.white=T, blank.lines.skip=F, comment.char="",
                 col.names=sprintf("V%d", seq(max_col)))

# variables to store info from input
library_files = list()
samples = list()

# read input file row by row
reading_section = "none"
for (ri in seq(nrow(input))) {
    # skip blank lines and commont lines
    if (all(input[ri,] == "")) { next }
    if (substr(input[ri,1],1,1)=="#") { next }
    # if (any(substr(input[ri,],1,1)=="#")) { next }

    # store which section is being read
    if (toupper(input[ri,1]) %in% c("SETTINGS","LIBRARIES","SAMPLES")) {
        reading_section=tolower(input[ri,1])
	next
    }

    # store info dependen on current section
    if (reading_section=="settings") {
         setting = tolower(input[ri,1])
         if (! setting %in% names(settings)) {
	     print(sprintf("ERROR: Setting '%s' is unknown",setting))
	     stopifnot(F)
	 }
	 settings[[setting]] = input[ri,2]
    } else if (reading_section=="libraries") {
        library = tolower(input[ri,1])
	if (library %in% names(library_files)) { print(sprintf("ERROR: Library '%s' already listeds",library)); stopifnot(F) }
	if (! file.exists(input[ri,2])) {
	    print(sprintf("ERROR: Cannot find library file '%s' listed in row %d", input[ri,2],ri))
	    stopifnot(file.exists(input[ri,2]))
	}
        library_files[[library]] = input[ri,2]
    } else if (reading_section=="samples") {
        library = gsub(" ", "_", tolower(input[ri,1]))
        selection = gsub(" ", "_", tolower(input[ri,2]))
        replicate = gsub(" ", "_", tolower(input[ri,3]))
	# Not used yet, and available in the variable named 'input'
        # description = input[ri,4]
	sample_list = c()
	for (ci in seq(5,ncol(input))) {
	    if (input[ri,ci] != "") {
	    	if (! file.exists(input[ri,ci])) {
		    print(sprintf("ERROR: Cannot find file '%s' listed in row %d column %d under library %s, selection %s, replica %s",
		                  input[ri,ci],ri,ci,library,selection,replicate))
		    stopifnot(file.exists(input[ri,ci]))
		}
	        # sample_list = c(input[ri,ci],sample_list)
	        sample_list = c(sample_list, input[ri,ci])
	    }
	}
	# make a word that combines selection and replicate separated by dash - need to remove dash from names first
	sel_rep = paste(gsub("-","_",selection), gsub("-","_",replicate), sep="-")
	# stopifnot(sel_rep != "-")
	if (! library %in% names(samples)) { samples[[library]] = list() }
	if (sel_rep %in% names(samples[[library]])) { print(sprintf("ERROR: Sample '%s' already assigned",sel_rep)); stopifnot(F) }
	if (length(sample_list) == 0) {
	    print(sprintf("ERROR: No files listed in row %d with library '%s' selection '%s' replica '%s'",ri,library,selection,replicate))
	    stopifnot(length(sample_list)>0)
	}
	samples[[library]][[sel_rep]] = sample_list
    }    
}

# Cast and check settings
for (setting in bool_settings) { settings[[setting]] = as.logical(settings[[setting]]) }
for (setting in num_settings) { settings[[setting]] = as.numeric(settings[[setting]]) }

# The library sequences should be either protein coding, barcodes or any other DNA
stopifnot(settings$library_is_dna +
          settings$library_is_protein +
          settings$library_is_barcode == 1)

if (settings$include_unmapped & ! settings$store_unmapped) {
    print("Setting store_unmapped to TRUE in order to include unmapped (setting include_unmapped)")
    settings$store_unmapped = TRUE
}


##
## Helper functions
##

codon_table = list()
# T          T                      C                    A                     G
codon_table[["TTT"]] = "F"; codon_table[["TCT"]] = "S"; codon_table[["TAT"]] = "Y"; codon_table[["TGT"]] = "C" # T
codon_table[["TTC"]] = "F"; codon_table[["TCC"]] = "S"; codon_table[["TAC"]] = "Y"; codon_table[["TGC"]] = "C" # C
codon_table[["TTA"]] = "L"; codon_table[["TCA"]] = "S"; codon_table[["TAA"]] = "*"; codon_table[["TGA"]] = "*" # A
codon_table[["TTG"]] = "L"; codon_table[["TCG"]] = "S"; codon_table[["TAG"]] = "*"; codon_table[["TGG"]] = "W" # G
# C
codon_table[["CTT"]] = "L"; codon_table[["CCT"]] = "P"; codon_table[["CAT"]] = "H"; codon_table[["CGT"]] = "R" # T
codon_table[["CTC"]] = "L"; codon_table[["CCC"]] = "P"; codon_table[["CAC"]] = "H"; codon_table[["CGC"]] = "R" # C
codon_table[["CTA"]] = "L"; codon_table[["CCA"]] = "P"; codon_table[["CAA"]] = "Q"; codon_table[["CGA"]] = "R" # A
codon_table[["CTG"]] = "L"; codon_table[["CCG"]] = "P"; codon_table[["CAG"]] = "Q"; codon_table[["CGG"]] = "R" # G
# A
codon_table[["ATT"]] = "I"; codon_table[["ACT"]] = "T"; codon_table[["AAT"]] = "N"; codon_table[["AGT"]] = "S" # A
codon_table[["ATC"]] = "I"; codon_table[["ACC"]] = "T"; codon_table[["AAC"]] = "N"; codon_table[["AGC"]] = "S" # C
codon_table[["ATA"]] = "I"; codon_table[["ACA"]] = "T"; codon_table[["AAA"]] = "K"; codon_table[["AGA"]] = "R" # A
codon_table[["ATG"]] = "M"; codon_table[["ACG"]] = "T"; codon_table[["AAG"]] = "K"; codon_table[["AGG"]] = "R" # G
# G
codon_table[["GTT"]] = "V"; codon_table[["GCT"]] = "A"; codon_table[["GAT"]] = "D"; codon_table[["GGT"]] = "G" # A
codon_table[["GTC"]] = "V"; codon_table[["GCC"]] = "A"; codon_table[["GAC"]] = "D"; codon_table[["GGC"]] = "G" # C
codon_table[["GTA"]] = "V"; codon_table[["GCA"]] = "A"; codon_table[["GAA"]] = "E"; codon_table[["GGA"]] = "G" # A
codon_table[["GTG"]] = "V"; codon_table[["GCG"]] = "A"; codon_table[["GAG"]] = "E"; codon_table[["GGG"]] = "G" # G

translate = function(dna, allow_truncation=FALSE, non_nat_char="x") {
    n = nchar(dna)
    codons = substring(dna, seq(1, n-2, by=3), seq(3, n, by=3))
    if (length(codons)*3 != n & ! allow_truncation) return(NA)
    aas = ifelse(codons %in% names(codon_table), codon_table[codons], non_nat_char)
    paste0(aas, collapse="")
}
# translate = function(dna, allow_truncation=FALSE) {
#     n = nchar(dna)
#     codons = substring(dna, seq(1, n-2, by=3), seq(3, n, by=3))
#     if (length(codons)*3 != n & ! allow_truncation) return(NA)
#     paste0(codon_table[codons], collapse="")
# }


# read a sequence file (.seq) - a text file with one sequence per line
# return a data.frame with 2 columns: name and dna
read.seqfile = function(filename) {
    df = read.table(filename)
    if (ncol(df) == 1) {
        colnames(df) = c("dna")
    } else if (ncol(df) == 2) {
        colnames(df) = c("name","dna")
    } else {
        print(sprintf("ERROR: Seq-file %s does not have 1 or 2 columns", filename))
	return(NULL)
    }
    return(df)
}

# read a csv file with sequences
# return a data.frame with 2 columns: name and dna
read.seqcsv = function(filename, name_colname="name", dna_colname="dna", other_colnames=c()) {
    df = read.csv2(filename)
    if (! name_colname %in% colnames(df)) {
        print(sprintf("ERROR: Cannot find name column '%s' in %s",name_colname,filename))
	return(NA)
    }
    if (! dna_colname %in% colnames(df)) {
        print(sprintf("ERROR: Cannot find DNA sequence column '%s' in %s",dna_colname,filename))
	return(NA)
    }
    ret_df = data.frame(name=df[,name_colname], dna=df[,dna_colname])

    for (cn in other_colnames) {
        stopifnot(cn %in% colnames(df))
        ret_df[,cn] = df[,cn]
    }
    return(ret_df)
}

###
###  Libraries
###

# based on library type, 
if (settings$library_is_protein) {
    lib_reader = function(file) { read.seqcsv(file, dna_colname="seq") }
} else if (settings$library_is_barcode) {
    lib_reader = function(file) { read.seqcsv(file, dna_colname="barcode", other_colnames=c("subst_aa","subst_dna","n_subst_aa","n_subst_dna")) }
} else if (settings$library_is_dna) {
    lib_reader = function(file) { read.seqcsv(file) }
} else { stopifnot(F) }

# Read library files
# These are made from which(lib$even1), which(lib$even2), etc and should be non-redundant and including controls
counts = list()
for (lib in names(library_files)) {
    if (lib=="controls") { next }
    counts[[lib]] = lib_reader(library_files[[lib]])
    
    # DNA names and sequences should be unique per library
    stopifnot(length(counts[[lib]][,"name"]) == length(unique(counts[[lib]][,"name"])))    
    stopifnot(length(counts[[lib]][,"dna"]) == length(unique(counts[[lib]][,"dna"])))    
}

# calculate Pearson correlations between all pairs of columns in a data frame
pair_cor = function(df) {
    n = ncol(df)
    rps = c()
    for (i in seq(0,n-1)) for (j in seq(i,n-1)) if (i!=j) {
        rp = cor(df[,i+1], df[,j+1], method="pearson")
	rps = c(rps,rp)
    }
    return(rps)
}

plot_cor = function(df, pch=".", ...) {
    cns = colnames(df)
    if (length(cns) == 2) {
	plot(df, xlab=cns[1], ylab=cns[2], pch=pch, ...)
    } else if (length(cns) > 2) {
        plot(df, labels=cns, upper.panel=NULL, pch=pch, ...)
    }
}


###
###  Control sequences
###

# Mostly relevant with more libraries where controls should be present in all libraries
ctrl = data.frame()

# If given, read controls from csv and add to all libraries
if ("controls" %in% library_files) {
    ctrl = read.csv(library_files[["controls"]])

    # expand control df with columns from libraries before merging
    ctrl_merge = ctrl
    for (cn in colnames(counts[[1]])) {
        if (! cn %in% c("name","dna")) {
	    ctrl_merge[,cn] = NA
	}
    }

    # add controls to all libraries
    for (lib in names(counts)) { counts[[lib]] = rbind(counts[[lib]],ctrl_merge) }

    print(sprintf("Added %d control sequences from %s top all libraries",nrow(ctrl),library_files[["controls"]]))
    
} else if (length(counts) > 1 & settings$auto_detect_controls) {
    # autodetect controls as sequences common to all libraries (if there are)
    ctrl_dna = counts[[1]][,"dna"]
    for (ri in 2:length(counts)) { ctrl_dna = intersect(ctrl_dna, counts[[ri]][,"dna"]) }
    # get names from first library
    ctrl = data.frame(names = counts[[1]][match(ctrl_dna,counts[[1]][,"dna"]),"name"], dna = ctrl_dna)
    print(sprintf("Using %d sequences common to all libraries as controls",nrow(ctrl)))
    
} else {
    print("Experiment has no controls!")
}


# translate coding seqiuences
if (settings$library_is_protein) {
    for (lib in names(counts)) {
        counts[[lib]][,"aa"] = sapply(counts[[lib]][,"dna"], translate, allow_truncation=settings$allow_truncated_translation)
        print(sprintf("Check if library %s contains unique AA sequences: %s",lib,length(counts[[lib]][,"aa"]) == length(unique(counts[[lib]][,"aa"]))))
    }
    ctrl$aa = sapply(ctrl$dna, translate, allow_truncation=settings$allow_truncated_translation)
    print(sprintf("Check if controls contains unique AA sequences: %s", length(ctrl[,"aa"]) == length(unique(ctrl[,"aa"]))))
# } else if (settings$library_is_barcode) {
#     counts[[lib]][,"aa"] = 
}

# check that all library files have the same number of columns
name_cols = colnames(counts[[1]])
stopifnot(all(sapply(counts,ncol) == length(name_cols)))


###
### Read count files
###
raw_read_counts = list()
bio_rep = list()
file_ids = list()

# this can take substantial memory so only store unmapped sequences if requested
if (settings$store_unmapped) { unmapped_counts = list() }

for (lib in names(samples)) {
    bio_rep[[lib]] = list()
    file_ids[[lib]] = list()
    for (sel_rep in names(samples[[lib]])) {
        sel_rep_vec = strsplit(sel_rep,"-")[[1]]
        selection = sel_rep_vec[1]
        replicate = sel_rep_vec[2]
	if (selection %in% names(bio_rep[[lib]])) {
	    bio_rep[[lib]][[selection]] = c(bio_rep[[lib]][[selection]],replicate)
	} else {
            bio_rep[[lib]][[selection]] = c(replicate)
        }

	print(sprintf("Library %s selection '%s' replicate %s with %d files", lib, selection, replicate, length(samples[[lib]][[sel_rep]])))
	
        # collect technical replicates here
        tech_reps = data.frame(dna=counts[[lib]][,"dna"])
	n_total_reads = 0
        for (file in samples[[lib]][[sel_rep]]) {
            path = strsplit(file, "/")[[1]]
	    fn = path[length(path)]
	    file_id = NULL
	    if (substr(fn,nchar(fn)-3,nchar(fn)) == ".txt") {
	        file_id = substr(fn,1,nchar(fn)-4)
                cf = read.table(file)
	    } else if ( substr(fn,nchar(fn)-6,nchar(fn)) == ".txt.gz" ) {
	        file_id = substr(fn,1,nchar(fn)-7)
		cf = read.table(gzfile(file))
	    } else {
	        print(sprintf("Unkown extension of %s",file))
	    }
            stopifnot(! is.null(file_id))
	    
	    if (ncol(cf) == 3) {
	        # old read counts file format: col 1 should be the DNA sequence, col 2 the counts
 		colnames(cf) = c("dna","counts","rpm")
	    } else if (ncol(cf) == 2) {
	        # direct output from uniq -c: col 1 is counts, col 2 the sequence
		colnames(cf) = c("counts","dna")
	    } else {
	        print(sprintf("ERROR: Unknown count file format with %d columns: %s", ncol(cf), fn))
	    }
            # map file read counts to library
            i_mapped = match(tech_reps[,"dna"],cf$dna)
            tech_reps[,file_id] = cf[i_mapped,"counts"]
    
            # set un-observed variants to zero counts
            i_na = is.na( tech_reps[,file_id] )
            tech_reps[i_na,file_id] = 0

            # store number of unique variants and total read counts
            n_unq_dna = nrow(cf)
            n_raw_counts = sum(cf$counts)
	    n_total_reads = n_total_reads + n_raw_counts
            raw_read_counts[[file_id]] = c(n_unq_dna, n_raw_counts)

            # Report
            n_mapped_dna = sum(tech_reps[,file_id] > 0)
            n_mapped_counts = sum(tech_reps[,file_id])
	    print(sprintf("    File %s has %d read counts of %d unique sequences", file_id, n_raw_counts, n_unq_dna))
	    print(sprintf("      mapped %d read counts (%.1f%% of sequenced) to %d members of library %s (%.1f%% of sequenced and %.1f%% of library)",
	                  n_mapped_counts, n_mapped_counts/n_raw_counts*100, n_mapped_dna, lib, n_mapped_dna/nrow(cf)*100, n_mapped_dna/nrow(tech_reps)*100))
	    
            # store unmapped (but substantial) reads for later analysis
            if (settings$store_unmapped) {
                i_above_threshold = which(cf$counts >= settings$unmapped_min_reads)
                unmapped_counts[[file_id]] = cf[setdiff(i_above_threshold,i_mapped),c("dna","counts")]
		n_unmapped = sum(unmapped_counts[[file_id]][,"counts"])
		print(sprintf("      store %d unique unmapped reads each with %d or more counts, total %d (%.1f%%) counts",
		              nrow(unmapped_counts[[file_id]]), settings$unmapped_min_reads, n_unmapped, n_unmapped/n_raw_counts*100))
            }
        }
	file_ids[[lib]][[sel_rep]] = colnames(tech_reps)[2:ncol(tech_reps)]
	
	# sum counts over technical replicates
	counts[[lib]][,sel_rep] = apply(tech_reps[,2:ncol(tech_reps)], MARGIN=1, sum)

	n_mapped_reads = sum(counts[[lib]][,sel_rep])
	print(sprintf("    Total read counts %d, mapped %d (%.1f%%)",
	              n_total_reads, n_mapped_reads, n_mapped_reads/n_total_reads*100))
	print(sprintf("    Coverage: %d of %d (%.1f%%) library members have any read counts, %d (%.1f%%) have 10 or more and %d (%.1f%%) have 20 or more",
	              sum(counts[[lib]][,sel_rep]>0), nrow(counts[[lib]]), sum(counts[[lib]][,sel_rep]>0)/nrow(counts[[lib]])*100,
		      sum(counts[[lib]][,sel_rep]>=10), sum(counts[[lib]][,sel_rep]>=10)/nrow(counts[[lib]])*100,
		      sum(counts[[lib]][,sel_rep]>=20), sum(counts[[lib]][,sel_rep]>=20)/nrow(counts[[lib]])*100))

	# Pearson's correlations of all pairs, pseudocounts only for log correlations
	pseudocounts = 1
	rps = pair_cor(tech_reps[,2:ncol(tech_reps)])
	rps_log = pair_cor(log(tech_reps[,2:ncol(tech_reps)]+1))
	print(sprintf("    Avg Pearson %.2f range %.2f - %.2f  (log avg %.2f range %.2f - %.2f with pseudo counts %d)",
	              mean(rps),min(rps),max(rps),mean(rps_log),min(rps_log),max(rps_log),pseudocounts))

	# # plot correlations of technical replicates
	# quartz()
	# main = =sprintf("Lib %s, %s rep %s, avg Pearson %.2f (log %.2f)",lib,selection,replicate,mean(rps),mean(rps_log))
	# plot_cor(tech_reps[,2:ncol(tech_reps)]+pseudocounts, log="xy", main=main)
    }
    
    # extend library if unmapped reads should be included
    if (settings$store_unmapped & settings$include_unmapped) {
        new_dna = c()
        for (sel_rep in names(file_ids[[lib]])) {
	    for (file_id in file_ids[[lib]][[sel_rep]]) {
	        new_dna = unique(c(new_dna, unmapped_counts[[file_id]][,"dna"]))
	    }
	}
	print(sprintf("Extending library %s with %d unmapped reads (min read count %d per file)", lib, length(new_dna), settings$unmapped_min_reads))
	if (length(new_dna) >= 1e6) { um_names = sprintf("unmap%09d",seq_along(new_dna)) } else { um_names = sprintf("unmap%06d",seq_along(new_dna)) }
	um_df = data.frame(name = um_names, dna = new_dna)
        if (settings$library_is_protein) {
            um_df$aa = sapply(um_df$dna, translate, allow_truncation=settings$allow_truncated_translation)
	}
	rownames(um_df) = NULL		   
        for (sel_rep in names(file_ids[[lib]])) {
	    um_df[,sel_rep] = 0
            for (file_id in file_ids[[lib]][[sel_rep]]) {
	        umi = match(unmapped_counts[[file_id]][,"dna"], um_df$dna)
		stopifnot(all( ! is.na(umi) ))
	        um_df[umi,sel_rep] = um_df[umi,sel_rep] + unmapped_counts[[file_id]][,"counts"]
	    }
	    print(sprintf("    Added %d read counts to %d unmapped reads for selection %s", sum(um_df[,sel_rep]), sum(um_df[,sel_rep] > 0), sel_rep))
	}
	stopifnot(all( colnames(um_df) == colnames(counts[[lib]]) ))
	counts[[lib]] = rbind(counts[[lib]], um_df)
    }
}


# Analyses of biological replicates
for (lib in names(bio_rep)) {
    for (selection in names(bio_rep[[lib]])) {
        n_rep = length(bio_rep[[lib]][[selection]])
        print(sprintf("Library %s selection '%s' has %d biological replicates",lib,selection,n_rep))
        if (n_rep < 2) { next }

	# column names of biological replicates
        cns = paste(selection, bio_rep[[lib]][[selection]], sep="-")
        rps = pair_cor(counts[[lib]][,cns]) 
        rps_log = pair_cor( log(counts[[lib]][,cns]+pseudocounts) ) 
	print(sprintf("    Avg Pearson %.2f range %.2f - %.2f  (log avg %.2f range %.2f - %.2f with pseudo counts %d)",
	              mean(rps),min(rps),max(rps),mean(rps_log),min(rps_log),max(rps_log),pseudocounts))

        # for protein and barcodes, report correlations on AA level
        # if (settings$library_is_barcode | settings$library_is_protein) {
        if (settings$library_is_barcode) {
	    # all the following relis on having a column names subst containing colon-separated amino acid substitutions
	    stopifnot("subst_aa" %in% colnames(counts[[lib]]))
            # aggregate barcode counts per amino acid substitutions
            agg = aggregate(counts[[lib]][,cns], by=list(counts[[lib]][,"subst_aa"]), sum)
            aavar = data.frame(subst=agg[,1])
            aavar[,cns] = agg[,cns]
	    aavar$n_subst = sapply(strsplit(aavar$subst,":"), length)
	    i_wt = which(aavar$subst=="")
	    i_single = which(aavar$n_subst == 1)

            for (cn in cns) {
                print(sprintf("    Barcodes observed in bio. rep. %s for %d of %d (%.2f%%) amino acid variants in library %s",
	                      cn, sum(aavar[,cn]>0), nrow(aavar), sum(aavar[,cn]>0)/nrow(aavar)*100, lib))
  	        print(sprintf("        Avg read count %.1f, median %.1f, range %d - %d", mean(aavar[,cn]), median(aavar[,cn]), min(aavar[,cn]), max(aavar[,cn])))
                print(sprintf("        Barcodes observed for %d of %d (%.2f%%) single amino acid variants in library %s",
	                      sum(aavar[i_single,cn]>0), length(i_single), sum(aavar[i_single,cn]>0)/length(i_single)*100, lib))
  	        print(sprintf("        Single variant avg read count %.1f, median %.1f, range %d - %d",
		              mean(aavar[i_single,cn]), median(aavar[i_single,cn]), min(aavar[i_single,cn]), max(aavar[i_single,cn])))
	        if (length(i_wt) > 0) {
		    ic_wt = which(counts[[lib]][,"subst_aa"]=="")
		    n_wt_obs = sum(counts[[lib]][ic_wt,cn] > 0)
		    print(sprintf("        WT amino acid sequence read counts: %d (%.2f%%), single variants %d (%.2f%%)",
		                  aavar[i_wt,cn], aavar[i_wt,cn]/sum(aavar[cn])*100, sum(aavar[i_single,cn]), sum(aavar[i_single,cn])/sum(aavar[cn])*100))
		    print(sprintf("        WT barcodes observed %d of %d (%.2f%%)", n_wt_obs, length(ic_wt), n_wt_obs/length(ic_wt)*100))
		}
	    }

	    # report Pearson correlation on amino acid level without WT
	    if (length(i_wt) > 0) {
	        print("    Excluding WT from correlations at amino acid level")
	        rps = pair_cor(aavar[-i_wt,cns])
                rps_log = pair_cor( log(aavar[-i_wt,cns]+pseudocounts) )
	    } else {
	        rps = pair_cor(aavar[,cns])
                rps_log = pair_cor( log(aavar[,cns]+pseudocounts) )
	    }
            print(sprintf("    Amino acid level avg Pearson %.2f range %.2f - %.2f  (log avg %.2f range %.2f - %.2f with pseudo counts %d)",
                      mean(rps),min(rps),max(rps),mean(rps_log),min(rps_log),max(rps_log),pseudocounts))

            # report Pearson correlation for single amino acid variants
	    rps = pair_cor(aavar[i_single,cns])
            rps_log = pair_cor( log(aavar[i_single,cns]+pseudocounts) )
            print(sprintf("    Single AA variant avg AA Pearson %.2f range %.2f - %.2f  (log avg %.2f range %.2f - %.2f with pseudo counts %d)",
                      mean(rps),min(rps),max(rps),mean(rps_log),min(rps_log),max(rps_log),pseudocounts))

            # plot correlations of biological replicates on AA level
	    filename = sprintf("correlation_%s_%s.jpg", lib, selection)
	    main = sprintf("AA seq in lib %s, %s, avg Pearson %.2f (log %.2f)",lib,selection,mean(rps),mean(rps_log))
	    jpeg(filename, width=9, height=9, units="cm", res=300, pointsize=7.5)
	    if (length(i_wt) > 0) {
                plot_cor(aavar[-i_wt,cns]+pseudocounts, log="xy", main=main)
	    } else {
                plot_cor(aavar[,cns]+pseudocounts, log="xy", main=main)
	    }
	    dev.off()
	    
        } else {
	    # plot correlations of biological replicates on DNA level
	    filename = sprintf("correlation_%s_%s.jpg", lib, selection)
	    main = sprintf("DNA seq in lib %s, %s, avg Pearson %.2f (log %.2f)",lib,selection,mean(rps),mean(rps_log))
	    jpeg(filename, width=9, height=9, units="cm", res=300, pointsize=7.5)
	    plot_cor(counts[[lib]][,cns]+pseudocounts, log="xy", main=main)
	    dev.off()
	}
    }
}


# Analyses of unmapped reads here or in separate script - to come

print("")
print("Saving counts.rda")

# save
save(counts, bio_rep, input, samples, library_files, ctrl, settings, file="counts.rda")

# if requested, write a csv file per library
if (settings$save_csv) {
    print("Saving counts_[lib].csv")
    for (lib in names(counts)) {
        filename = sprintf("counts_%s.csv", lib)
        write.csv(counts[[lib]], file=filename, row.names=F, quote=F)
    }
}

# data on unmapped reads
if (settings$store_unmapped) {
    print("Saving unmapped counts to counts_unmapped.rda")
    save(unmapped_counts, raw_read_counts, file="counts_unmapped.rda")
}
