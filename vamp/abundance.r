options(width=160)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    counts_file = "counts_prkn_bc.csv"
} else if (length(args) < 1) {
    print("")
    print("usage: Rscript abundance.r  <counts.csv>")
    quit(save="no")
} else {
    counts_file = args[1]
}
print(sprintf("Input counts file: %s",counts_file))

# Settings
settings = list(file=counts_file)

# Samples
settings$bio_reps = list()
settings$bio_reps[["ctrl"]] = sprintf("bio%d", seq(2))
settings$bio_reps[["b6ko"]] = sprintf("bio%d", seq(2))
settings$n_facs_bins = 2
settings$idx_facs_bins = seq(settings$n_facs_bins)
settings$pop_facs_bins = list()
settings$pop_facs_bins[["ctrl"]] = c(0.17, 0.83)
settings$pop_facs_bins[["b6ko"]] = c(0.05, 0.95)

# require this number of replicates to keep abundance score
settings$min_rep = 1

stopifnot( length(settings$idx_facs_bins) == settings$n_facs_bins )
stopifnot( length(settings$pop_facs_bins) == settings$n_facs_bins )

# names of coulmns with counts
replicas = c()
cn_counts = c()
for (cell in names(settings$bio_reps)) {
    for (bio_rep in settings$bio_reps[[cell]]) {
        replicas = c(replicas, sprintf("%s_%s",cell,bio_rep))
        for (facs_bin in seq(settings$n_facs_bins)) {
            cn_counts =c(cn_counts, sprintf("%s_bin%d.%s",cell,facs_bin,bio_rep))
        }
    }
}

# needed for making output files
wt = list()
wt[["name"]] = "PRKN"
wt[["aa"]] = "MIVFVRFNSSHGFPVEVDSDTSIFQLKEVVAKRQGVPADQLRVIFAGKELRNDWTVQNCDLDQQSIVHIVQRPWRKGQEMNATGGDDPRNAAGGCEREPQSLTRVDLSSSVLPGDSVGLAVILHTDSRKDSPPAGSPAGRSIYNSFYVYCKGPCQRVQPGKLRVQCSTCRQATLTLTQGPSCWDDVLIPNRMSGECQSPHCPGTSAEFFFKCGAHPTSDKETSVALHLIATNSRNITCITCTDVRSPVLVFQCNSRHVICLDCFHLYCVTRLNDRQFVHDPQLGYSLPCVAGCPNSLIKELHHFRILGEEQYNRYQQYGAEECVLQMGGVLCPRPGCGAGLLPEPDQRKVTCEGGNGLGCGFAFCRECKEAYHEGECSAVFEASGTTTQAYRVDERAAEQARWEAASKETIKKTTKPCPRCHVPVEKNGGCMHMKCPQPQCRLEWCWNCGCEWNRVCMGDHWFDV"


##
##  Build a data frame of coutns per amino acid variant
##
# consider keeping syn. WT here but still merging different DNA variants of the same aa_variants
counts = read.csv(counts_file)
stopifnot(all( c("subst_aa","subst_dna",cn_counts) %in% colnames(counts) ))

# First, find all single amino acid variants
ic_aa1 = which(counts$n_subst_aa==1)
print(sprintf("Aggregating read counts for %d of %d (%.1f%%) barcodes of single amino acid variants",
              length(ic_aa1), nrow(counts), length(ic_aa1)/nrow(counts)*100))
for (cn in cn_counts) {
    print(sprintf("   using %8d of %8d (%.1f%%) read counts for %s", sum(counts[ic_aa1,cn]), sum(counts[,cn]), sum(counts[ic_aa1,cn])/sum(counts[,cn])*100, cn))
}

# aggregate barcode counts per amino acid substitution
df = counts[ic_aa1,c("subst_aa",cn_counts)]
agg = aggregate(df[,cn_counts], by=list(df$subst_aa), sum)
aa1 = data.frame(subst=agg[,1], wt=NA, resi=NA, mut=NA)
aa1[,cn_counts] = agg[,cn_counts]

# extract amino acids and residue number from substitution string
nc = nchar(aa1$subst)
aa1$wt = substr(aa1$subst,1,1)
aa1$resi = as.numeric( substr(aa1$subst,2,nc-1) )
aa1$mut = substr(aa1$subst,nc,nc)

# re-oreder variants
aa1 = aa1[order(aa1$resi,aa1$mut),]
rownames(aa1) = NULL

# Second, synonymous wild-type variants with a single DNA substitutions
ic_syn1 = which(counts$n_subst_aa == 0 & counts$n_subst_dna == 1)
print(sprintf("Aggregating read counts for %d of %d (%.1f%%) barcodes of single nucl synonymous wild types",
              length(ic_syn1), nrow(counts), length(ic_syn1)/nrow(counts)*100))
for (cn in cn_counts) {
    print(sprintf("   using %8d of %8d (%.1f%%) read counts for %s", sum(counts[ic_syn1,cn]), sum(counts[,cn]), sum(counts[ic_syn1,cn])/sum(counts[,cn])*100, cn))
}

# aggregate barcode counts per nucleotide substitution
df = counts[ic_syn1,c("subst_dna",cn_counts)]

# calculate amino acid numbering and make substitution string
df$resi = (as.numeric( substr(df$subst_dna,2,nchar(df$subst_dna)-1) )-1) %/% 3 + 1
df$wt = strsplit(wt$aa,"")[[1]][df$resi]
df$subst_aa = paste0(df$wt, df$resi, "=")

agg = aggregate(df[,cn_counts], by=list(df$subst_aa), sum)
syn1 = data.frame(subst=agg[,1], wt=NA, resi=NA, mut=NA)
syn1[,cn_counts] = agg[,cn_counts]

# extract amino acids and residue number from substitution string
nc = nchar(syn1$subst)
syn1$wt = substr(syn1$subst,1,1)
syn1$resi = as.numeric( substr(syn1$subst,2,nc-1) )
syn1$mut = substr(syn1$subst,nc,nc)

# re-oreder variants
syn1 = syn1[order(syn1$resi),]
rownames(syn1) = NULL

# Third, merge barcodes for wildtype
ic_wt = which(counts$n_subst_aa == 0 & counts$n_subst_dna == 0)
print(sprintf("Aggregating read counts for %d of %d (%.1f%%) barcodes of wild type (DNA level)",
              length(ic_wt), nrow(counts), length(ic_wt)/nrow(counts)*100))
for (cn in cn_counts) {
    print(sprintf("   using %8d of %8d (%.1f%%) read counts for %s", sum(counts[ic_wt,cn]), sum(counts[,cn]), sum(counts[ic_wt,cn])/sum(counts[,cn])*100, cn))
}
wtrow = data.frame(subst="WT", wt="", resi=0, mut="")
wtrow[,cn_counts] = apply(counts[ic_wt,cn_counts], MARGIN=2, sum)

# Merge into a single data frame
raw = rbind(wtrow, aa1, syn1)

# All substitutions should be unique
stopifnot( nrow(raw) == length(unique(raw$subst)) )

# Each barcode should only be used once
ic = c(ic_wt, ic_aa1, ic_syn1)
stopifnot( length(ic) == length(unique(ic)) )

print(sprintf("Collected data for %d substitutions using %d barcodes (%.1f%%)", nrow(raw), length(ic), length(ic)/nrow(counts)*100))


##
##  Calculate abundance scores
##

# Calculate RPM with psudo counts 
settings$pseudocounts = 1
calc_rpm = function(v) { v/sum(v)*10^6 }
raw_rpm = data.frame(subst=raw$subst)
for (cn in cn_counts) {
    raw_rpm[,cn] = calc_rpm(raw[,cn] + settings$pseudocounts)
}

# Function to calculate PSI based on rows in a data frame
protein_stability_index = function(df, name, indices=NA, populations=NA) {
    # Each FACS gate should have an index
    if (all(is.na(indices))) {
        indices = seq(ncol(df))
    } else {
        stopifnot( length(indices) == ncol(df) )
    }

    calc_psi = function(rpm,idx) {
        rpm_sum = sum(rpm)
	psi = sum( rpm/rpm_sum * idx )
	return( c(rpm_sum, psi) )
    }
    calc_psi_pop = function(rpm,idx,pop) {
        rpm_sum = sum(rpm * pop)
	psi = sum( rpm/rpm_sum * idx * pop )
	return( c(rpm_sum, psi) )
    }

    # Each FACS gate should have a population
    if (all(is.na(populations))) {
        # populations = rep(1/ncol(df), times=ncol(df))
	print(sprintf("    Calc. PSI for %s using indices %s",name,paste0(indices,collapse=",")))
	psi = t(apply(df, MARGIN=1, calc_psi, idx=indices))
    } else {
        stopifnot( length(populations) == ncol(df) )
	print(sprintf("    Calc. PSI for %s using indices %s and populations %s",name,paste0(indices,collapse=","),paste0(populations,collapse=",")))
	psi = t(apply(df, MARGIN=1, calc_psi_pop, idx=indices, pop=populations))
    }
    
    # put in a data frame and give column name
    ret_df = data.frame(psi)
    colnames(ret_df) = paste0(name,c("_rpm_sum","_psi"))
    return(ret_df)
}


# Required number of reads in all FACS bins needed to trust PSI of a variant
settings$threshold_counts_per_rep = 50

# Calculate PSI
raw_psi = data.frame(subst=raw$subst)
for (cell in names(settings$bio_reps)) {
    for (bio_rep in settings$bio_reps[[cell]]) {
        replica = sprintf("%s_%s",cell,bio_rep)
        cn_gates = sprintf("%s_bin%d.%s",cell,seq(settings$n_facs_bins),bio_rep)
        cn_reads = sprintf("%s_read_sum",replica)
        print(sprintf("Calculate PSI for %s",replica))
    
        # Total number of reads per variant over all gates
        raw_psi[,cn_reads] = apply(raw[,cn_gates], MARGIN=1, sum)

        # PSI per variant
        raw_psi = cbind(raw_psi, protein_stability_index(raw_rpm[,cn_gates],  name = replica,
	                                                 indices = settings$idx_facs_bins,
							 populations = settings$pop_facs_bins[[cell]]))

        i_rm = which(raw_psi[,cn_reads] < settings$threshold_counts_per_rep)
	if (length(i_rm) > 0) {
            print(sprintf("    Removing %d PSI scores with less than %d total reads",length(i_rm),settings$threshold_counts_per_rep))
            raw_psi[i_rm,sprintf("%s_psi",replica)] = NA
	}
    }
}

# plot(raw_psi$bio1_facs1_psi, raw_psi$bio1_facs2_psi)
# i = which(raw_psi$bio1_facs1_read_sum < 100 | raw_psi$bio1_facs2_read_sum < 100)
# points(raw_psi[i,"bio1_facs1_psi"], raw_psi[i,"bio1_facs2_psi"], pch=20, col=2)


# Calculate a mean of replica
for (cell in names(settings$bio_reps)) {
    replica_cns = sprintf("%s_%s_psi",cell,settings$bio_reps[[cell]])
    raw_psi[,paste0(cell,"_mean")] = apply(raw_psi[,replica_cns], MARGIN=1, mean, na.rm=T)
    raw_psi[,paste0(cell,"_sd")]   = apply(raw_psi[,replica_cns], MARGIN=1, sd, na.rm=T)
}


##
##  Per residue data
##

# coverage, which single amino acid subst have calculated scores
aa_one = strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
stopifnot(all(raw$subst == raw_psi$subst))

# only count mutations that have both conditions measured
i_aa1 = which(raw$mut %in% aa_one & ! is.na(raw_psi$delta_psi))
agg = aggregate(raw[i_aa1,"mut"], by=list(raw[i_aa1,"resi"]), paste0, collapse="")

residue = data.frame(resi=seq(nchar(wt[["aa"]])), wt=strsplit(wt[["aa"]],"")[[1]])
residue[,"mut"] = agg[match(residue$resi,agg[,1]),2]
residue[which(is.na(residue$mut)),"mut"] = ""
residue$coverage = sapply(residue$mut, nchar)

agg = aggregate(raw_psi[i_aa1,"ctrl_mean"], by=list(raw[i_aa1,"resi"]), median, na.rm=T)
residue$med_ctrl = agg[match(residue$resi,agg[,1]),2]

agg = aggregate(raw_psi[i_aa1,"b6ko_mean"], by=list(raw[i_aa1,"resi"]), median, na.rm=T)
residue$med_bag6 = agg[match(residue$resi,agg[,1]),2]

in_ctrl = which(raw$mut == "*" & ! is.na(raw_psi$ctrl_mean))
in_ctrl_match = in_ctrl[match(residue$resi,raw[in_ctrl,"resi"])]
residue$nons_ctrl = raw_psi[in_ctrl_match,"ctrl_mean"]

in_bag6 = which(raw$mut == "*" & ! is.na(raw_psi$b6ko_mean))
in_bag6_match = in_bag6[match(residue$resi,raw[in_bag6,"resi"])]
residue$nons_bag6 = raw_psi[in_bag6_match,"b6ko_mean"]

is_ctrl = which(raw$mut == "=" & ! is.na(raw_psi$ctrl_mean))
is_ctrl_match = is_ctrl[match(residue$resi,raw[is_ctrl,"resi"])]
residue$syno_ctrl = raw_psi[is_ctrl_match,"ctrl_mean"]

is_bag6 = which(raw$mut == "=" & ! is.na(raw_psi$b6ko_mean))
is_bag6_match = is_bag6[match(residue$resi,raw[is_bag6,"resi"])]
residue$syno_bag6 = raw_psi[is_bag6_match,"b6ko_mean"]

# report coverage
print(sprintf("Data has %d delta_PSI measurements,  %d of %d (%.1f%%) single variants",
              sum(! is.na(raw_psi$delta_psi)),
              length(i_aa1), 19*nchar(wt[["aa"]]), 100*length(i_aa1)/(19*nchar(wt[["aa"]])) ))
print(sprintf("  For ctrl cells,  %d (%.1f%%) nonsense, %d (%.1f%%) synonymous and WT",
	      length(in_ctrl), 100*length(in_ctrl)/nchar(wt[["aa"]]),
	      length(is_ctrl), 100*length(is_ctrl)/nchar(wt[["aa"]]) ))
print(sprintf("  For bag6 cells,  %d (%.1f%%) nonsense, %d (%.1f%%) synonymous and WT",
	      length(in_bag6), 100*length(in_bag6)/nchar(wt[["aa"]]),
	      length(is_bag6), 100*length(is_bag6)/nchar(wt[["aa"]]) ))


#
# Dump 
#

# new data frame without missing scores
i = which(! (is.na(raw_psi$ctrl_mean) & is.na(raw_psi$b6ko_mean)))
prkn_b6 = raw_psi[i,c("subst","ctrl_mean","ctrl_sd","b6ko_mean","b6ko_sd","delta_psi")]
colnames(prkn_b6) = c("var","ctrl_mean","ctrl_sd","b6ko_mean","b6ko_sd","delta_psi")
not_both_measured = (is.na(prkn_b6$ctrl_mean) + is.na(prkn_b6$b6ko_mean)) != 0
stopifnot(all( not_both_measured == is.na(prkn_b6$delta_psi) ))
print(sprintf("Final score set has %d variants, incl. %d with one celltype missing, %d nonsense and %d synonymous)",
              nrow(prkn_b6), sum(not_both_measured), sum(raw[i,"mut"]=='*'), sum(raw[i,"mut"]=='=')))

# Dump CSV with substitutions
write.csv(prkn_b6, file="prkn_bag6.csv", row.names=F, quote=F)

# Dump CSV with per-residues data
write.csv(residue, file="prkn_bag6_residues.csv", row.names=F, quote=F)

# save everything in R data file
save(raw, raw_psi, residue, prkn_b6, wt, settings, cn_counts, replicas, file="prkn_bag6.rda")

print("")
print(settings)


#
# Plots and analyses
#

# plot distributions
hc = hist(prkn_b6$ctrl_mean, breaks=30, plot=F)
hb = hist(prkn_b6$b6ko_mean, breaks=30, plot=F)
quartz(width=8, height=5)
plot(0,0,col=0, xlim=c(1,2), ylim=c(0,1), xlab="PSI", ylab="Density")
lines(hc$mids, hc$density, lwd=2, col=1)
lines(hb$mids, hb$density, lwd=2, col=2)
legend("top", c("Ctrl strain","BAG6 KO"), lty=1, lwd=2, col=c(1,2))
quartz.save("prkn_bag6_distributions.png", type="png")

# plot distributions
i_nons = which(grepl("*",prkn_b6$var,fixed=T))
i_syno = which(grepl("=",prkn_b6$var,fixed=T))
breaks = seq(1,2,0.02)
hca = hist(prkn_b6$ctrl_mean, breaks=breaks, plot=F)
hcn = hist(prkn_b6[i_nons,"ctrl_mean"], breaks=breaks, plot=F)
hcs = hist(prkn_b6[i_syno,"ctrl_mean"], breaks=breaks, plot=F)
hba = hist(prkn_b6$b6ko_mean, breaks=breaks, plot=F)
hbn = hist(prkn_b6[i_nons,"b6ko_mean"], breaks=breaks, plot=F)
hbs = hist(prkn_b6[i_syno,"b6ko_mean"], breaks=breaks, plot=F)

quartz(width=12, height=5)
par(mfrow=c(1,2))
plot(0,0,col=0, xlim=c(1,2), ylim=c(0,3), xlab="PSI", ylab="Density", main="Ctrl strain")
lines(hca$mids, hca$density, lwd=2, col=1)
lines(hcn$mids, hcn$density, lwd=2, col=2)
lines(hcs$mids, hcs$density, lwd=2, col=3)
legend("top", c("All","Nonsense","Synonymous"), lty=1, lwd=2, col=c(1,2,3))
plot(0,0,col=0, xlim=c(1,2), ylim=c(0,3), xlab="PSI", ylab="Density", main="BAG6 KO strain")
lines(hba$mids, hba$density, lwd=2, col=1)
lines(hbn$mids, hbn$density, lwd=2, col=2)
lines(hbs$mids, hbs$density, lwd=2, col=3)
legend("top", c("All","Nonsense","Synonymous"), lty=1, lwd=2, col=c(1,2,3))
quartz.save("prkn_bag6_nons_syn.png", type="png")


# plot replica PSI correlations
quartz(width=10, height=5)
par(mfrow=c(1,2))
rpc = cor(raw_psi$ctrl_bio1_psi, raw_psi$ctrl_bio2_psi, use="complete.obs", method="pearson")
plot(raw_psi$ctrl_bio1_psi, raw_psi$ctrl_bio2_psi, pch=20, cex=.2, main=sprintf("Ctrl. cells, Pearson %.2f",rpc),
     xlab="PSI, biol. rep. 1", ylab="PSI, biol. rep. 2")
rpb = cor(raw_psi$b6ko_bio1_psi, raw_psi$b6ko_bio2_psi, use="complete.obs", method="pearson")
plot(raw_psi$b6ko_bio1_psi, raw_psi$b6ko_bio2_psi, pch=20, cex=.2, main=sprintf("BAG6 KO cells, Pearson %.2f",rpb),
     xlab="PSI, biol. rep. 1", ylab="PSI, biol. rep. 2")
quartz.save("prkn_rep_cor.png", type="png")

