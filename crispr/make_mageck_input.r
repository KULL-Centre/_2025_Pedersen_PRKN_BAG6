
load("counts.rda")
df = data.frame(sgRNA = counts$GUIDE_ID, gene = counts$gene,
	       T0 = counts[,"T0_S1"],
	       T12 = counts[,"T12-No-Facs_S2"],
	       T12.low = counts[,"5-Lower_S4"],
	       T12.high = counts[,"5-Upper_S3"])
write.table(df, "counts.txt", quote=F, row.names=F, sep="\t")

