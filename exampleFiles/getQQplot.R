# Example R script to produce QQ plot for SLPs
# Based on script written by Chris Rayner

# Assumes there is a summary file containing all SLPs called e.g. ajsz.ct08.rare.summ.txt
# Substitute your own file name below.
# If no recessive analysis has been performed (so each .sao file contains only one SLP)
# then a summary file can be produced with a command like this:

# grep "SLP =" *.sao > ajsz.ct08.rare.summ.txt

# The first three lines of the summary file might then look like this:
# ajsz.ct08.rare.A1BG-AS1.sao:SLP =     0.00 (signed log10(p), positive if cases score higher than controls)
# ajsz.ct08.rare.A1BG.sao:SLP =    -0.68 (signed log10(p), positive if cases score higher than controls)
# ajsz.ct08.rare.A1CF.sao:SLP =     0.08 (signed log10(p), positive if cases score higher than controls)

# The SLP is then in the third column and the following commands should work to produce a plot of SLP against expected SLP

tab=read.table("ajsz.ct08.rare.summ.txt",header=FALSE)
names(tab) <- c("saoFile","eq","SLP")
attach (tab)
rankSLP <- rank(SLP, na.last=TRUE, ties.method="first")
nRows=nrow(tab)
midrank <- (nRows+1)/2
rMinusMr <- (rankSLP-midrank)
absDiff <- abs(rMinusMr/midrank)
pVal <- 1-absDiff
logP <- log10(pVal)
eSLP <- sign(rankSLP - midrank) * -logP
plot(SLP,eSLP)

