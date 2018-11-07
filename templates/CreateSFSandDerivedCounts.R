#!/usr/bin/env Rscript
# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2017/08/01

## Generate SFS

pdf("${sfs_plot}")
frqFile <- read.delim2("${frq_input}", skip=1, stringsAsFactors = FALSE, header=FALSE)
derFrqCol <- as.vector(frqFile[,6])
derFrq <- as.numeric(sapply(strsplit(derFrqCol, ":"), "[[", 2))
hist(round(as.numeric(derFrq), digits=2), xlab="Frequency", ylab="Number of variants", main=paste("SFS", "${popID}", sep=" "), col="red")
dev.off()

