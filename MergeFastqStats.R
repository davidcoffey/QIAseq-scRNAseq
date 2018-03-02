#!/usr/bin/env Rscript
# Merge fastq stats files
# David Coffey dcoffey@fredhutch.org
# January 5, 2017

library(data.table)
library(plyr)

# Import environmental variables
fastqstats.directory = Sys.getenv("FASTQSTATS_DIRECTORY")
system(command = "mkdir $FASTQSTATS_DIRECTORY/Summary")

# Define fastq paths
paths.fastq = list.files(fastqstats.directory, full.names = TRUE, all.files = FALSE, recursive = FALSE, include.dirs = FALSE)
paths.fastq = paths.fastq[grep(paths.fastq, pattern = ".fastq.stats.txt")]
file.info.fastq = file.info(paths.fastq)
paths.fastq = rownames(file.info.fastq)[file.info.fastq$size > 0]

# Merge fastq stats
fstats.list = list()
dup.list = list()

p = 1
for(p in 1:length(paths.fastq)){
  fstats = fread(file = paths.fastq[p], fill = TRUE)
  dup = fstats[13:22,c(4,3)]
  fstats = fstats[c(1:12,23:34),1:2]
  names(dup) = c("sequence", "count")
  names(fstats) = c("measure", gsub(x = basename(paths.fastq[p]), pattern = ".fastq.stats.txt", replacement = ""))
  fstats.list = c(fstats.list, list(fstats))
  dup.list = c(dup.list, list(dup))
  names(dup.list)[p] = gsub(x = basename(paths.fastq[p]), pattern = ".fastq.stats.txt", replacement = "")
}

fstats.df =  Reduce(function(...) merge(..., all = TRUE), fstats.list)
fstats.df = stats::setNames(data.frame(t(fstats.df[,-1])), fstats.df$measure)
fstats.df$samples = rownames(fstats.df)
rownames(fstats.df) = NULL
fstats.df = fstats.df[,c("samples", setdiff(names(fstats.df), "samples"))]
write.csv(fstats.df, file = paste(fastqstats.directory, "Summary/Fastq stats.csv", sep = "/"), row.names = FALSE)

dup.df = ldply(dup.list, data.frame)
names(dup.df)[1] = "samples"
write.csv(dup.df, file = paste(fastqstats.directory, "Summary/Duplication stats.csv", sep = "/"), row.names = FALSE)

# Define fastx paths
paths.fastx = list.files(fastqstats.directory, full.names = TRUE, all.files = FALSE, recursive = FALSE, include.dirs = FALSE)
paths.fastx = paths.fastx[grep(paths.fastx, pattern = ".fastx.stats.txt")]
file.info.fastx = file.info(paths.fastx)
paths.fastx = rownames(file.info.fastx)[file.info.fastx$size > 0]

# Merge fastx files
fastx.list = list()

p = 1
for(p in 1:length(paths.fastq)){
  fastx = fread(file = paths.fastx[p], fill = TRUE)
  fastx.list = c(fastx.list, list(fastx))
  names(fastx.list)[p] = gsub(x = basename(paths.fastx[p]), pattern = ".fastx.stats.txt", replacement = "")
}

fastx.df = ldply(fastx.list, data.frame)
names(fastx.df)[1] = "samples"
write.csv(fastx.df, file = paste(fastqstats.directory, "Summary/Fastx stats.csv", sep = "/"), row.names = FALSE)