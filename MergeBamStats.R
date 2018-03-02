#!/usr/bin/env Rscript
# Merge bam stats files
# David Coffey dcoffey@fredhutch.org
# January 5, 2017

library(data.table)
library(plyr)

# Import environmental variable and make destination directory
bamstats.directory = Sys.getenv("BAMSTATS_DIRECTORY")
system(command = "mkdir $BAMSTATS_DIRECTORY/Summary")

# Define bam paths
paths.bam = list.files(bamstats.directory, full.names = TRUE, all.files = FALSE, recursive = FALSE, include.dirs = FALSE)
paths.bam = paths.bam[grep(paths.bam, pattern = "bam.stats.txt")]
file.info.bam = file.info(paths.bam)
paths.bam = rownames(file.info.bam)[file.info.bam$size > 0]

# Merge bam stats
bstats.list = list()

p = 1
for(p in 1:length(paths.bam)){
  bstats = fread(file = paths.bam[p], fill = TRUE)
  names(bstats) = c("measure", gsub(x = basename(paths.bam[p]), pattern = ".bam.stats.txt", replacement = ""))
  bstats.list = c(bstats.list, list(bstats))
}

bstats.df = Reduce(function(...) merge(..., all = TRUE), bstats.list)
bstats.df = stats::setNames(data.frame(t(bstats.df[,-1])), bstats.df$measure)
bstats.df$samples = rownames(bstats.df)
rownames(bstats.df) = NULL
bstats.df = bstats.df[,c("samples", setdiff(names(bstats.df), "samples"))]
write.csv(bstats.df, file = paste(bamstats.directory, "Summary/Bam stats.csv", sep = "/"), row.names = FALSE)

# Define coverage paths
paths.coverage = list.files(bamstats.directory, full.names = TRUE, all.files = FALSE, recursive = FALSE, include.dirs = FALSE)
paths.coverage = paths.coverage[grep(paths.coverage, pattern = "coverage.txt")]
file.info.coverage = file.info(paths.coverage)
paths.coverage = rownames(file.info.coverage)[file.info.coverage$size > 0]

# Merge coverage stats
coverage.list = list()

p = 1
for(p in 1:length(paths.coverage)){
  coverage = fread(file = paths.coverage[p], fill = TRUE)
  coverage = coverage[order(coverage$V1, decreasing = TRUE),]
  names(coverage) = c("transcript", "length", "count"," percent coverage", "skewness", "coverage CV", "signature")
  coverage.list = c(coverage.list, list(coverage))
  names(coverage.list)[p] = gsub(x = basename(paths.coverage[p]), pattern = ".coverage.txt", replacement = "")
}

coverage.df = ldply(coverage.list, data.frame)
names(coverage.df)[1] = "samples"
write.csv(coverage.df, file = paste(bamstats.directory, "Summary/Coverage dedup stats.csv", sep = "/"), row.names = FALSE)