#!/bin/bash
# Generate basic fastq statistics
# David Coffey
# January 5, 2018

# Specify variables
#export PREFIX1="Sample_1.R1"
#export PREFIX2="Sample_1.R1"
#export READ1="Sample_1.R1.fastq.gz"
#export READ2="Sample_1.R1.fastq.gz"
#export FASTQSTATS_DIRECTORY="$ROOT/Fastq/Stats"
#export FASTQSTATS="/home/dcoffey/Apps/ea-utils/clipper/fastq-stats"

# Start time
START=`date +%s`
echo Begin FastqStats.sh for file $FILE on `date +"%B %d, %Y at %r"`

# Make directory for deduplicated files
mkdir -p $FASTQSTATS_DIRECTORY

# Create fastq stats for read 1
$FASTQSTATS $READ1 \
-x $FASTQSTATS_DIRECTORY/$PREFIX1.fastx.stats.txt \
-b $FASTQSTATS_DIRECTORY/$PREFIX1.base.quality.txt > $FASTQSTATS_DIRECTORY/$PREFIX1.fastq.stats.txt

# Create fastq stats for read 2
$FASTQSTATS $READ2 \
-x $FASTQSTATS_DIRECTORY/$PREFIX2.fastx.stats.txt \
-b $FASTQSTATS_DIRECTORY/$PREFIX2.base.quality.txt > $FASTQSTATS_DIRECTORY/$PREFIX2.fastq.stats.txt

END=`date +%s`
MINUTES=$(((END-START)/60))
