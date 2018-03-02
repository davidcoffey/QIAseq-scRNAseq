#!/bin/bash
# Generate basic bam statistics
# David Coffey
# January 5, 2018

# Specify variables
#export PREFIX="Sample_1"
#export BAM="file.BAM.gz"
#export BAMSTATS_DIRECTORY="$ROOT/BAM/Stats"
#export BAMSTATS="/home/dcoffey/Apps/ea-utils/clipper/sam-stats"

# Start time
START=`date +%s`
echo Begin BamStats.sh for file $FILE on `date +"%B %d, %Y at %r"`

# Make directory for deduplicated files
mkdir -p $BAMSTATS_DIRECTORY

# Create BAM stats
$BAMSTATS $BAM \
-DBz \
-R $BAMSTATS_DIRECTORY/$PREFIX.coverage.txt > $BAMSTATS_DIRECTORY/$PREFIX.bam.stats.txt

END=`date +%s`
MINUTES=$(((END-START)/60))