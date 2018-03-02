#!/bin/bash
# Deduplicate alignment and count genes
# Written by David Coffey dcoffey@fhcrc.org
# Updated December 27, 2017

# Specify variables
#export BAM="Alignment.bam"
#export DEDUPLICATED_DIRECTORY="$ROOT/Aligned/deduplicated"

# Start time
START=`date +%s`
echo Begin Deduplicate.sh for sample $SAMPLE cell $CELL on `date +"%B %d, %Y at %r"`

# Make directory for deduplicated files
mkdir -p $DEDUPLICATED_DIRECTORY

# Deduplicate BAM file
module load Python/3.6.3-foss-2016b-fh2
umi_tools dedup \
-I $BAM \
--output-stats=$DEDUPLICATED_DIRECTORY/$SAMPLE.$CELL.deduplicated.stats \
-S $DEDUPLICATED_DIRECTORY/$SAMPLE.$CELL.deduplicated.bam

# Extract feature counts 
$SUBREAD/featureCounts -a $HG19_GTF -o $DEDUPLICATED_DIRECTORY/$SAMPLE.$CELL.gene.assigned -R BAM $DEDUPLICATED_DIRECTORY/$SAMPLE.$CELL.deduplicated.bam -T 4

# Sort and index feature counts BAM
module purge
module load SAMtools/1.6-foss-2016b
samtools sort $DEDUPLICATED_DIRECTORY/$SAMPLE.$CELL.deduplicated.bam.featureCounts.bam -o $DEDUPLICATED_DIRECTORY/$SAMPLE.$CELL.deduplicated.featureCounts.bam
rm $DEDUPLICATED_DIRECTORY/$SAMPLE.$CELL.deduplicated.bam.featureCounts.bam

samtools index $DEDUPLICATED_DIRECTORY/$SAMPLE.$CELL.deduplicated.featureCounts.bam
mv $DEDUPLICATED_DIRECTORY/$SAMPLE.$CELL.deduplicated.featureCounts.bam.bai $DEDUPLICATED_DIRECTORY/$SAMPLE.$CELL.deduplicated.featureCounts.bai

# Count genes
module load Python/3.6.2-foss-2016b-fh2
umi_tools count \
--per-gene \
--gene-tag=XT \
-I $DEDUPLICATED_DIRECTORY/$SAMPLE.$CELL.deduplicated.featureCounts.bam \
-S $DEDUPLICATED_DIRECTORY/$SAMPLE.$CELL.deduplicated.counts.tsv

END=`date +%s`
MINUTES=$(((END-START)/60))
