#!/bin/bash
# Align RNAseq to reference genome
# Written by David Coffey dcoffey@fhcrc.org
# Updated December 26, 2017


## Prerequisites (see Software_installation.sh)
# Download and install STAR aligner
# Create STAR genome using STAR_genome.sh script

## Variables
# export GENOME="Hg19"
# export GENOME_DIRECTORY=".../STAR_genomes/$GENOME"
# export SAMPLE="..."
# export READ1=".../R1.fastq.gz" # Fastq file must be gzipped
# export READ2=".../R2.fastq.gz" # Fastq file must be gzipped
# export PREFIX=".../STAR_alignment/$GENOME/$SAMPLE/$SAMPLE"
# export ALIGNMENT_DIRECTORY=".../STAR_alignment/$GENOME/$SAMPLE"
# export SAMTOOLS=".../samtools"

START=`date +%s`
echo Begin STAR_alignment.sh for sample $SAMPLE on `date +"%B %d, %Y at %r"`

# Align FASTQ file to reference genome
mkdir -p $ALIGNMENT_DIRECTORY

module load STAR/2.5.0a-foss-2015b

STAR \
--genomeDir $GENOME_DIRECTORY \
--readFilesIn $READ1 \
--readFilesCommand zcat \
--runThreadN 4 \
--quantMode GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix $PREFIX.

# Index BAM file
module load SAMtools/1.6-foss-2016b

samtools index $PREFIX.Aligned.sortedByCoord.out.bam
mv -f $PREFIX.Aligned.sortedByCoord.out.bam.bai $PREFIX.Aligned.sortedByCoord.out.bai

# Clean up
rm -R $PREFIX._STARtmp

END=`date +%s`
MINUTES=$(((END-START)/60))
echo End STAR_alignment.sh for sample $SAMPLE.  The run time was $MINUTES minutes.