#!/bin/bash
# UMI extraction
# Written by David Coffey dcoffey@fhcrc.org
# Updated January 14, 2017

# Specify variables
#export UMI_EXTRACTED_DIRECTORY="$ROOT/fastq_demultiplexed"
#export SAMPLE="Sample.ID"
#export CELL="Cell.ID"


# Start time
START=`date +%s`
echo Begin ExtractUMI.sh for sample $SAMPLE cell $CELL on `date +"%B %d, %Y at %r"`

# Make directory for the demultiplexed fastq files
mkdir -p $UMI_EXTRACTED_DIRECTORY

# Load python which has installed umi_tools version 0.5.2 and cutadapt version 1.15
module load Python/3.6.3-foss-2016b-fh2

# Remove reads less than 12 nt in read 2
cutadapt -m 12 -o $DEMULTIPLEX_DIRECTORY/$SAMPLE.$CELL.R2.filtered.fastq.gz -p $DEMULTIPLEX_DIRECTORY/$SAMPLE.$CELL.R1.filtered.fastq.gz $DEMULTIPLEX_DIRECTORY/$SAMPLE.$CELL.R2.fastq.gz $DEMULTIPLEX_DIRECTORY/$SAMPLE.$CELL.R1.fastq.gz

# Extract UMI from read 2
umi_tools extract \
--stdin=$DEMULTIPLEX_DIRECTORY/$SAMPLE.$CELL.R1.filtered.fastq.gz \
--read2-in=$DEMULTIPLEX_DIRECTORY/$SAMPLE.$CELL.R2.filtered.fastq.gz \
--stdout=$UMI_EXTRACTED_DIRECTORY/$SAMPLE.$CELL.R1.UMI.fastq.gz \
--read2-out=$UMI_EXTRACTED_DIRECTORY/$SAMPLE.$CELL.R2.UMI.fastq.gz \
--extract-method=string \
--bc-pattern2=NNNNNNNNNNNN

END=`date +%s`
MINUTES=$(((END-START)/60))
