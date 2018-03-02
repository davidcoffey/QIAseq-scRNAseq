#!/bin/bash
# De-multiplex fastq files
# David Coffey
# January 14, 2018

# Specify variables
#export DEMULTIPLEX_DIRECTORY="$ROOT/fastq_demultiplexed"
#export SAMPLE="Sample.ID"

# Start time
START=`date +%s`
echo Begin Demultiplex.sh for sample $SAMPLE on `date +"%B %d, %Y at %r"`

# Make directory for the demultiplexed fastq files
mkdir -p $DEMULTIPLEX_DIRECTORY

# Trim reverse compliment 3' adapter, TS oligo, and UMI from read 1
module load Python/3.6.3-foss-2016b-fh2
cutadapt -m 12 -O 20 -e 0 -a file:$ROOT/Adapters/reverse.adapters.fasta -o $DEMULTIPLEX_DIRECTORY/$SAMPLE.trimmed.R1.fastq.gz -p $DEMULTIPLEX_DIRECTORY/$SAMPLE.trimmed.R2.fastq.gz $READ1 $READ2

# Demultiplex
$FASTQMULTX -m 2 -b \
-B $ROOT/Adapters/adpaters.txt \
$DEMULTIPLEX_DIRECTORY/$SAMPLE.trimmed.R2.fastq.gz \
$DEMULTIPLEX_DIRECTORY/$SAMPLE.trimmed.R1.fastq.gz \
-o $DEMULTIPLEX_DIRECTORY/$SAMPLE.%.R2.fastq.gz \
-o $DEMULTIPLEX_DIRECTORY/$SAMPLE.%.R1.fastq.gz

END=`date +%s`
MINUTES=$(((END-START)/60))
