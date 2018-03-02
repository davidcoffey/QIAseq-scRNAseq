#!/bin/bash
# Peform QC analysis on fastq files
# Adapted from https://www.broadinstitute.org/gatk/guide/article?id=3891
# Written by David Coffey dcoffey@fhcrc.org
# Updated February 17, 2016

## Prerequisites (see Software_installation.sh)
# Download and install FastQC

## Variables
# export FASTQC=".../fastqc" if using the app
# export FASTQC_DIRECTORY="..."
# export FASTQ_FILES="..."
# export SAMPLE="..."

# Load fastqc on rhino
module load FastQC/0.11.5-Java-1.8.0_92

START=`date +%s`
echo Begin FastQC/0.11.5-Java-1.8.0_92 for sample $SAMPLE on `date +"%B %d, %Y at %r"`

mkdir -p $FASTQC_DIRECTORY/$SAMPLE

# Run fastqc
zcat $FASTQ_FILES/*$SAMPLE* | fastqc stdin --outdir=$FASTQC_DIRECTORY/$SAMPLE

# Change filename in fastqc report
find $FASTQC_DIRECTORY/$SAMPLE/stdin_fastqc.html -type f -print0 | xargs -0 sed -i "s/stdin/$SAMPLE/g"

# Rename html file
mv $FASTQC_DIRECTORY/$SAMPLE/stdin_fastqc.html $FASTQC_DIRECTORY/$SAMPLE.fastqc.html

# Clean up
rm -R $FASTQC_DIRECTORY/$SAMPLE

END=`date +%s`
MINUTES=$(((END-START)/60))
echo End GATK.sh for sample $SAMPLE.  The run time was $MINUTES minutes.

