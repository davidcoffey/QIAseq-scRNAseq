#!/bin/bash
# scRNAseq pipeline
# Written by David Coffey dcoffey@fhcrc.org
# Updated December 26, 2017

# Variables
export ROOT=".../Root/Folder"
export MIXCR_DIRECTORY="$ROOT/MiXCR"
export GENOME_DIRECTORY=".../hSTAR_genomes/hg19"
export HG19_GTF=".../hg19.gtf"
export DEMULTIPLEX_DIRECTORY="$ROOT/Fastq/Demultiplexed"
export ALIGNMENT_DIRECTORY="$ROOT/Aligned/Original/"
export DEDUPLICATED_DIRECTORY="$ROOT/Aligned/Deduplicated/"
export UMI_EXTRACTED_DIRECTORY="$ROOT/Fastq/UMI_extracted"
export FASTQSTATS_DIRECTORY="$ROOT/Fastq/Stats"
export FASTQMULTX=".../ea-utils/clipper/fastq-multx"
export FASTQSTATS="..../ea-utils/clipper/fastq-stats"
export BAMSTATS=".../ea-utils/clipper/sam-stats"
export SUBREAD=".../subread-1.6.0-Linux-x86_64/bin"

export SAMPLES="
Plate1
Plate2..."

export CELLS="
Cell1
Cell2..."

# Trim and demultiplex
for S in ${SAMPLES}; do
    echo ${S}
    export SAMPLE=${S}
    export READ1="$ROOT/Fastq/Original/${S}_R1_001.fastq.gz"
    export READ2="$ROOT/Fastq/Original/${S}_R2_001.fastq.gz"
    sbatch -n 1 -c 4 -t 3-0 --job-name="DEMULTI" --output=$ROOT/Logs/Demultiplex.${S}.log $ROOT/Scripts/Demultiplex.sh
done

DEMULTI=$(squeue -o "%A" -h -u dcoffey -n "DEMULTI" -S i | tr "\n" ":")

# Extract antigen receptor sequences with MiXCR
for S in ${SAMPLES}; do
    echo ${S}
    export SAMPLE=${S}
    for C in ${CELLS}; do
        echo ${C}
        export CELL=${C}
        export READ1="$DEMULTIPLEX_DIRECTORY/$SAMPLE.$CELL.R1.filtered.fastq.gz"
        sbatch -n 1 -c 4 -t 3-0 --job-name="MiXCR" --dependency=afterany:${DEMULTI%?} --output=$ROOT/Logs/MiXCR.${S}.${C}.log $ROOT/Scripts/MiXCR.sh
    done
done

# Extract UMI
for S in ${SAMPLES}; do
    echo ${S}
    export SAMPLE=${S}
    for C in ${CELLS}; do
        echo ${C}
        export CELL=${C}
        sbatch -n 1 -c 4 -t 3-0 --job-name="EXTRACT" --dependency=afterany:${DEMULTI%?} --output=$ROOT/Logs/Extract.UMI.${S}.${C}.log $ROOT/Scripts/ExtractUMI.sh
    done
done

EXTRACT=$(squeue -o "%A" -h -u dcoffey -n "EXTRACT" -S i | tr "\n" ":")

# STAR alignment to hg19
for S in ${SAMPLES}; do
    echo ${S}
    export SAMPLE=${S}
    for C in ${CELLS}; do
        echo ${C}
        export CELL=${C}
        export PREFIX="$ALIGNMENT_DIRECTORY/$SAMPLE.$CELL"
        export READ1="$UMI_EXTRACTED_DIRECTORY/$SAMPLE.$CELL.R1.UMI.fastq.gz"
    sbatch -n 1 -c 4 -t 3-0 --job-name="STAR" --dependency=afterany:${EXTRACT%?} --output=$ROOT/Logs/STAR.${S}.${C}.log $ROOT/Scripts/STAR_alignment.sh
    done
done

STAR=$(squeue -o "%A" -h -u dcoffey -n "STAR" -S i | tr "\n" ":")

# Deduplicate alignment and count genes
for S in ${SAMPLES}; do
    echo ${S}
    export SAMPLE=${S}
    for C in ${CELLS}; do
        echo ${C}
        export CELL=${C}
        export BAM="$ALIGNMENT_DIRECTORY/$SAMPLE.$CELL.Aligned.sortedByCoord.out.bam"
    sbatch -n 1 -c 4 -t 3-0 --job-name="DEDUP" --dependency=afterany:${STAR%?} --output=$ROOT/Logs/Deduplicate.${S}.${C}.log $ROOT/Scripts/Deduplicate.sh
    done
done

DEDUP=$(squeue -o "%A" -h -u dcoffey -n "DEDUP" -S i | tr "\n" ":")

# Calculate fastq file statistics
for S in ${SAMPLES}; do
    echo ${S}
    export SAMPLE=${S}
    for C in ${CELLS}; do
        echo ${C}
        export CELL=${C}
        export PREFIX1="$SAMPLE.$CELL.R1"
        export READ1="$UMI_EXTRACTED_DIRECTORY/$SAMPLE.$CELL.R1.UMI.fastq.gz"
        export PREFIX2="$SAMPLE.$CELL.R2"
        export READ2="$UMI_EXTRACTED_DIRECTORY/$SAMPLE.$CELL.R2.UMI.fastq.gz"
        sbatch -n 1 -c 4 -t 3-0 --job-name="FSTATS" --dependency=afterany:${EXTRACT%?} --output=$ROOT/Logs/FastqStats.${S}.${C}.log $ROOT/Scripts/FastqStats.sh
    done
done

FSTATS=$(squeue -o "%A" -h -u dcoffey -n "FSTATS" -S i | tr "\n" ":")

# Merge fastq stat files
module load R/3.4.2-foss-2016b-fh1
sbatch -n 1 -c 4 -t 3-0 --job-name="FMERGE" --output=$ROOT/Logs/FastqStatsMerge.log $ROOT/Scripts/MergeFastqStats.R

# Calculate bam file statistics
for S in ${SAMPLES}; do
    echo ${S}
    export SAMPLE=${S}
    for C in ${CELLS}; do
        echo ${C}
        export CELL=${C}
        export BAMSTATS_DIRECTORY="$ROOT/Aligned/Stats/Deduplicated"
        export PREFIX="$SAMPLE.$CELL"
        export BAM="$DEDUPLICATED_DIRECTORY/$SAMPLE.$CELL.deduplicated.featureCounts.bam"
        sbatch -n 1 -c 4 -t 3-0 --job-name="BSTATS1" --dependency=afterany:${DEDUP%?} --output=$ROOT/Logs/BamStats.${S}.${C}.log $ROOT/Scripts/BamStats.sh
    done
done

BSTATS1=$(squeue -o "%A" -h -u dcoffey -n "BSTATS1" -S i | tr "\n" ":")

# Calculate bam file statistics
sbatch -n 1 -c 4 -t 3-0 --job-name="BMERGE1" --dependency=afterany:${BSTATS1%?} --output=$ROOT/Logs/FastqStatsMerge.log $ROOT/Scripts/MergeBamStats.R

for S in ${SAMPLES}; do
    echo ${S}
    export SAMPLE=${S}
    for C in ${CELLS}; do
        echo ${C}
        export CELL=${C}
        export BAMSTATS_DIRECTORY="$ROOT/Aligned/Stats/Original"
        export PREFIX="$SAMPLE.$CELL"
        export BAM="$ALIGNMENT_DIRECTORY/$SAMPLE.$CELL.Aligned.sortedByCoord.out.bam"
        sbatch -n 1 -c 4 -t 3-0 --job-name="BSTATS2" --dependency=afterany:${STAR%?} --output=$ROOT/Logs/BamStats.${S}.${C}.log $ROOT/Scripts/BamStats.sh
    done
done

BSTATS2=$(squeue -o "%A" -h -u dcoffey -n "BSTATS2" -S i | tr "\n" ":")

# Merge bam stat files
sbatch -n 1 -c 4 -t 3-0 --job-name="BMERGE2" --dependency=afterany:${BSTATS2%?} --output=$ROOT/Logs/BamStatsMerge.log $ROOT/Scripts/MergeBamStats.R
