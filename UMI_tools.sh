#!/bin/bash
# UMI extraction
# Written by David Coffey dcoffey@fhcrc.org
# Updated December 22, 2017

module load Python/3.6.2-foss-2016b-fh2
module load SAMtools/1.6-foss-2016b
module load STAR/2.5.0a-foss-2015b

umi_tools whitelist \
--stdin Plate-1_S1_L001.trimmed.R2.fastq.gz \
--bc-pattern='(?P<cell_1>.{33})(?P<umi_1>.{12})' \
--extract-method=regex \
--set-cell-number=94 \
--plot-prefix=$SAMPLE \
--plot-prefix \
--log2stderr > whitelist.txt

umi_tools extract \
--bc-pattern='(?P<cell_1>.{33})(?P<umi_1>.{12})' \
--extract-method=regex \
--stdin Plate-1_S1_L001.trimmed.R2.fastq.gz \
--stdout Plate-1_S1_L001.extracted.R2.fastq.gz \
--read2-in Plate-1_S1_L001.trimmed.R1.fastq.gz \
--read2-out=Plate-1_S1_L001.extracted.R1.fastq.gz \
--filter-cell-barcode \
--whitelist=whitelist.txt \
--log=extractumi.log

STAR --runThreadN 4 \
--genomeDir $GENOME_DIRECTORY \
--readFilesIn Plate-1_S1_L001.extracted.R2.fastq.gz \
--readFilesCommand zcat \
--outFilterMultimapNmax 1 \
--outSAMtype BAM SortedByCoordinate

featureCounts -a geneset.gtf -o gene_assigned -R BAM Aligned.sortedByCoord.out.bam -T 4
samtools sort Aligned.sortedByCoord.out.bam.featureCounts.bam -o assigned_sorted.bam
samtools index assigned_sorted.bam

umi_tools count --per-gene --gene-tag=XT --per-cell -I assigned_sorted.bam -S counts.tsv.gz