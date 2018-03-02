# sqRNAseq Pipeline v1.0
#### Updated March 1, 2018
These are a series of shell and R scripts used to process single cell RNAseq files using the QIAseq UPX 3’ transcriptome library construction kit.  To use this workflow, modify the environmental variables defined at the header of the scRNAseq_Pipeline.sh script.  These scripts were written to run on a cluster using [slurm work load manager](https://slurm.schedmd.com).  Our cluster is designed to load software using the `module load` command.  Since this is unique to our server, these lines of code can be removed and you will want to make sure that the direct path to the software has been specified or and environmental variable pointing to that software program is created.

The general sequence that the files are run in are as follows:

1. Demultiplex.sh
2. MiXCR.sh
3. ExtractUMI.sh
4. STAR_alignment.sh
5. Deduplicate.sh
6. FastqStats.sh
7. MergeFastqStats.R
8. BamStats.sh
9. MergeBamStats.R

