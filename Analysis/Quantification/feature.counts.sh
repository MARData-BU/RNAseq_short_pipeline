#!/bin/bash
#SBATCH -p bigmem            # Partition to submit to
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu 20Gb     # Memory in MB
#SBATCH -J counts           # job name
#SBATCH -o logs/counts.%j.out    # File to which standard out will be written
#SBATCH -e logs/counts.%j.err    # File to which standard err will be written


#Júlia Perera
#Generate the table of counts with featurecounts
#10/2/2020

module purge
module load subread/1.6.4

###################################################################################
############################COUNTS#################################################
BAMDIR=$1
OUTDIR=$2

ANNOTGENE=/bicoh/MARGenomics/AnalysisFiles/Annot_files_GTF/Human

featureCounts -T $SLURM_CPUS_PER_TASK -p -t exon --largestOverlap -g gene_name -a $ANNOTGENE/gencode.v29.primary_assembly.annotation.gtf -o $OUTDIR/CountsTable.txt \
$BAMDIR/0101-1_S1.bam \
$BAMDIR/0102-2_S2.bam \
$BAMDIR/0201-3_S3.bam \
$BAMDIR/0202-4_S4.bam \
$BAMDIR/0401-rep-11_S4.bam \
$BAMDIR/0402-rep-12_S5.bam \
$BAMDIR/0601-solar_S12.bam \
$BAMDIR/0602-solar_S13.bam \
$BAMDIR/0701-9_S5.bam \
$BAMDIR/0702-15_S6.bam \
$BAMDIR/0801-rep-21_S9.bam \
$BAMDIR/0802-rep-22_S10.bam \
$BAMDIR/0901-16_S2.bam \
$BAMDIR/0902-17_S3.bam \
$BAMDIR/1001-23_S8.bam \
$BAMDIR/1002-rep-24_S9.bam \
$BAMDIR/1201-27_S6.bam \
$BAMDIR/1202-28_S7.bam \
$BAMDIR/HC0101-rep-18_S13.bam \
$BAMDIR/HC01201_S11.bam \
$BAMDIR/HC01301_S12.bam \
$BAMDIR/HC0401-rep-20_S8.bam \
$BAMDIR/HC0501-10_S7.bam \
$BAMDIR/HC0601-8_S1.bam \
$BAMDIR/HC0801-rep-30_S10.bam \
$BAMDIR/HC0901-rep-31_S11.bam
