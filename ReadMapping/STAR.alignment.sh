#!/bin/bash
#SBATCH -p normal            # Partition to submit to
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu 8Gb           # Memory in MB
#SBATCH -J STAR           # job name
#SBATCH -o logs/STAR.%j.out    # File to which standard out will be written
#SBATCH -e logs/STAR.%j.err    # File to which standard err will be written

module load STAR/2.6.0a


ANNOTGENE=/bicoh/MARGenomics/Analysis_Files/Annot_files_GTF
GNMIDX=/bicoh/MARGenomics/Analysis_Files/Index_Genomes_STAR/Idx_Gencode_hg38_readlength75

lane=$1
OUTDIR=$2
name=`basename $lane`

R1=_R1_001.fastq.gz
R2=_R2_001.fastq.gz
######################################################################################################
#####################################ALIGNMENT########################################################

STAR --runThreadN $SLURM_CPUS_PER_TASK\
 --genomeDir $GNMIDX --readFilesIn $lane$R1 $lane$R2 --readFilesCommand zcat --outFileNamePrefix\
 $OUTDIR/$name --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outFilterType BySJout\
 --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999\
 --outFilterMismatchNoverLmax 0.05 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
 --sjdbGTFfile $ANNOTGENE/gencode.v29.primary_assembly.annotation.gtf
