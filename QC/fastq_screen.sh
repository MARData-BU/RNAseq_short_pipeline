#!/bin/bash
#SBATCH -p long,normal,short            # Partition to submit to
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu 5Gb     # Memory in MB
#SBATCH -J FastQScreen           # job name
#SBATCH -o logs/FastQScreen.%j.out    # File to which standard out will be written
#SBATCH -e logs/FastQScreen.%j.err    # File to which standard err will be written

#-------------------------------

module purge  ## Why? Clear out .bashrc /.bash_profile settings that might interfere
module load fastq_screen/0.14.0
module load Bowtie2/2.3.5.1		# Required for Fastqscreen


#------------------------
# Prapare folders

PROJECT=$1
BATCH=$2

path=/bicoh/MARGenomics
DIR=${path}/${PROJECT}
FASTQDIR=$DIR/rawData/${BATCH}

mkdir -p $DIR/QC/${BATCH}/FastqScreen
OUTDIR=$DIR/QC/${BATCH}/FastqScreen

#--------------------
# Prapare input files

FASTQFILES=($(ls -1 $FASTQDIR/*.fastq.gz)) 
i=$(($SLURM_ARRAY_TASK_ID - 1)) ## bash arrays are 0-based
INFILE=${FASTQFILES[i]}


##############################################
#############  FASTQ SCREEN ##################

config=/bicoh/MARGenomics/AnalysisFiles/Index_Genomes_Bowtie2/fastq_screen.conf

fastq_screen --threads $SLURM_CPUS_PER_TASK --conf $config --outdir $OUTDIR $INFILE


