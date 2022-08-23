#!/bin/bash
#SBATCH -p short            # Partition to submit to
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu 4Gb     # Memory in MB
#SBATCH -J FastQC           # job name
#SBATCH -o logs/FastQC.%A_%a.out    # File to which standard out will be written
#SBATCH -e logs/FastQC.%A_%a.err    # File to which standard err will be written



module purge 
module load FastQC/0.11.5-Java-1.7.0_80        

#------------------------
# Prapare folders

PROJECT=$1
BATCH=$2

path=/bicoh/MARGenomics
DIR=${path}/${PROJECT}
FASTQDIR=$DIR/rawData/${BATCH}

mkdir -p $DIR/QC/${BATCH}/FastQC
OUTDIR=$DIR/QC/${BATCH}/FastQC

#--------------------
# Prapare input files

FASTQFILES=($(ls -1 $FASTQDIR/*.fastq.gz)) 
i=$(($SLURM_ARRAY_TASK_ID - 1)) ## bash arrays are 0-based
INFILE=${FASTQFILES[i]}

##############################################
#############  FASTQC ########################

fastqc --outdir $OUTDIR --threads $SLURM_CPUS_PER_TASK $INFILE

