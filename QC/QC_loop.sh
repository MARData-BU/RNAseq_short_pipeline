#!/bin/bash
#SBATCH -p lowmem            # Partition to submit to
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu 6Gb     # Memory in MB
#SBATCH -J QC_loop           # job name
#SBATCH -o logs/QC_loop.%A_%a.log    # File to which standard out will be written
#SBATCH -e logs/QC_loop.%A_%a.log    # File to which standard err will be written

# QC Analysis of RNASeq samples of project: 

PROJECT=$1
BATCH=$2
# Prepare variables
#------------------

path=/bicoh/MARGenomics
DIR=${path}/${PROJECT}
FASTQDIR=${DIR}/rawData
suffix=_R1.fastq.gz

# Prepare folders
#------------------
mkdir $DIR/QC/logs

#============#
#   FASTQC   #
#============#
#get the number of files with fastq.gz extention:
length_files=$(ls -lR $FASTQDIR/${BATCH}/*.fastq.gz | wc -l)
#run the job (QC_loop) with array mode:
sbatch --array=1-$length_files $DIR/QC/fastqc.sh $PROJECT $BATCH


#=================#
#   FASTQSCREEN   #
#=================#
#get the number of files with fastq.gz extention:
length_files=$(ls -lR $FASTQDIR/${BATCH}/*.fastq.gz | wc -l)
#run the job (QC_loop) with array mode:
sbatch --array=1-$length_files $DIR/QC/fastq_screen.sh $PROJECT $BATCH


#=================#
#      FASTP      #
#=================#

#mkdir $DIR_out/QC/${RUN}/fastp
#OUTDIR=$DIR_out/QC/${RUN}/fastp

#for i in $(ls $FASTQDIR/*${suffix} | sed "s/$suffix//")
#	do
#	echo $i
#	sbatch $DIR/QC/fastp.sh $i $OUTDIR $suffix
#	sleep 1
#done



