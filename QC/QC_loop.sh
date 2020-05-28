#!/bin/bash
#SBATCH -p lowmem            # Partition to submit to
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu 6Gb     # Memory in MB
#SBATCH -J QC_loop           # job name
#SBATCH -o QC_loop.%j.out    # File to which standard out will be written
#SBATCH -e QC_loop.%j.err    # File to which standard err will be written

# QC Analysis of RNASeq samples of project: 

PROJECT=20190826_AGimenez_FIS_UFred

# Prepare variables
#------------------

path=/bicoh/MARGenomics
DIR=${path}/${PROJECT}
FASTQDIR=${DIR}/rawData

# Prepare folders
#------------------
mkdir $DIR/QC/logs

#============#
#   FASTQC   #
#============#
mkdir $DIR/QC/FastQC

for i in $(ls $FASTQDIR/*fastq.gz) 
	do
	echo $i
	outdir=$DIR/QC/FastQC
	sbatch $DIR/QC/fastqc.sh $i $outdir
	sleep 1
done


#=================#
#   FASTQSCREEN   #
#=================#
mkdir $DIR/QC/FastqScreen

for i in $(ls $FASTQDIR/*fastq.gz) 
	do
	echo $i
	outdir=$DIR/QC/FastqScreen
	sbatch $DIR/QC/fastq_screen.sh $i $outdir
	sleep 1
done


