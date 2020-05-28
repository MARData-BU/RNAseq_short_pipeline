#!/bin/bash

## INPUT
project=$1

## BUILD VARS AND CREATE FOLDERS
PROJECT=$project
path=/bicoh/MARGenomics
DIR=${path}/${PROJECT}

mkdir $DIR/Analysis/ReadMapping/
mkdir $DIR/Analysis/ReadMapping/BAM_Files
FASTQDIR=${DIR}/rawData
OUTDIR=$DIR/Analysis/ReadMapping/BAM_Files

## SEND ALIGNMENT JOBS
for i in $(ls $FASTQDIR/*_read1.fastq.gz | sed 's/_read1.fastq.gz//')
	do
	echo $i
	sbatch $DIR/Analysis/ReadMapping/STAR.alignment.sh $i $OUTDIR
	sleep 1
done
