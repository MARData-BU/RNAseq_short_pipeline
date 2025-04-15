#!/bin/bash

PROJECT="project_name"

# Prepare variables
#------------------

path=/bicoh/MARGenomics
DIR=${path}/${PROJECT}


# Prepare folders 
#------------------
FASTQDIR=${DIR}/rawData #folder with input fastq files
mkdir $DIR/logs #all stdout and stderr will be stored here

mkdir $DIR/01_QC
mkdir $DIR/02_ReadMapping
mkdir $DIR/03_Quantification


# Prepare to start pipeline 
#------------------
cd $DIR # Move to wd


#################################### 01 Quality check:  #################################################

#============#
#   FASTQC   #
#============#

#get the number of files with fastq.gz extention:
length_files=$(ls -lR $FASTQDIR/*.fastq.gz | wc -l)
#run the job with array mode:
sbatch --array=1-$length_files $DIR/functions/fastqc.sh $DIR



#################################### 02 Read Alignement:  #################################################

#============#
#   STAR   #
#============#

# Prepare variables
#------------------
R1=_R1.fastq.gz #specify suffix R1 of the fastq files
R2=_R2.fastq.gz #specify suffix R2 of the fastq files (leave empty if single)
END=PAIRED #or SINGLE


#get the number of files with $R1 extention (it has to run once per sample, not per file. Important in case of paired data):
length_files=$(ls -1 $DIR/rawData/*$R1| wc -l)
#run the job with array mode:
sbatch --array=1-$length_files $DIR/functions/STAR.alignment.sh .sh $DIR $END $R1 $R2

#################################### 03 Quantification:  #################################################

#====================#
#   FeatureCouns MM  #
#====================#
#run the job just once, as featurecounts runs for multiple input BAM files:
sbatch $DIR/functions/feature.counts.sh $DIR



#################################### MultiQC:  #################################################

#=================#
#   MultiQC       #
#=================#

module load Python/3.5.2-foss-2016b

multiqc . -o $DIR/01_QC -f




