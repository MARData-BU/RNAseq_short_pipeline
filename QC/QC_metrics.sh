#!/bin/bash
#SBATCH -p lowmem            # Partition to submit to
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu 6Gb     # Memory in MB
#SBATCH -J QC_metrics          # job name
#SBATCH -o QC_metrics.%j.out    # File to which standard out will be written
#SBATCH -e QC_metrics.%j.err    # File to which standard err will be written

# QC Analysis of RNASeq samples of project: 

PROJECT=20190826_AGimenez_FIS_UFred

# Prepare variables
#------------------

path=/bicoh/MARGenomics
DIR=${path}/${PROJECT}
QC=$DIR/QC

# Load modules
#------------------
module load Python/3.5.2-foss-2016b
module load R/4.0.0

# Prepare folders
#------------------
mkdir $DIR/QC/multiQC

# Move to QC folder
cd $QC
#------------------

#=================#
#   MultiQC       #
#=================#
multiqc . -f -o $QC/multiQC


#=================================#
# Create table 4 QC presentation  #
#=================================#
LANES=4
R=2 # paired end
COLOR=FALSE #color for duplications

Rscript $DIR/QC/table4QCpresentation.R $QC $LANES $R $COLOR # QC dir, numer of lanes, paired end, color for duplications