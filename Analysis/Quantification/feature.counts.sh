#!/bin/bash
#SBATCH -p bigmem,normal            # Partition to submit to
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu 13Gb     # Memory in MB
#SBATCH -J feturecounts           # job name
#SBATCH -o logs/feturecounts.%j.out    # File to which standard out will be written
#SBATCH -e logs/feturecounts.%j.err    # File to which standard err will be written


#Júlia Perera
#Generate the table of counts with featurecounts
#5/8/2020

module purge
module load subread/1.6.4

#------------------------
# Prapare folders

PROJECT=$1
BATCH=$2

path=/bicoh/MARGenomics
DIR=${path}/${PROJECT}
BAMDIR=$DIR/Analysis/ReadMapping/BAM_Files/${BATCH}


mkdir $DIR/Analysis/Quantification/CountFiles/${BATCH}
OUTDIR=$DIR/Analysis/Quantification/CountFiles/${BATCH}


# Save in $INPUT a string of all input files separated by space
INPUT=`ls $BAMDIR/*.bam | paste -sd " " -`
echo $INPUT
###################################################################################
########################    COUNTS    #############################################
ANNOTGENE=/bicoh/MARGenomics/AnalysisFiles/Annot_files_GTF/Mouse

featureCounts -T $SLURM_CPUS_PER_TASK -p -t gene --largestOverlap -g gene_name -a $ANNOTGENE/gencode.vM28.primary_assembly.annotation.gtf -o $OUTDIR/CountsTable.txt $INPUT

###################################################################################
######################    MULTIQC     #############################################

#Copy CountsTable.txt.summary to multiqc folder
cp $DIR/Analysis/Quantification/CountFiles/${BATCH}/CountsTable.txt.summary $DIR/QC/${BATCH}/multiQC

cd $DIR/QC/${BATCH}
module load Python/3.5.2-foss-2016b

multiqc . -o $DIR/QC/${BATCH}/multiQC -f


