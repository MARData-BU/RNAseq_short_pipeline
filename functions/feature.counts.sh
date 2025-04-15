#!/bin/bash
#SBATCH -p bigmem,normal            # Partition to submit to
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu 13Gb     # Memory in MB
#SBATCH -J feturecounts           # job name
#SBATCH -o logs/feturecounts.%j.out    # File to which standard out will be written
#SBATCH -e logs/feturecounts.%j.err    # File to which standard err will be written


module purge
module load subread/2.0.34

#------------------------
# Prapare folders

DIR=$1
END=$2
STRAND=$3

BAMDIR=$DIR/02_ReadMapping/BAM_Files


mkdir -p $DIR/03_Quantification/CountFiles
OUTDIR=$DIR/03_Quantification/CountFiles


# Save in $INPUT a string of all input files separated by space
INPUT=`ls $BAMDIR/*.bam | paste -sd " " -`
echo $INPUT
###################################################################################
########################    COUNTS    #############################################
ANNOTGENE=/bicoh/MARGenomics/AnalysisFiles/Annot_files_GTF/Human/gencode.v41.primary_assembly.annotation.gtf


if [ $END == PAIRED ]
	then
	# Paired end
	  featureCounts -T $SLURM_CPUS_PER_TASK -s $STRAND -p -t exon --countReadPairs --largestOverlap -g gene_name -a $ANNOTGENE -o $OUTDIR/CountsTable.txt $INPUT
elif [ $END == SINGLE ]
	then
	# NOT paired end
	  featureCounts -T $SLURM_CPUS_PER_TASK -s $STRAND -t exon --largestOverlap -g gene_name -a $ANNOTGENE -o $OUTDIR/CountsTable.txt $INPUT
	fi

