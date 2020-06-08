#!/bin/bash
#SBATCH -p long            # Partition to submit to
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu 8Gb     # Memory in MB
#SBATCH -J merge_BAM                # job name
#SBATCH -o logs/merge_BAM.%j.out    # File to which standard out will be written
#SBATCH -e logs/merge_BAM.%j.err    # File to which standard err will be written


#Júlia Perera
#10/2/2020

module purge
module load SAMtools/1.8-foss-2016b

line=$1
name=`basename $line`
OUTDIR=$2


#Definim els Lanes
LANE1=L001Aligned.sortedByCoord.out.bam
LANE2=L002Aligned.sortedByCoord.out.bam
LANE3=L003Aligned.sortedByCoord.out.bam
LANE4=L004Aligned.sortedByCoord.out.bam


#######################################################################
#Samtools's merge expects each input BAM to be sorted by algorithm X,  
#and efficiently produces an output file also sorted by algorithm X
#No cal fer sorting després del merge!


samtools merge $OUTDIR/$name.bam ${line}_${LANE1} ${line}_${LANE2} ${line}_${LANE3} ${line}_${LANE4}


#####################################################################################################
##################################### RNA METRICS ########################################################
ANNOTGENE=/bicoh/MARGenomics/AnalysisFiles/Annot_files_GTF
 
module purge  ## Why? Clear out .bashrc /.bash_profile settings that might interfere
module load picard/2.2.4-Java-1.8.0_92


# With FIRST_READ_TRANSCRIPTION_STRAND we get almost all (97%) reads classified as "INCORRECT STRAND READS"
# SECOND_READ_TRANSCRIPTION_STRAND (about 98% correct)
# Can try still NONE (no results)

java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics \
	I=$OUTDIR/$name.bam \
	REF_FLAT=$ANNOTGENE/gencode.v29.flatFile \
    	RIBOSOMAL_INTERVALS=$ANNOTGENE/gencode.v29.ribosomal.interval_list \
    	STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
	O=${OUTDIR}/${name}.RNA_Metrics \
	CHART=${OUTDIR}/${name}.RNA_Metrics.pdf
	


