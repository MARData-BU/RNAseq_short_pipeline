#!/bin/bash
#SBATCH -p normal,long,bigmem   # Partition to submit to
#SBATCH --cpus-per-task=10	#change to 5 if no rush or cluster is full
#SBATCH --mem-per-cpu 9Gb      # Memory in MB
#SBATCH -J STAR                # job name
#SBATCH -o logs/STAR.%A_%a.out    		# Sdt out file name: %j SLURM_JOB_ID; %A SLURM_ARRAY_JOB_ID, %a SLURM_ARRAY_TASK_ID
#SBATCH -e logs/STAR.%A_%a.err    		# Sdt err file name: %j SLURM_JOB_ID; %A SLURM_ARRAY_JOB_ID, %a SLURM_ARRAY_TASK_ID

module load STAR/2.7.1a-foss-2016b

#------------------------
# Prapare folders

PROJECT=$1
BATCH=$2
suffix=$3

path=/bicoh/MARGenomics
DIR=${path}/${PROJECT}
FASTQDIR=$DIR/rawData/${BATCH}


mkdir -p $DIR/Analysis/ReadMapping/BAM_Files/${BATCH}
OUTDIR=$DIR/Analysis/ReadMapping/BAM_Files/${BATCH}

#--------------------
# Prapare input files

FASTQFILES=($(ls -1 $FASTQDIR/*${suffix} | sed "s/$suffix//")) 
i=$(($SLURM_ARRAY_TASK_ID - 1)) ## bash arrays are 0-based
INFILE=${FASTQFILES[i]}

name=`basename $INFILE`

R1=$suffix
R2=`echo $suffix | sed "s/R1/R2/"`

#--------------------
# Paths to index files
ANNOTGENE=/bicoh/MARGenomics/AnalysisFiles/Annot_files_GTF/Mouse # Mouse ---> For human use: ~/Human
GNMIDX=/bicoh/MARGenomics/AnalysisFiles/Index_Genomes_STAR/Idx_Gencode_v28_mm10_readlength75 # Mouse ---> For human use: ~/Idx_Gencode_v36_hg38_readlength75

######################################################################################################
##################################### ALIGNMENT ########################################################

STAR --runThreadN $SLURM_CPUS_PER_TASK\
 --genomeDir $GNMIDX --readFilesIn $INFILE$R1 $INFILE$R2 --readFilesCommand zcat --outFileNamePrefix\
 $OUTDIR/$name --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outFilterType BySJout\
 --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999\
 --outFilterMismatchNoverLmax 0.05 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
 --sjdbGTFfile $ANNOTGENE/gencode.vM28.primary_assembly.annotation.gtf # Mouse ---> For human use: gencode.v36.primary_assembly.annotation.gtf

 ######################################################################################################
##################################### Create index (.bai)###############################################

module purge
module load SAMtools/1.8-foss-2016b


samtools index ${OUTDIR}/${name}Aligned.sortedByCoord.out.bam ${OUTDIR}/${name}Aligned.sortedByCoord.out.bai
 

######################################################################################################
##################################### RNA METRICS ########################################################
 
module purge  ## Why? Clear out .bashrc /.bash_profile settings that might interfere
module load picard/2.2.4-Java-1.8.0_92

# With FIRST_READ_TRANSCRIPTION_STRAND we get almost all (97%) reads classified as "INCORRECT STRAND READS"
# SECOND_READ_TRANSCRIPTION_STRAND (about 98% correct)
# Can try still NONE (no results)

java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics \
	I=${OUTDIR}/${name}Aligned.sortedByCoord.out.bam \
	REF_FLAT=$ANNOTGENE/gencode.vM28.flatFile \
    	RIBOSOMAL_INTERVALS=$ANNOTGENE/gencode.vM28.ribosomal.interval_list \
    	STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
	O=${OUTDIR}/${name}.RNA_Metrics 

# Mouse ---> For human use v36 instead of vM28
