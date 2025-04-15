#!/bin/bash
#SBATCH -p normal,long,bigmem   # Partition to submit to
#SBATCH --cpus-per-task=10	#change to 5 if no rush or cluster is full
#SBATCH --mem-per-cpu 9Gb      # Memory in MB
#SBATCH -J STAR                # job name
#SBATCH -o logs/STAR.%A_%a.out    		# Sdt out file name: %j SLURM_JOB_ID; %A SLURM_ARRAY_JOB_ID, %a SLURM_ARRAY_TASK_ID
#SBATCH -e logs/STAR.%A_%a.err    		# Sdt err file name: %j SLURM_JOB_ID; %A SLURM_ARRAY_JOB_ID, %a SLURM_ARRAY_TASK_ID

module load STAR/2.7.8a-GCC-10.2.0

#------------------------
# Prapare folders

DIR=$1
END=$2
R1=$3
R2=$4
FASTQDIR=$DIR/rawData


mkdir -p $DIR/02_ReadMapping/BAM_Files
OUTDIR=$DIR/02_ReadMapping/BAM_Files

#--------------------
# Prapare input files

FASTQFILES=($(ls -1 $FASTQDIR/*${R1} | sed "s/$R1//")) 
i=$(($SLURM_ARRAY_TASK_ID - 1)) ## bash arrays are 0-based
INFILE=${FASTQFILES[i]}

name=`basename $INFILE`


#--------------------
# Paths to index files
ANNOTGENE=/bicoh/MARGenomics/AnalysisFiles/Annot_files_GTF/Human/gencode.v41.primary_assembly.annotation.gtf # For mouse use: ~/Mouse/gencode.vM33.primary_assembly.annotation.gtf
GNMIDX=/bicoh/MARGenomics/AnalysisFiles/Index_Genomes_STAR/Idx_Gencode_v41_hg38_readlength75 # For mouse use: Idx_Gencode_v33_mm10_readlength75

######################################################################################################
##################################### ALIGNMENT ########################################################

if [ $END == SINGLE ]
	then

	STAR --runThreadN $SLURM_CPUS_PER_TASK\
	 --genomeDir $GNMIDX --readFilesIn $INFILE$R1 --readFilesCommand zcat --outFileNamePrefix\
	 $OUTDIR/$name --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outFilterType BySJout\
	 --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999\
	 --outFilterMismatchNoverLmax 0.05 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
	 --sjdbGTFfile $ANNOTGENE


elif [ $END == PAIRED ]
	then
	
	STAR --runThreadN $SLURM_CPUS_PER_TASK\
	 --genomeDir $GNMIDX --readFilesIn $INFILE$R1 $INFILE$R2 --readFilesCommand zcat --outFileNamePrefix\
	 $OUTDIR/$name --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outFilterType BySJout\
	 --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999\
	 --outFilterMismatchNoverLmax 0.05 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
	 --sjdbGTFfile $ANNOTGENE
	 
	fi
 ######################################################################################################
##################################### Create index (.bai)###############################################

module purge
module load SAMtools/1.12-GCC-10.2.0


samtools index ${OUTDIR}/${name}Aligned.sortedByCoord.out.bam ${OUTDIR}/${name}Aligned.sortedByCoord.out.bai
 

