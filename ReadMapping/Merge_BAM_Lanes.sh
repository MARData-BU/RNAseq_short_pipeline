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


