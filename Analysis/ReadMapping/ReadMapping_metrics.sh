#!/bin/bash
#SBATCH -p normal,long            # Partition to submit to
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu 18Gb     # Memory in MB
#SBATCH -J ReadMapping_metrics          # job name
#SBATCH -o logs/ReadMapping_metrics.%j.out    # File to which standard out will be written
#SBATCH -e logs/ReadMapping_metrics.%j.err    # File to which standard err will be written

# Prepare variables
#------------------
## INPUT
PROJECT=$1
BATCH=$2

path=/bicoh/MARGenomics
DIR=${path}/${PROJECT}
QC=${DIR}/QC/${BATCH}


# Load modules
#------------------
module load Python/3.5.2-foss-2016b
module load R/4.0.0


# Prepare folders
#------------------
mkdir $QC/STAR

# Move to QC folder
cd $QC

#------------------

#=================#
#  Manual summary
#=================#

OUTDIR=$DIR/Analysis/ReadMapping/BAM_Files
touch $OUTDIR/TotalCounts_Alignment

for i in $OUTDIR/*.final.out; do basename $i >> $OUTDIR/TotalCounts_Alignment; \
grep "Number of input reads" "$i" >> $OUTDIR/TotalCounts_Alignment; grep "Uniquely mapped reads" "$i"\
 >> $OUTDIR/TotalCounts_Alignment; grep "Average mapped length" "$i" >> $OUTDIR/TotalCounts_Alignment;\
 grep "reads mapped to too many loci" "$i" >> $OUTDIR/TotalCounts_Alignment; grep "too many mismatches" "$i"\
 >> $OUTDIR/TotalCounts_Alignment; grep "too short" "$i" >> $OUTDIR/TotalCounts_Alignment;\
 grep "other" "$i" >> $OUTDIR/TotalCounts_Alignment; done


#=================#
#   MultiQC       #
#=================#

# Copy **Log.final.out to multiqc folder
cp $DIR/Analysis/ReadMapping/BAM_Files/${BATCH}/*Log.final.out $QC/STAR
cp $DIR/Analysis/ReadMapping/BAM_Files/${BATCH}/*RNA_Metrics $QC/STAR

cd $QC
multiqc . -f -o $QC/multiQC


#=================================#
# Create ReadMapping QC table     #
#=================================#
LANES=1
ReadMapping=$DIR/Analysis/ReadMapping
starDir=$DIR/Analysis/ReadMapping/BAM_Files/${BATCH}
rnametricsDir=$starDir

Rscript $DIR/Analysis/ReadMapping/ReadMappingQC.R $ReadMapping $LANES $starDir $rnametricsDir 

#====================================#
# Add new RUN metrics to BATCH excel #
#====================================#

#Rscript $DIR/QC/excelQC_addRUNtoBATCH.R ${DIR}/QC/${BATCH} $BATCH $RUN


