# RNASeq pipeline

The most simple RNASeq pipeline that includes the 3 basic steps:

1. QC (FastQC)
2. Alignment with STAR
3. Quantification with FeatureCounts

Scripts are prepared to run in a cluster with SLURM workload manager. and software manged with modules. In particular, it is optimized for the [GRIB](https://grib.upf.edu/) CPU cluster. 

STAR genome index can be provided upon request (mardata-bu@researchmar.net). However, the code to generate the index is included in the scripts.
