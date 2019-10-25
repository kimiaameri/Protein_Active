#!/bin/sh
#SBATCH --time=100:00:00   # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=64gb     # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=bed
#SBATCH --error=bed.%J.err
#SBATCH --output=bed.%J.out

#------------------bedfile
python3 GenomeBedPull.py $WORK/PA_reference_genome 
export GENOME_BED_PATH="$WORK/PA_reference_genome/"
python3 pythonVcfbed.py ./InputFiles.csv $MINICONDA_HOME
sh vcfBed.sh
python3 pythonIntersections.py ./InputFiles.csv $GENOME_BED_PATH $MINICONDA_HOME
sh mapVCF-to-Bed.sh
export INTERSECTIONS_PATH="$WORK/PA-outputs/intersection/"
export OUTPUT_PATH="$WORK/PA-outputs/"


cd $WORK/PA-outputs
export SOURCE_DIR="$WORK/PA"

cd $WORK/PA_reference_genome
cat  nctc8325.bed | tail -n+2 > nctc8325-1.bed 
cd $WORK/PA/  
Rscript mappinghmmerhit.R $SOURCE_DIR $GENOME_BED_PATH $INTERSECTIONS_PATH ./InputFiles.csv bigtable.csv tableWeight.csv 



Rscript maincode.R $SOURCE_DIR $GENOME_BED_PATH $INTERSECTIONS_PATH ./InputFiles.csv bigtable.csv tableWeight.csv 
