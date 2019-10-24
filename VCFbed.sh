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
