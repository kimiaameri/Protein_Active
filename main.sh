#!/bin/sh
#SBATCH --time=100:00:00   # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=64gb     # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=main
#SBATCH --error=main.%J.err
#SBATCH --output=main.%J.out

export INTERSECTIONS_PATH="$WORK/Protein_Active-outputs/intersections/"
export OUTPUT_PATH="$WORK/Protein_Active-outputs/"
export GENOME_BED_PATH="$WORK/Protein_Active_reference_genome/"
export $SOURCE_DIR= "$WORK/Protein_Active/"


Rscript mappinghmmerhit.R $SOURCE_DIR $GENOME_BED_PATH $INTERSECTIONS_PATH 
