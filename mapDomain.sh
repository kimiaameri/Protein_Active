#!/bin/sh
#SBATCH --time=100:00:00   # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=64gb     # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=mapdomain$1
#SBATCH --error=mapdomain$1.%J.err
#SBATCH --output=mapdomain$1.%J.out

export INTERSECTIONS_PATH="$WORK/Protein_Active-outputs/intersections/"
export OUTPUT_PATH="$WORK/Protein_Active-outputs/"
export GENOME_BED_PATH="$WORK/Protein_Active_reference_genome/"
export SOURCE_DIR="$WORK/Protein_Active/"
export VARIANT_PATH="$WORK/Protein_Active-outputs/VariantPosition/"
export DOMAIN_ISOLATES="$WORK/Protein_Active-outputs/DomainIsolates/"

Rscript HammrMapping.R $SOURCE_DIR $GENOME_BED_PATH $INTERSECTIONS_PATH $OUTPUT_PATH $VARIANT_PATH/$1 $DOMAIN_ISOLATES/$1 $DOMAIN_ISOLATES/$1.norm




