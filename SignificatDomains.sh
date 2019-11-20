#!/bin/sh
#SBATCH --time=80:00:00   # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=16384     # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=SAVEA
#SBATCH --error=SAVEA.%J.err
#SBATCH --output=SAVEA.%J.out


export SOURCE_PATH="$WORK/Protein_Active/Source/"
export DOMAIN_ISOLATES="$WORK/Protein_Active-outputs/DomainIsolates/"

export OUTPUT_PATH="$WORK/Protein_Active-outputs/"


Rscript ./Source/significatDomains.R $SOURCE_PATH $DOMAIN_ISOLATES $OUTPUT_PATH/significantDomain.csv
cp $OUTPUT_PATH/significantDomain.csv ./Outputs/significantDomain.csv
