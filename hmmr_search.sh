#!/bin/sh
#SBATCH --time=100:00:00   # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=64gb     # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=hmmrsearch
#SBATCH --error=hmmrsearch.%J.err
#SBATCH --output=hmmrsearch.%J.out
cd $WORK/independent2
 ../hmmer/hmmer-3.2.1/src/hmmsearch Pfam-A.hmm trans.txt > ./hmmerfile/s.Areuse.txt
