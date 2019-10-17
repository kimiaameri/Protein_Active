#!/bin/sh
#SBATCH --time=100:00:00   # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=64gb     # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=hmmrsearch
#SBATCH --error=hmmrsearch.%J.err
#SBATCH --output=hmmrsearch.%J.out
/hmmer/hmmer-3.2.1/easel/miniapps/
make
cd $WORK/independent2
#../hmmer/hmmer-3.2.1/src/hmmsearch Pfam-A.hmm trans.txt > ./hmmerfile/s.Areuse.txt


../hmmer/hmmer-3.2.1/src/hmmsearch -A S.Areuse.sto Pfam-A.hmm trans.fasta 
../hmmer/hmmer-3.2.1/easel/miniapps/esl-reformat fasta S.Areuse.sto > S.Areuse.fa
../hmmer/hmmer-3.2.1/easel/esl-sftech --index trans.fasta 
