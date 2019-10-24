#!/bin/sh
#SBATCH --time=100:00:00   # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=64gb     # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=hmmrsearch
#SBATCH --error=hmmrsearch.%J.err
#SBATCH --output=hmmrsearch.%J.out
/hmmer/hmmer-3.2.1/easel/miniapps/
make
cd $WORK/independent2
../hmmer/hmmer-3.2.1/src/hmmsearch Pfam-A.hmm trans.fasta > ./hmmerfile/s.Areuse.txt
../hmmer/hmmer-3.2.1/src/hmmsearch --tblout ./hmmerfile/S.Areuse.tbl Pfam-A.hmm trans.fasta


../hmmer/hmmer-3.2.1/src/hmmsearch -A ./hmmerfile/S.Areuse.sto Pfam-A.hmm trans.fasta 
../hmmer/hmmer-3.2.1/easel/miniapps/esl-reformat fasta ./hmmerfile/S.Areuse.sto > ./hmmerfile/S.Areuse.fa
../hmmer/hmmer-3.2.1/easel/miniapps/esl-sfetch --index trans.fasta 

../hmmer/hmmer-3.2.1/src/hmmsearch --domtblout ./hmmerfile/S.Areuse.dtbl Pfam-A.hmm trans.fasta
cd hmmerfile
cat S.Areuse.dtbl |grep -v "^#" | awk '{print $1" "$4" "$18" "$19}'  > S.areuse.hit.fa
sort S.areuse.hit.fa > sort.hit.S.areus.csv
