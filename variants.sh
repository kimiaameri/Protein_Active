Rscript filename.R ../Protein_Active-outputs/VariantPosition/ ./SRA.txt

 for x in `cat SRA.txt`; do echo $x; sbatch mapDomian.sh $x; done
