# Protein_Active
### The purpose of this project is to identify active protein domains in **S***taphylococcus* **A***ureus* isolates that are resistance or suseptible to Gentimicine.
##### for this purpose, first 363 Staphyloccuse Areuse isolates are downloaded from NCBI. The selected dataset has 174 resistant isolate and 189 suspetiple 
-----------------------------------------------------------------
 ####
 list of SRA is availbel at ![SRA list](https://github.com/kimiaameri/Protein_Active/blob/master/SRA_list.txt)
 - to downlaod data from NCBI the latset version of SRAtoolkit from Anaconda
 ## Download SRA from NCBI
```bash
for x in `cat SRA.txt`; do 
fastq-dump --split-files $x ;  
sleep 5;
done
```
This bash code can also run as:
```bash
sbatch SRA.sh  
```
####
----------------------------------------------------------------- 
## Pre-processing and variant analysis of the whole genome sequences. 
#### For variant analysis we built the "snpvariant environment" for bioinformatic tools we usein this step:

## preprocessing step by step:
```bash
~/miniconda3/envs/snpvariant/bin/bwa mem $WORK/reference_genome/Staphylococcus_aureus_NCTC_8325/NCBI/2006-02-13/Sequence/BWAIndex/genome.fa $WORK/outputs/trimmomatic/SrA.Accession-R1.paired.fq $WORK/outputs/trimmomatic/SrA.Accession-R2.paired.fq >$WORK/outputs/samfiles/SrA.Accession.sam
~/miniconda3/envs/snpvariant/bin/samtools view -bt $WORK/reference_genome/Staphylococcus_aureus_NCTC_8325/NCBI/2006-02-13/Sequence/BWAIndex/genome.fa $WORK/outputs/samfiles/SrA.Accession.sam >$WORK/outputs/bamfiles/SrA.Accession.bam
~/miniconda3/envs/snpvariant/bin/samtools flagstat $WORK/outputs/bamfiles/SrA.Accession.bam > $WORK/outputs/flagsam/SrA.Accession.flagstat.log
~/miniconda3/envs/snpvariant/bin/samtools sort $WORK/outputs/bamfiles/SrA.Accession.bam -O bam -o $WORK/outputs/sortsam/SrA.Accession.sorted.bam
~/miniconda3/envs/snpvariant/bin/samtools stats $WORK/outputs/sortsam/SrA.Accession.sorted.bam >$WORK/outputs/stats/SrA.Accession.txt 
~/miniconda3/envs/snpvariant/bin/samtools depth -a $WORK/outputs/sortsam/SrA.Accession.sorted.bam > $WORK/outputs/depth/SrA.Accession.depth
~/miniconda3/envs/snpvariant/bin/picard MarkDuplicates I=$WORK/outputs/sortsam/SrA.Accession.sorted.bam O=$WORK/outputs/picard/SrA.Accession.picard.bam M=$WORK/outputs/picard/picardlog/SrA.Accession.picard.log
~/miniconda3/envs/snpvariant/bin/freebayes -f $WORK/reference_genome/Staphylococcus_aureus_NCTC_8325/NCBI/2006-02-13/Sequence/WholeGenomeFa
~/miniconda3/envs/snpvariant/bin/vcffilter -f "DP > $DEPTH" $WORK/outputs/freebayesoutput/SrA.Accession.vcf > $WORK/outputs/vcffilter-dp/SrA.Accession.vcf
~/miniconda3/envs/snpvariant/bin/bcftools view -Ob $WORK/outputs/vcffilter-dp/SrA.Accession.vcf > $WORK/outputs/bcfoutput/SrA.Accession.vcf.gz
~/miniconda3/envs/snpvariant/bin/bcftools index $WORK/outputs/bcfoutput/SrA.Accession.vcf.gz

```
