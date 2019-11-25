# Protein_Active
### The purpose of this project is to identify active protein domains in **S***taphylococcus* **A***ureus* isolates that are resistance or suseptible to Gentimicine.
##### for this purpose, first 363 Staphyloccuse Areuse isolates are downloaded from NCBI. The selected dataset has 174 resistant isolate and 189 suspetiple 
-----------------------------------------------------------------
 ####
 list of SRA is available at ![SrA_Accession.txt](https://github.com/kimiaameri/Protein_Active/blob/master/Inputs/SrA.Accession.txt)
 - to downlaod fastq files from NCBI the latset version of SRAtoolkit is downloaded from Anaconda.
 ## Download SRA from NCBI
```bash
for x in `cat SrA.Accession.txt`; do 
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
#### For variant analysis we built the "snpvariant environment" for bioinformatic tools we use for this step:

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
~/miniconda3/envs/snpvariant/bin/snpEff -v Staphylococcus_aureus_subsp_aureus_nctc_8325 $WORK/SNP-outputs/vcffilter-dp/SrA.Accession.vcf > $WORK/SNP-outputs/snpEff/SrA.Accession.ann.vcf 

```
### This code is genereted by ![pythonVariantAnalysis.py](https://github.com/kimiaameri/Protein_Active/tree/master/Source) Python code for each isolate 
------------------------------------------------------------------------------------------------------
## Hmmr Search
##### For search in Hmmr database, we need the sequence reads of proteins in ***Staphylococcus*** ***Areues***. This list is gereneretaed by ![translation.R](https://github.com/kimiaameri/Protein_Active/blob/master/Source/translation.R) R code.

```R
gb = readGenBank("./Inputs/sequence.gbk")
tr <- transcripts(gb)
proteinlist<-tr$gene_id
l<- length(proteinlist)
translations<- matrix(NA,ncol = 2,nrow = l)

for(i in 1:l)
{
      seq<- as.character(proteins$translation[i,2])
      translations[i,1]<- paste0(">",proteinlist[i])
      translations[i,2]<- seq
  print(i)
}
colnames(translations)<- c("Gene_id","translation")
write.table(x=translations,file="./Outputs/trans.fasta",quote =FALSE,sep="\n",eol = "\n", row.names = FALSE, col.names = FALSE)
```
#### Hmmr v.3.2.1
```bash
../hmmer/hmmer-3.2.1/src/hmmsearch Pfam-A.hmm trans.fasta > ./hmmerfile/s.Areuse.txt
../hmmer/hmmer-3.2.1/src/hmmsearch --tblout ./hmmerfile/S.Areuse.tbl Pfam-A.hmm trans.fasta

../hmmer/hmmer-3.2.1/src/hmmsearch -A ./hmmerfile/S.Areuse.sto Pfam-A.hmm trans.fasta 
../hmmer/hmmer-3.2.1/easel/miniapps/esl-reformat fasta ./hmmerfile/S.Areuse.sto > ./hmmerfile/S.Areuse.fa
../hmmer/hmmer-3.2.1/easel/miniapps/esl-sfetch --index trans.fasta 

../hmmer/hmmer-3.2.1/src/hmmsearch --domtblout ./hmmerfile/S.Areuse.dtbl Pfam-A.hmm trans.fasta
cd hmmerfile
cat S.Areuse.dtbl |grep -v "^#" | awk '{print $1" "$4" "$18" "$19}'  > S.areuse.hit.fa
sort S.areuse.hit.fa > sort.hit.S.areus.csv
```
------------------------------------------------------------------------------------------------------
#### Read hmmrfiles 
```R
hmmrHit <- read.table(paste0(outputPath,"hmmrfile/sort.hit.S.areus.csv"),header=F,sep=" ",stringsAsFactors = F)
length.hmmrHit <-  nrow(hmmrHit)
intersections<- list.files(intersectionspath)
Domain.length<- as.numeric(hmmrHit[,4]) - as.numeric(hmmrHit[,3])
hmmrHit<-cbind(hmmrHit,Domain.length)
colnames(hmmrHit)<- c("Gene.Id","Domain.Name","Start.DomPos","End.DomPos","Domain.Length")
uniq.Domains<- length(unique(hmmrHit[,2]))
```
------------------------------------------------------------------------------------------------------

### Find mutation position in each gene for isolates

```R
length.intersection= nrow(intersections)
    variant.matrix<- matrix(NA, ncol=5, nrow=length.intersection)
    colnames(variant.matrix)<- c("Gene.Id","Variant.start","Variant.end","Gene.length","Chromosome.Length")
      for (k in 1:length.intersection) 
        for (j in 1 :length.genome)
         if (intersections[k,2] >= reference_Genome[j,2] & intersections[k,2] <= reference_Genome[j,3]) 
        { 
           variant.matrix[k,1] = as.character(reference_Genome [j,4])
           z<-round(abs(intersections[k,2] - reference_Genome [j,2])/3)
           variant.matrix[k,2] = z
           variant.matrix[k,3] = z+1
           m<-reference_Genome [j,3]-reference_Genome [j,2]
          variant.matrix[k,4] = m
         variant.matrix[k,5] = round(m/3)
        }
 
```
------------------------------------------------------------------------------------------------------
### Find number of mutation in each domain for isolates

```R
for (j in 1 :length.domain)
{ 
  for (i in 1:length.varinats)

    if (variant[i,1]== hmmrHit[j,1] & variant[i,2]>=as.numeric( hmmrHit[j,3]) & variant[i,3] <= as.numeric(hmmrHit[j,4]))
      var =append(var,hmmrHit[j,2])
print(j)
}
domain.per.isolate <- table(var)
Domain.Isolate[1,names(domain.per.isolate)] <-  as.numeric(domain.per.isolate) 

```
------------------------------------------------------------------------------------------------------
 ### Merge all domain isolates in one matrix 
 ```R
z=list.files("../Protein_Active-outputs/DomainIsolates/",full.names = T)

myMergedData <-  do.call(rbind,lapply(z,function(x) read.csv(x)))
```
---------------------------------------------------------------------------------------------------
### prermuattion test
 ```R
  x<- t(myMergedData)
  y<-c(rep(1,R),rep(0,S))
  n = 1000
  
  ndist<- matrix(ncol=n, nrow=nrow(x))
  set.seed(1)
  
  np.value<- rep(NA,nrow(x))
  for (i in 1:nrow(x))
  {
    ndist<- replicate(n,diff(by(x[i,],sample(y,length(y),FALSE),mean)))
    print(i)
    np.value[i] <- sum(abs(ndist) > abs(diff(by(x[i,],y, mean))))/n
  }
  
  np.adjusted <- p.adjust(np.value,method="fdr")
  np.adjusted.significant <- which(np.adjusted < 0.0000000005)
  sig.domains<-myMergedData[,np.adjusted]
  ```
 ---------------------------------------------------------------------------------------------------
### Plot the line for average number of mutations for each group (resistanc, suseptible)
```R
plot(hist(ndist))
abline(v=abs(diff(by(x[i,],sample(y,length(y),FALSE), mean))),col="blue")
```
