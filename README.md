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

## SANVA preprocessing step by step:
* #### In the first step, low-quality reads were trimmed using Trimmomatic 0.38 that trims the reads with length less than 100 and will cut the adapter in reads, if there is any.
* #### BWA-MEM 0.7.17 was used to align sequencing reads against the reference *S. aureus* .
* #### [NCTC 8325](https://www.ncbi.nlm.nih.gov/nuccore/NC_007795.1 ) complete genome. BWA-MEM is introduced to perform local alignment more accurate and faster than other BWA algorithms for high-quality queries. 
* #### Samtools 1.5 can transform the SAM file to a binary version of the same input (view `<-bt>` and `<flagstat>`parameters) and sort alignments by leftmost coordinates, or by read name (`<sort -o>` parameter). 
* #### Picard MarkDuplicate 2.9.0 can remove duplicate reads by locating and tag them in a BAM file, where duplicate reads are defined as originating from a single fragment of DNA.
* #### Freebayes 0.9.10 (with `<-f>` parameter that specifies reference FASTA file) was used to determine the sequence variants present in the samples. Freebayes generates a VCF file that includes all predicted variants (reference + observed allele) and a brief report on their quality scores. 
* #### VCFTools 4.1 (with vcffilter `<-f QUAL>` and `<DP>` parameters for filtering based on quality and depth, respectively) was used to filter the quality scores and depth of each read. Low quality and low depth (reads with quality and depth less than the average of reads) variants were removed. 
* #### We used BCFtools 1.8 to merge all variants in all isolates (view `<-Ob>` parameter was used to subset and filter VCF files by position and filtering expression, index parameter was created index file and merge `<--force>` was used to merge all variants in all samples in one file). 
* #### Variants were identified and analyzed using SpfEff 4.3, which predicts the effects of a variant using an interval forest approach. 
