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
This bash code can aslo run as:
```bash
sbatch SRA.sh  
```
####
 
 ## SANVA Pipeline step by step:
