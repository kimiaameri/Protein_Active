install.packages("devtools")

##BiocManager::install("genbankr")

install.packages("backports")
library(devtools)
install_git("https://git.bioconductor.org/packages/genbankr")
library(genbankr)

##source('https://bioconductor.org/biocLite.R')
##BiocManager::install("BiocUpgrade")
##BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
library("Biostrings")
library(bio3d)
#-----------------------------read protein translations ----------------------------------#
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
write.table(x=translations,file="./Outputs/translations.txt",quote =FALSE,sep="\n",eol = "\n", row.names = FALSE, col.names = FALSE)

