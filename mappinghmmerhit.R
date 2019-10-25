argv <- commandArgs(trailingOnly = TRUE)
sourcePath <- argv[1]
GenomebedPath <- argv[2]
intersectionspath <- argv[3]
outputPath <- argv[4]
bigtableFile <- argv[5]
bigtableWeightFile<-argv[6]
#---------------------------------------------------------------------------------------------------#
#                 Please Change the Setwd to the location you download the folder                   #
#---------------------------------------------------------------------------------------------------#
#source("geneDomain.R")
source(paste0(sourcePath,"MutationPosition.R"))
#-----------------------------------------------------------------------#
#                             read hmmrfiles                             #
#-----------------------------------------------------------------------#
hmmrHit <- as.matrix(read.table(paste0(outputPath,"/hmmrfile/sort.hit.S.areus.csv",header=F,sep=" ",stringsAsFactors = F)))
length.hmmrHit= nrow(hmmrHit)
intersections<- list.files(paste0(intersectionspath,"/intersections/"))
Domain.length<- as.numeric(hmmrHit[,4]) - as.numeric(hmmrHit[,3])
hmmrHit<-cbind(hmmrHit,Domain.length)
colnames(hmmrHit)<- c("Gene.Id","Domain.Name","Start.DomPos","End.DomPos","Domain.Length")
uniq.Domains<- length(unique(hmmrHit[,2]))

#------------------------------------------------------------------------------------#
#                                                                                    #
#                     find the mutation position in each isolate                      #
#                                                                                    #
#------------------------------------------------------------------------------------#
reference_Genome<- as.matrix(read.table(paste0(GenomebedPath,"/nctc8325.bed",header=F,sep="\t",stringsAsFactors = F)))
MutationPosition(reference_Genome,intersections,outputPath)
