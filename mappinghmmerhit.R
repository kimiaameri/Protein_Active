#---------------------------------------------------------------------------------------------------#
#                 Please Change the Setwd to the location you download the folder                   #
#---------------------------------------------------------------------------------------------------#
#source("geneDomain.R")
source("MutationPosition.R")
#-----------------------------------------------------------------------#
#                             read hmmrfiles                             #
#-----------------------------------------------------------------------#
hmmrHit <- as.matrix(read.table("/hmmrfile/sort.hit.S.areus.csv",header=F,sep=" ",stringsAsFactors = F))
length.hmmrHit= nrow(hmmrHit)
intersections<- list.files("/intersection/")
Domain.length<- as.numeric(hmmrHit[,4]) - as.numeric(hmmrHit[,3])
hmmrHit<-cbind(hmmrHit,Domain.length)
colnames(hmmrHit)<- c("Gene.Id","Domain.Name","Start.DomPos","End.DomPos","Domain.Length")
uniq.Domains<- length(unique(hmmrHit[,2]))

#------------------------------------------------------------------------------------#
#                                                                                    #
#                     find the mutation position in each isolate                      #
#                                                                                    #
#------------------------------------------------------------------------------------#
inputs<- list.files("/intersections/")
reference_Genome<- as.matrix(read.table("/nctc8325.bed",header=F,sep="\t",stringsAsFactors = F))
MutationPosition(reference_Genome,inputs)
