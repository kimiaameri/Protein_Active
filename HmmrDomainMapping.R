argv <- commandArgs(trailingOnly = TRUE)
sourcePath <- argv[1]
GenomebedPath <- argv[2]
intersectionspath <- argv[3]
outputPath <- argv[4]
variantPath <- argv[5]
DomainIsolates<-argv[6]
Domain.Isolate<-argv[7]
Domain.Isolate.norm<- argv[8]
#---------------------------------------------------------------------------------------------------#
#                 Please Change the Setwd to the location you download the folder                   #
#---------------------------------------------------------------------------------------------------#
#source(paste0(sourcePath,"geneDomain.R"))
source(paste0(sourcePath,"MutationPosition.R"))
#-----------------------------------------------------------------------#
#                             read hmmrfiles                             #
#-----------------------------------------------------------------------#
hmmrHit <- read.table(paste0(outputPath,"hmmrfile/sort.hit.S.areus.csv"),header=F,sep=" ",stringsAsFactors = F)
length.hmmrHit <-  nrow(hmmrHit)
intersections<- list.files(intersectionspath)
Domain.length<- as.numeric(hmmrHit[,4]) - as.numeric(hmmrHit[,3])
hmmrHit<-cbind(hmmrHit,Domain.length)
colnames(hmmrHit)<- c("Gene.Id","Domain.Name","Start.DomPos","End.DomPos","Domain.Length")
uniq.Domains<- length(unique(hmmrHit[,2]))

#------------------------------------------------------------------------------------#
#                                                                                    #
#                     find the mutation position in each isolate                      #
#                                                                                    #
#------------------------------------------------------------------------------------#
reference_Genome<- read.table(paste0(GenomebedPath,"/nctc8325.bed"),header=F,sep="\t",stringsAsFactors = F)
#MutationPosition(reference_Genome,intersections,intersectionspath,outputPath)

#-----------------------------------------------------------------------#
#            calculate number of mutations per each gene                #
#-----------------------------------------------------------------------#
variant.matrix<- read.table(variantPath,header=T,sep=",",stringsAsFactors = F)

Domain.Isolate <- matrix(0,nrow=1,ncol=uniq.Domains) 
colnames(Domain.Isolate) <- unique(hmmrHit[,2])
Domain.Isolate.norm<-Domain.Isolate
length.varinats= nrow(variant.matrix)
length.domain= nrow(hmmrHit)
var<-NULL 
for (j in 1 :length.domain)
{ 
  for (i in 1:length.varinats)
    if (variant.matrix[i,1]== hmmrHit[j,1] & variant.matrix[i,2]>= hmmrHit[j,3] & hmmrHit[j,4] <=variant.matrix[i,3]) 
      var =append(var,hmmrHit[j,2])
print(j)
}
domain.per.isolate <- table(var)
Domain.Isolate[1,names(domain.per.isolate)] <-  as.numeric(domain.per.isolate) 
Domain.Isolate.norm[1,names(domain.per.isolate)] <-  as.numeric(domain.per.isolate) / as.numeric(hmmrHit[i,5])
write.csv(x =Domain.Isolate,file=DomainIsolates, row.names = FALSE, quote=FALSE )
write.csv(x =Domain.Isolate.norm,file=Domain.Isolate.norm, row.names = FALSE, quote=FALSE )

