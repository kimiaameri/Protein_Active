argv <- commandArgs(trailingOnly = TRUE)
sourcePath <- argv[1]
GenomebedPath <- argv[2]
intersectionspath <- argv[3]
outputPath <- argv[4]
variantPath <- argv[5]
Domain.Isolate<-argv[6]
Domain.Isolate.norm<- argv[7]
#---------------------------------------------------------------------------------------------------#
#                 Please Change the Setwd to the location you download the folder                   #
#---------------------------------------------------------------------------------------------------#
source(paste0(sourcePath,"geneDomain.R"))
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
#------------------------------------------------------------------------------------#
#                                                                                    #
#                     find the mutaion domain   in each isolate                      #
#                                                                                    #
#------------------------------------------------------------------------------------#
variants<- list.files(variantPath)
Domain.Isolate <- matrix(0,nrow=length(variants ),ncol=uniq.Domains)
isolates<-gsub(pattern = ".csv",replacement = "",variants, perl = T)
rownames(Domain.Isolate) <- isolates
colnames(Domain.Isolate) <- unique(hmmrHit[,2])
Domain.Isolate.norm<-Domain.Isolate
variants<- list.files(variantPath)
for (i in 1:length(variants))
{
  print(i)
  domain.per.isolate <- geneDomain(variants[i],hmmrHit, variantPath)
  Domain.Isolate[i,names(domain.per.isolate)] <-  as.numeric(domain.per.isolate) 
  Domain.Isolate.norm[i,names(domain.per.isolate)] <-  as.numeric(domain.per.isolate) / as.numeric(hmmrHit[i,5])
}
write.csv(x =Domain.Isolate,file=paste0(outputPath,"Domain.Isolate.csv"), row.names = FALSE, quote=FALSE )
write.csv(x =Domain.Isolate.norm,file=paste0(outputPath,"Domain.Isolate.norm.csv")), row.names = FALSE, quote=FALSE )


