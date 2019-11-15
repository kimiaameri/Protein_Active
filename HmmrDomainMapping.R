argv <- commandArgs(trailingOnly = TRUE)
sourcePath <- argv[1]
GenomebedPath <- argv[2]
intersectionspath <- argv[3]
outputPath <- argv[4]
variantPath <- argv[5]
DomainIsolates<-argv[6]
#---------------------------------------------------------------------------------------------------#
#                 Please Change the Setwd to the location you download the folder                   #
#---------------------------------------------------------------------------------------------------#
#source(paste0(sourcePath,"geneDomain.R"))
#source(paste0(sourcePath,"MutationPosition.R"))
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
MutationPosition(reference_Genome,intersections,intersectionspath,outputPath)
length.genome<- nrow(reference_Genome)
intersections <- read.table(intersectionspath,header=F,sep="\t",stringsAsFactors = F)
   length.intersection= nrow(intersections)
    variant<- matrix(NA, ncol=5, nrow=length.intersection)
    colnames(variant.matrix)<- c("Gene.Id","Variant.start","Variant.end","Gene.length","Chromosome.Length")
      for (k in 1:length.intersection) 
        for (j in 1 :length.genome)
         if (intersections[k,2] >= reference_Genome[j,2] & intersections[k,2] <= reference_Genome[j,3]) 
        { 
           variant[k,1] = as.character(reference_Genome [j,4])
           z<-round(abs(intersections[k,2] - reference_Genome [j,2])/3)
           variant[k,2] = z
           variant[k,3] = z+1
           m<-reference_Genome [j,3]-reference_Genome [j,2]
          variant[k,4] = m
         variant[k,5] = round(m/3)
        }
inputs<- gsub(pattern = ".bed",replacement = ".csv",intersectionspath, perl = T)
variant<- variant[complete.cases(variant),]
write.csv(x=variant,file = paste0(variantPath ,"inputs")), row.names = FALSE, quote=FALSE )

#-----------------------------------------------------------------------#
#            calculate number of mutations per each Domain              #
#-----------------------------------------------------------------------#
variant<- read.table(variantPath,header=T,sep=",",stringsAsFactors = F)

Domain.Isolate <- matrix(0,nrow=1,ncol=uniq.Domains) 
colnames(Domain.Isolate) <- unique(hmmrHit[,2])
length.varinats= nrow(variant)
length.domain= nrow(hmmrHit)
var<-NULL 
for (j in 1 :length.domain)
{ 
  for (i in 1:length.varinats)

    if (variant[i,1]== hmmrHit[j,1] & variant[i,2]>=as.numeric( hmmrHit[j,3]) & variant[i,3] <= as.numeric(hmmrHit[j,4]))
      var =append(var,hmmrHit[j,2])
print(j)
}
domain.per.isolate <- table(var)
Domain.Isolate[1,names(domain.per.isolate)] <-  as.numeric(domain.per.isolate) 
write.csv(x =Domain.Isolate,file=DomainIsolates, row.names = FALSE, quote=FALSE )

