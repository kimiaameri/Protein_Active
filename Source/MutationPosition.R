argv <- commandArgs(trailingOnly = TRUE)
MutationPosition<- function(reference_Genome,inputs,intersectionspath,outputPath)
{
  length.genome<- nrow(reference_Genome)
  for (i in 1:length(inputs))
  {
    intersections <- read.table(paste(intersectionspath,inputs[i],sep=""),header=F,sep="\t",stringsAsFactors = F)
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
      
    print(paste("i=",i))
    inputs[i]<- gsub(pattern = ".bed",replacement = "",inputs[i], perl = T)
    variant.matrix<- variant.matrix[complete.cases(variant.matrix),]
    write.csv(x=variant.matrix,file = paste0(variantPath,paste0(inputs[i],".csv")), row.names = FALSE, quote=FALSE )
  }
}
