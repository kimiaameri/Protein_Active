argv <- commandArgs(trailingOnly = TRUE)

geneDomain.vcf <- function(i, reference,variantPath)
{
    varinat <- read.table(paste0(variantPath,i,sep=""),header=T,sep=",",stringsAsFactors = F)
    length.varinats= nrow(variant)
    length.domain= nrow(hmmrHit)
    var<-NULL 
    for (j in 1 :length.domain)
    { 
      for (i in 1:length.varinats)
        if (variant[i,1]== hmmrHit[j,1] & variant[i,2]>= hmmrHit[j,3] & hmmrHit[j,4] <=variant[i,3]) 
     var =append(var,variant[i,2])
    }
    
    return(table(var))
  }
