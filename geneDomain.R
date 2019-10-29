argv <- commandArgs(trailingOnly = TRUE)

geneDomain.vcf <- function(i, hmmrhit,outputPath)
{
  varinats <- read.table(paste0(outputPath,varinats[i],sep=""),header=F,sep="\t",stringsAsFactors = F)
  length.varinats= nrow(varinats)
  length.domain= nrow(hmmrhit)
  for (j in 1 :length.domain)
  { 
    for (i in 1:length.varinats)
      if (varinats[i,1]]= hmmrhit[j,1] & varinats[i,2]]= hmmrhit[j,3] & varinats[i,3]]= hmmrhit[j,4]) 
        varinats[i,4] = reference [j,4]
  }
  
  return(table(varinats[,4]))
}
