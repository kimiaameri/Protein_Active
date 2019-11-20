argv <- commandArgs(trailingOnly = TRUE)
sourcePath <- argv[1]
outputPath<- argv[2]
DomainIsolates <- argv[3]


source(paste0(sourcePath,"/permutationTest.R"))


memory.size(9000000)
memory.limit(9000000)
z=list.files(DomainIsolates,full.names = T)

myMergedData <-  do.call(rbind,lapply(z,function(x) read.csv(x)))

rname <- gsub(".*\\/","",z)
rownames(myMergedData) <- gsub(".csv","",rname)

myMergedData<- myMergedData[,colSums(myMergedData)!=0]
write.csv(myMergedData,paste0(outputPath,"merge.csv" ), quote=FALSE)

R<- 174
S<- 189

sig.domains<- permutationTest(myMergedData, R,S)
data<- t(sig.domains)
sig.domains.act<- myMergedData[,np.adjusted1]
write.csv(data, paste0(outputPath,"significantDomain.csv", row.names = T, quote=FALSE)
