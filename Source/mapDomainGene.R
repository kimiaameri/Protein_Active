#---------------------------------------------------------------------
#         mapping domians with gene ids for Staphyloccuce areuse             
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#         find the gene related to significant domains                
#---------------------------------------------------------------------
hmmrhit<- read.csv("~/Dropbox/INDEPENDENT sTYDY 2/sort.hit.S.areus.csv", sep = " ", header = FALSE)
genelist<-as.matrix(unique(hmmrhit[which(hmmrhit[,2] %in% significan.domians==TRUE),1:2]))

mutated.gene.list<-as.matrix(unique(genelist[which(genelist[,1] %in% mutated.genes==TRUE),1:2]))

write.table(unique(mutated.gene.list[,1]),"~/Dropbox/INDEPENDENT sTYDY 2/data/genelist.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
#---------------------------------------------------------------------
list.of.significant.genes<-unique(mutated.gene.list[,1])
