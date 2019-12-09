#---------------------------------------------------------------------
#         mapping gene ids with gene name  for Staphyloccuce areuse            
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#         find gene names from uniprot taxonomy            
#---------------------------------------------------------------------
library(stringr)
uniprot<- read.csv("~/Dropbox/INDEPENDENT sTYDY 2/data/uniprot-taxonomy_93061.tab", header = TRUE, sep="\t")
gene_ids<-as.matrix(uniprot[,"Gene.names"])
#uniIds<-unlist(strsplit(gene_ids[,1], ' '))
p.id<-str_split(gene_ids, " S", simplify = TRUE)
uniprot.ids<- gsub("AOUHSC_", "SAOUHSC_", pid)
uniprot.ids<- gsub("SS", "S", uniprot.ids)
#uniprot.ids<- gsub(" ", "/", uniprot.ids)
uniprot.Ids<-unlist(strsplit(uniprot.ids[,1],' '))

write.csv(uniprot.ids,"~/Dropbox/INDEPENDENT sTYDY 2/data/uniprot.Ids.csv", row.names = TRUE, quote = FALSE)
write.table(uniprot.Ids,"~/Dropbox/INDEPENDENT sTYDY 2/data/Uniprot_ids.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
