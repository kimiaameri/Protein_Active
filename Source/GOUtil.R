#---------------------------------------------------------------------
#         find the gene related to significant domains                
#---------------------------------------------------------------------
hmmrhit<- read.csv("../Dropbox/INDEPENDENT sTYDY 2/sort.hit.S.areus.csv", sep = " ", header = FALSE)
genelist<-as.matrix(unique(hmmrhit[which(hmmrhit[,2] %in% rownames(z)==TRUE),1:2]))
write.csv(unique(genelist[,1]),"../Dropbox/INDEPENDENT sTYDY 2/gene.txt", row.names = FALSE, quote = FALSE)
#---------------------------------------------------------------------
listof.mutated.genes<-unique(genelist[,1])
#---------------------------------------------------------------------
#         find the gene ids in uniprot               
#---------------------------------------------------------------------
library(stringr)
uniprot<- read.csv("../Dropbox/INDEPENDENT sTYDY 2/uniprot-taxonomy_93061.tab", header = TRUE, sep="\t")
gene_ids<-as.matrix(uniprot[,"Gene.names"])
#uniIds<-unlist(strsplit(gene_ids[,1], ' '))
pid<-str_split(gene_ids, " S", simplify = TRUE)
uniprot.ids<- gsub("AOUHSC_", "SAOUHSC_", pid)
uniprot.ids<- gsub("SS", "S", uniprot.ids)
#uniprot.ids<- gsub(" ", "/", uniprot.ids)
uniIds<-unlist(strsplit(uniprot.ids[,1],' '))

write.table(uniIds,"../Dropbox/INDEPENDENT sTYDY 2/uniprot_ids.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
#-------------------------------------------------------------------------
output<- read.csv("../Dropbox/INDEPENDENT sTYDY 2/output.txt",header = FALSE )
nlist<-as.matrix(output[,1])
plistid<-str_split(nlist, "\t", simplify = TRUE)
output.list<- plistid[ which((plistid[,1]%in% uniIds) ==TRUE),]

write.table(output.list,"../Dropbox/INDEPENDENT sTYDY 2/background.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)



