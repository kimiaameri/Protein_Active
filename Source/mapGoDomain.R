go.output<- read.csv("~/Dropbox/INDEPENDENT sTYDY 2/enrich-result/target5-output.txt", sep= "\t", header = FALSE)
enrich.genes<- as.matrix(go.output[,ncol(go.output)])
Go.genes<-str_split(enrich.genes, " ", simplify = TRUE)

go.list.gene<-unlist(strsplit(enrich.genes[,1],' '))
length(unique(go.list.gene))
Go.list.gene<-unique(go.list.gene)
#---------------------------------------------------------------------
#         find the gene locuse related to genes from GO enrichmet analysis                
#---------------------------------------------------------------------
Go.gene.id<- uniprot.ids[which (Go.list.gene %in% uniprot.ids[,1]),2]

#---------------------------------------------------------------------
#         find the Domiain related to each Go enriched gene                
#---------------------------------------------------------------------
Go.domain<- hmmrhit[which(Go.gene.id %in% hmmrhit[,1]),2]
Go.Domain.list<- unique(Go.domain)
length(Go.Domain.list)
