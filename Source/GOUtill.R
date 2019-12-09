#-------------------------------------------------------------------------#
#                           annotation file list                          #
#-------------------------------------------------------------------------
annotation<- read.csv("~/Dropbox/INDEPENDENT sTYDY 2/data/ann.txt",header = FALSE )
annlist<-as.matrix(annotation[,1])
protein.ann.listid<-str_split(annlist, "\t", simplify = TRUE)
#-------------------------------------------------------------------------#
#                           background.txt                                #
#-------------------------------------------------------------------------
background.list<- protein.ann.listid[which((protein.ann.listid[,1]%in% uniprot.Ids) ==TRUE),1]

write.table(background.list,"~/Dropbox/INDEPENDENT sTYDY 2/data/background.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
#-------------------------------------------------------------------------#
#                          target.txt                                     #
#-------------------------------------------------------------------------#
target.gene.list<- uniprot.ids[which((uniprot.ids[,2]%in% list.of.significant.genes) ==TRUE),1]

write.table(target.gene.list,"~/Dropbox/INDEPENDENT sTYDY 2/data/target.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)

