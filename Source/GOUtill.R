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

#-------------------------------------------------------------------------#
#               annotation file list for domains                          #
#-------------------------------------------------------------------------
mappdb2Pfam<- read.csv("~/Dropbox/INDEPENDENT sTYDY 2/pfam/hmmer_pdb_all.csv", header = TRUE , sep = ",")
mapGO2Pfam<- read.csv("~/Dropbox/INDEPENDENT sTYDY 2/pfam/pfam_gomf_most_specific.txt", header = TRUE, sep=";")

pfams<-str_split(as.character(mappdb2Pfam$PFAM_ACC), "\\.", simplify = TRUE)
mappdb2Pfam$PFAM_ACC<-pfams[,1]
colnames(mappdb2Pfam)[5] <-"PFAM"
mapGo2PDB <- mapGO2Pfam[which(mapGO2Pfam$PFAM %in% mappdb2Pfam$PFAM),]


combined <- sort(union(levels(mapGO2Pfam$PFAM), levels(mappdb2Pfam$PFAM)))
mapGo2PDB<- inner_join(mutate(mapGO2Pfam, PFAM=factor(PFAM, levels=combined)),
                       mutate(mappdb2Pfam, PFAM=factor(PFAM, levels=combined)))

mapGo2Pname<- unique(mapGo2PDB[,c("GO" ,"PFAM_Name")])
pfamName<- hmmrhit[,2]
goPfam<-NA
mapGo2Pname<- mapGo2Pname[which(mapGo2Pname[,2] %in%pfamName ),]
ann<- aggregate(mapGo2Pname[,1], by=list(mapGo2Pname[,2]), toString )
write.table(ann,"~/Dropbox/INDEPENDENT sTYDY 2/domains/ann.domain.txt",sep="\t",col.names=F,row.names=F)  
#-------------------------------------------------------------------------#
#                           background for domains.txt                   #
#-------------------------------------------------------------------------

write.table(hmmrhit[,2] ,"~/Dropbox/INDEPENDENT sTYDY 2/domains/background.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
#-------------------------------------------------------------------------#
#                          target.txt                                     #
#-------------------------------------------------------------------------#

write.table(significan.domians,"~/Dropbox/INDEPENDENT sTYDY 2/domains/target.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)

