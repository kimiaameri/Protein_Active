#-------------------------------------------------------------------------#
#                          Active sites                                   #
#-------------------------------------------------------------------------#
uniprot<- read.csv("~/Dropbox/INDEPENDENT sTYDY 2/uniprot-yourlist.tab", sep="\t")

actives<-as.matrix(uniprot[,c("Gene.names","Active.site")])
actives.list<-str_split(actives[,2], "\\}.;", simplify = TRUE)
activesite_list<-gsub(".*ACT_SITE ([0-9]+).*","\\1",actives.list)
#activesite_list<-gsub(".*ACT_SITE ([0-9]+ [0-9]+).*","\\1",actives.list)
gene.locuse <- read.csv("~/Dropbox/INDEPENDENT sTYDY 2/data/gene_locus.csv", header = T, sep=",")

Activesite<-apply(activesite_list,2,as.numeric)
rownames(activesite)<-gene.locuse[,1]
Activesite[is.na(Activesite)]<-0
Activesite<- Activesite[which(rowSums(Activesite)!=0),]
#-------------------------------------------------------------------------#
#                          Binding sites                                  #
#-------------------------------------------------------------------------#

Binding<-as.matrix(uniprot[,c("Gene.names","Binding.site")])
binding.list<-str_split(Binding[,2], "\\}.;", simplify = TRUE)
bindingsite_list<-gsub(".*BINDING ([0-9]+).*","\\1",binding.list)
#activesite_list<-gsub(".*ACT_SITE ([0-9]+ [0-9]+).*","\\1",actives.list)

Bindingsite<-apply(bindingsite_list,2,as.numeric)
rownames(Bindingsite)<-gene.locuse[,1]
Bindingsite[is.na(Bindingsite)]<-0
Bindingsite<- Bindingsite[which(rowSums(Bindingsite)!=0),]
#-------------------------------------------------------------------------#
#             if mutaion  happened in Binding site                         #
#-------------------------------------------------------------------------#
variantlist<- list.files("~/Dropbox/INDEPENDENT sTYDY 2/data/VariantPosition/")
for (m in 1:length(variantlist))
{
  variants=read.csv(paste0("~/Dropbox/INDEPENDENT sTYDY 2/data/VariantPosition/",variantlist[m]), header = T, sep = ",")
  matchBindings<-NULL
  for (i in 1:nrow(variants))
  {
    for (j in 1:nrow(Bindingsite))
    {
      if (rownames(Bindingsite)[j]==variants[i,1])
      {
        for (k in 1:ncol(Bindingsite))
          if (Bindingsite[j,k] == variants[i,2])
          {
            matchBindings<- rbind(matchBindings, c(rownames(Bindingsite)[j],variants[i,2]))
            break(k)
          }
      }
    }
    print(i)
  }
  write.csv(matchBindings , file = paste0("~/Dropbox/INDEPENDENT sTYDY 2/BindingSite/",variantlist[m]), row.names = F, quote = F)
  #-------------------------------------------------------------------------#
  #             if mutaion  happened in active site                         #
  #-------------------------------------------------------------------------#
  matchActives<-NULL
  for (i in 1:nrow(variants))
  {
    for (j in 1:nrow(Activesite))
    {
      if (rownames(Activesite)[j]==variants[i,1])
      {
        for (k in 1:ncol(Activesite))
          if (Activesite[j,k] == variants[i,2])
          {
            matchActives<- rbind(matchActives, c(rownames(Activesite)[j],variants[i,2]))
            break(k)
          }
      }
    }
    print(i)
  }  
  write.csv(matchActives , file = paste0("~/Dropbox/INDEPENDENT sTYDY 2/ActiveSite/", variantlist[m]), row.names = F, quote = F)
  
}
