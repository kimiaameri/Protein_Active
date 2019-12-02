memory.limit(90000000)
memory.size(90000000)
#load("~/Dropbox/snp.matrix.bin")
l.domain<- read.csv("~/Dropbox/INDEPENDENT sTYDY 2/data/mergewithlable.csv", header = TRUE)
listdata<- l.domain[,1]
rownames(l.domain)<- listdata
l.domain<- l.domain[,-1]

#---------------------------------------------------
#           keep non zero domains for each group
#---------------------------------------------------
just.r.group<- r.group[,which(colSums(r.group)!=0)]
just.s.group<- s.group[,which(colSums(s.group)!=0)]
uniqu.domains.resistance<- setdiff(colnames(just.r.group), colnames(just.s.group))
uniqu.domains.suseptible<- setdiff(colnames(just.s.group), colnames(just.r.group))
#---------------------------------------------------
#           keep  domains with more than 1.qu each group
#---------------------------------------------------
colsums.domian<- colSums(l.domain[,-ncol(l.domain)])
summary(colsums.domian)

domains<-cbind( l.domain[,which(colSums(l.domain[,-ncol(l.domain)])>200)], l.domain[,ncol(l.domain)])

colnames(domains) <- c(colnames(domains[-ncol(domains)]),"lable")
suseptible.group<- domains[which(domains$lable == "Susceptible-gentamicin"),]
Resistant.group<- domains[which(domains$lable == "Resistant-gentamicin"),]
s.group<- suseptible.group[,-ncol(suseptible.group)]
r.group<- Resistant.group[,-ncol(Resistant.group)]
#--------------------------------------------#
#------------wilcox.test
#--------------------------------------------#
p.value.wilcox<- NULL
for (i in 1: ncol(s.group))
{
  res<- wilcox.test(s.group[,i], r.group[,i])
  p.value.wilcox<- append(p.value.wilcox,res$p.value)
  print(i)
}
p.adj.wilcox<- p.adjust(p.value.wilcox,method = "bonferroni")
length(which(p.adj.wilcox<0.0000000005))
sig.index.wilcox<- which(p.adj.wilcox<0.0000000005)
significat.domians.wilcox<- domains[,sig.index.wilcox]
write.csv(significat.domians.wilcox, "./Outputs/sig.wilcox.csv")

#--------------------------------------------#
#---------- chi-squre test
#--------------------------------------------#
r.x<- colSums(r.group)
s.x<- colSums(s.group)
chi.pvalue<-NULL
for (i in 1 : ncol(s.group))
{
  x<- matrix(c(r.x[i], s.x[i], sum(r.x[-i]), sum(s.x[-i])), ncol = 2)
  chitest<- chisq.test(x, y = NULL, correct = TRUE,
             p = rep(1/length(x), length(x)), rescale.p = FALSE,
             simulate.p.value = FALSE, B = 2000)
  chi.pvalue<- append(chi.pvalue,chitest$p.value)
  print(i)
}
p.adj.chitest<- p.adjust(chi.pvalue,method = "bonferroni")
length(which(p.adj.chitest<=0.000000000000000000005))
chi.sig.index<- which(p.adj.chitest<0.000000000000000000005)
significat.domians.chitest<- domains[,chi.sig.index]
write.csv(significat.domians.chitest, "./Outputs/sig.chitest.csv")
