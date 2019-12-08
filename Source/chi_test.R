memory.limit(90000000)
memory.size(90000000)
#-----------------------------------------------------------------------------------------#
#                 merge all mutations in the domian for each isolate                      #
#-----------------------------------------------------------------------------------------#
z=list.files("~/Dropbox/INDEPENDENT sTYDY 2/data/DomainIsolates/",full.names = T)
inside.Domain.Data <-  do.call(rbind,lapply(z,function(x) read.csv(x)))
rname <- gsub(".*\\/","",z)
rownames(inside.Domain.Data) <- gsub(".csv","",rname)
inside.Domain.Data<- inside.Domain.Data[,colSums(inside.Domain.Data)!=0]
write.csv(inside.Domain.Data,"~/Dropbox/INDEPENDENT sTYDY 2/data/inside.domian.csv" , quote=FALSE)
#-----------------------------------------------------------------------------------------#
#                        find resistance type for each isolate                            #
#-----------------------------------------------------------------------------------------#
sample<-read.csv("~/Dropbox/INDEPENDENT sTYDY 2/data/sample.csv", header = TRUE, stringsAsFactors = T)
in.domain<- cbind(inside.Domain.Data, sample[,2])
colnames(in.domain) <- c(colnames(in.domain[-ncol(in.domain)]),"lable")
write.csv(in.domain,"~/Dropbox/INDEPENDENT sTYDY 2/data/inside.domian.labeled.csv" , quote=FALSE)
#in.domain<- read.csv("~/Dropbox/INDEPENDENT sTYDY 2/data/mergewithlable.csv", header = TRUE)
#listdata<- in.domain[,1]
#rownames(in.domain)<- listdata
#in.domain<- in.domain[,-1]
#-----------------------------------------------------------------------------------------#
#             merge all mutations outside the domian for each isolate                     #
#-----------------------------------------------------------------------------------------#
z1=list.files("~/Dropbox/INDEPENDENT sTYDY 2/NotDomainIsolates/",full.names = T)
outside.Domain.Data <-  do.call(rbind,lapply(z1,function(x) read.csv(x)))
rname1 <- gsub(".*\\/","",z1)
rownames(outside.Domain.Data) <- gsub(".csv","",rname1)
outside.Domain.Data<- outside.Domain.Data[,colSums(outside.Domain.Data)!=0]
out.domain<-cbind(outside.Domain.Data, sample[,2])
colnames(out.domain) <- c(colnames(out.domain[-ncol(out.domain)]),"lable")

#-----------------------------------------------------------------------------------------#
#                    find antimicrobial type for each isolate                             #
#-----------------------------------------------------------------------------------------#
inside.suseptible.group<- in.domain[which(in.domain$lable == "Susceptible-gentamicin"),]
inside.Resistant.group<- in.domain[which(in.domain$lable == "Resistant-gentamicin"),]
s.in.group<- inside.suseptible.group[,-ncol(inside.suseptible.group)]
r.in.group<- inside.Resistant.group[,-ncol(inside.Resistant.group)]
#-----------------------------------------------------------------------------------------#
outside.suseptible.group<- out.domain[which(out.domain$lable == "Susceptible-gentamicin"),]
outside.Resistant.group<- out.domain[which(out.domain$lable == "Resistant-gentamicin"),]
s.out.group<- outside.suseptible.group[,-ncol(outside.suseptible.group)]
r.out.group<- outside.Resistant.group[,-ncol(outside.Resistant.group)]

#-----------------------------------------------------------------------------------------#
#                                         chi-squre test                                  #
#-----------------------------------------------------------------------------------------#
in.r.x<- colSums(r.in.group)
in.s.x<- colSums(s.in.group)
out.r.x<- colSums(r.out.group)
out.s.x<- colSums(s.out.group)
chi.pvalue<-NULL
chi.dif<-NULL
for (i in 1 : ncol(s.in.group))
{
  x<- matrix(c(in.r.x[i], in.s.x[i], sum(in.r.x[-i],out.r.x), sum(in.s.x[-i],out.s.x)), ncol = 2)
  chitest<- chisq.test(x, y = NULL, correct = TRUE,
                       p = rep(1/length(x), length(x)), rescale.p = FALSE,
                       simulate.p.value = FALSE, B = 2000)
  chi.pvalue<- append(chi.pvalue,chitest$p.value)
  dif<- (x[1,1]/x[2,1])/(x[1, 2]/x[2, 2])
  chi.dif<- append(chi.dif,dif)
  print(i)
}
p.adj.chitest<- p.adjust(chi.pvalue,method = "bonferroni")
#domians.chitest<- rbind(domains, chi.dif)

length(which(p.adj.chitest<=0.000000000000000000005))
chi.sig.index<- which(p.adj.chitest<0.000000000000000000005)
significat.domians.chitest<- in.domain[,chi.sig.index]

write.csv(significat.domians.chitest, "~/Dropbox/INDEPENDENT sTYDY 2/sig.chitest.csv")
