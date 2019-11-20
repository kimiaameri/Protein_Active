permutationTest<- function(dataset,R,S)
{
  ##------------------ find significant domains for each isolate -------##
  memory.size(9000000)
  memory.limit(9000000)
  x<- t(dataset)
  y<-c(rep(1,R),rep(0,S))
  n = 1000
  
  ndist<- matrix(ncol=n, nrow=nrow(x))
  set.seed(1)
  
  np.value<- rep(NA,nrow(x))
  for (i in 1:nrow(x))
  {
    ndist<- replicate(n,diff(by(x[i,],sample(y,length(y),FALSE),mean)))
    print(i)
    np.value[i] <- sum(abs(ndist) > abs(diff(by(x[i,],y, mean))))/n
  }
  
  np.adjusted <- p.adjust(np.value,method="fdr")
  np.adjusted.significant <- which(np.adjusted < 0.0005)
  sig.domains<-dataset[,np.adjusted]
}
