require(igraph)
require(ecp)
require(Rcpp)
require("kerSeg") 
source("../../functions/binary_segmentation.R") 
source("../../functions/objTest_fctns.R")
source("../../functions/depth_CPD_func.R") 
source("../../functions/ecp_distmat_input.R")
source("../../functions/kcp_distmat_input.R")
sourceCpp('../../functions/energyChangePoint.cpp')


iter<-as.numeric(commandArgs(TRUE)[1])

n_nodes<-300

n1=n2=n3=n4=100
num_permut<-1000

m1<-c()
for (i in 1:(n1)){
  pm <- cbind(c(0.2,0.001,0.001),c(0.001,0.2,0.001),c(0.001,0.001,0.2))
  g1 <- sample_sbm(n_nodes, pref.matrix = pm, block.sizes = c(100,100,100))
  m<-c(as.matrix(as.matrix(g1, "adjacency")))
  m1<-rbind(m1,m)
}
for (i in 1:n2){
  pm <- cbind(c(0.8,0.001,0.001),c(0.001,0.2,0.001),c(0.001,0.001,0.8))
  g2 <- sample_sbm(n_nodes, pref.matrix = pm, block.sizes = c(100,100,100))
  m<-c(as.matrix(as.matrix(g2, "adjacency")))
  m1<-rbind(m1,m)
}

for (i in 1:n3){
  pm <- cbind(c(0.8,0.001,0.001),c(0.001,0.2,0.001),c(0.001,0.001,0.8))
  g3 <- sample_sbm(n_nodes, pref.matrix = pm, block.sizes = c(200,50,50))
  m<-c(as.matrix(as.matrix(g3, "adjacency")))
  m1<-rbind(m1,m)
}
for (i in 1:n4){
  pm <- cbind(c(0.5,0.01),c(0.01,0.5))
  g4 <- sample_sbm(n_nodes, pref.matrix = pm, block.sizes = c(200,100))
  m<-c(as.matrix(as.matrix(g4, "adjacency")))
  m1<-rbind(m1,m)
}

distmat<-as.matrix(dist( m1, method = 'manhattan' ))
minLen<-20
seed_I<-seeded.intervals(nrow(distmat),decay = sqrt(2))
seed_I<-seed_I[(seed_I[,2]-seed_I[,1])>minLen,]
nrow(seed_I)

c<-0.1
result<-c()
for (i in 1:nrow(seed_I)){
  distmat_seed<-distmat[seed_I[i,1]:seed_I[i,2],seed_I[i,1]:seed_I[i,2]]
  r<-depth_CPD_permutated(distmat_seed,num_permut,c=0.1)
  r$loc=r$loc+seed_I[i,1]-1
  result<-rbind(result,r)
}

colnames(result) <- c("p_val",'loc','observed_test_statistics','permuted_test_statistics')

D<-distmat
#ecp
result_ecp<-e.divisive_distmat(D=distmat,sig.lvl=.05,R=999,k=NULL,min.size=30,alpha=1)

#kcp
result_kcp<-kcpa_distmat(D=distmat,L=nrow(D)/sqrt(log(nrow(D))))

Result<-list()
Result[["depth"]]<-result
Result[["ecp"]]<-result_ecp
Result[["kcp"]]<-result_kcp

path<-paste("SBM_iter_",iter,'.Rdata',sep="")
save(Result, file=path)