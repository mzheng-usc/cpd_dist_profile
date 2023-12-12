#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../../functions/binary_segmentation.R") 
source("../../functions/objTest_fctns.R")
source("../../functions/depth_CPD_func.R") 
require(usedist)
geo_dist<-function(x,y){
  acos(sum(sqrt(x)*sqrt(y)))
}

elec_data<-read.csv("../elec_processed.csv")
D<-as.matrix(dist_make(elec_data,geo_dist))


minLen<-10
seed_I<-seeded.intervals(nrow(D),decay = sqrt(2))
seed_I<-seed_I[(seed_I[,2]-seed_I[,1])>minLen,]
nrow(seed_I)

i<-as.numeric(commandArgs(TRUE)[1])
num_permut<-1000
distmat_seed<-D[seed_I[i,1]:seed_I[i,2],seed_I[i,1]:seed_I[i,2]]
r<-depth_CPD_permutated(distmat_seed,num_permut,c=0.3)
r$loc=r$loc+seed_I[i,1]-1

path<-paste("Elec_iter_",i,'.Rdata',sep="")
save(r, file=path)