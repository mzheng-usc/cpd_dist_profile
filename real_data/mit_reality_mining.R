require(gtools)
library(boot)
library(magrittr)
library(MASS)
#library(energy)
library(Hotelling)
library(gTests)
library(matrixStats)
library(ade4)
library(mvtnorm)
require(usedist)
require(reshape2)
library(ggplot2)
library(igraph)
##### new package

require("gSeg")
require(Rcpp)
require("kerSeg") 
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../functions/objTest_fctns.R")
source("../functions/depth_CPD_func.R") 
source("../functions/ecp_distmat_input.R")
source("../functions/kcp_distmat_input.R")
sourceCpp('../functions/energyChangePoint.cpp')


####### first load the dataset and compress it to daily data
####### since original dataset has time interval of 4-hour
####### then calculate the laplacian matrix
MIT<-get(load("reality_mining_1392.RData"))
G_list<-list()
for (i in 1:(dim(MIT)[3]/6)){
  G_list[[i]]<-rowSums(MIT[,,(i*6-5):(i*6)],dims=2)
}
Data<-c()
for (i in 1:(dim(MIT)[3]/6)){
  g<-as.matrix(laplacian_matrix(graph_from_adjacency_matrix(G_list[[i]])))
  Data<-rbind(Data,c(2*g[upper.tri(g)],diag(g)))
}

####### calculate the distance matrix 
distmat<-as.matrix(dist(Data , method = 'manhattan' ))

####### run depth_CPD function, and get the estimated change point location
num_permut<-1
depth_result<-depth_CPD(distmat,num_permut =num_permut,c=0.1)
depth_result$loc

#calculate the scan statistics sequence and plot it 
n=nrow(distmat);c<-0.1
test<-c()
for (cp in seq(ceiling(n*c),n-ceiling(n*c),1)){
  testStat <- getT( distmat = distmat, indices = 1:n, n = cp, m = n-cp, cut_off = 0 )
  test<-cbind(test,testStat)
}


scan_stat<-c(rep(0,ceiling(n*c)-1),test,rep(0,ceiling(n*c)))
date_seq<-seq(as.Date("2004/9/14"), as.Date("2004-09-14")+dim(distmat)[1]-1, "days")

df<-data.frame(scan_stat=scan_stat,day=date_seq)
ggplot(df,aes(x=day,y=scan_stat))+geom_line()+xlab("Day") + ylab("Scan statistic")+
  geom_vline(aes(xintercept = as.numeric(date_seq[depth_result$loc])-1,col='change point location \n 2004-12-15'))+
  scale_color_manual(name = "Dist-CP", values = c("change point location \n 2004-12-15" = "red"))+
  scale_x_date(date_breaks = "22 days", date_labels =  "%d %b %Y")+
  theme_bw()+theme(legend.position=c(0.85, 0.9))+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(panel.grid.minor = element_line(linetype="dashed"))+theme(panel.grid.major = element_line(linetype="dashed"))
ggsave("mit_scan_stat.pdf",width = 7, height = 5,units="in")
 
