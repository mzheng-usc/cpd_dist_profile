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


##### define the geodesic distance for compositional data 
geo_dist<-function(x,y){
  acos(sum(sqrt(x)*sqrt(y)))
}

##### load the preprocessed dataset
X<-read.csv("elec_processed.csv")
distmat<-as.matrix(dist_make(X,geo_dist))
elec_data<-X
##### run depth_CPD function, to get the estimated change point location
num_permut<-1
depth_result<-depth_CPD(distmat,num_permut =num_permut,c=0.1)
depth_result$loc

##### we combine the result and formulate into a long dataframe form
df_elec<-as.data.frame(elec_data)
names(df_elec)<-c("coal","petroleum","gas","nuclear",'conventional hydro',"renewable",'solar')
df_elec$time<-seq(as.Date("2001/1/1"), as.Date("2022/12/1"), "months")
df <- melt(df_elec ,  id.vars = 'time', variable.name = 'Resources')

##### plot the preprocess data with each composition of each resource,
##### along with a red vertical line indicating the estiamted change point location
ggplot(df, aes(time,value)) + 
  geom_line(aes(colour = Resources,lty=Resources))+
  xlab("Month") + ylab("Composition")+geom_vline(xintercept = as.numeric(df$time[170]), color = "black",lty="dashed")+
  theme_bw()+
  scale_x_date(date_breaks = "40 month", date_labels =  "%b %Y")+
  theme(panel.grid.minor = element_line(linetype="dashed"))+theme(panel.grid.major = element_line(linetype="dashed"))+
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave("elec_cpd.pdf",width = 7, height = 5,units="in")


##### calculate the scan statistics on interval I_c and plot it
n=nrow(distmat);c<-0.1
test<-c()
for (cp in seq(ceiling(n*c),n-ceiling(n*c),1)){
  testStat <- getT( distmat = distmat, indices = 1:n, n = cp, m = n-cp, cut_off = 0 )
  test<-cbind(test,testStat)
}

scan_stat<-c(rep(0,ceiling(n*c)-1),test,rep(0,ceiling(n*c)))
date_seq<-seq(as.Date("2001/1/1"), as.Date("2022/12/1"), "months")
df<-data.frame(scan_stat=scan_stat,day=date_seq)
ggplot(df,aes(x=day,y=scan_stat))+geom_line()+xlab("Month") + ylab("Scan statistic")+
  geom_vline(aes(xintercept = as.numeric(date_seq[depth_result$loc-1]),col='change point location \n Feb 2015'))+
  scale_color_manual(name = "Dist-CP", values = c("change point location \n Feb 2015" = "red"))+
  scale_x_date(date_breaks = "22 month", date_labels =  "%b %Y")+
  theme_bw()+theme(legend.position=c(0.85, 0.9))+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(panel.grid.minor = element_line(linetype="dashed"))+theme(panel.grid.major = element_line(linetype="dashed"))
ggsave("elec_scan_stat.pdf",width = 7, height = 5,units="in")


# 
# Date<-c()
# date_seq<-seq(as.Date("2001/1/1"), as.Date("2022/12/1"), "months")
# for (x in 1:length(date_seq)){
#   Date<-c(Date,format(date_seq[x], "%Y-%m"))
# }
# 
# library(lubridate)
# (d <- ymd("2001-01-01"))
# d %m+% months(170)
# 
