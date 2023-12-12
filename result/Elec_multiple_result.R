require(gtools)
require(usedist)
require(reshape2)
library(ggplot2)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../functions/binary_segmentation.R")
MIT_path<-"/Users/mzheng/Desktop/Research/Prof_Dubey/elec_result"
file_list = mixedsort(list.files(path = MIT_path ,full.names=TRUE,recursive=TRUE))
file_list
minLen<-10
seed_I<-seeded.intervals(264,decay = sqrt(2))
seed_I<-seed_I[(seed_I[,2]-seed_I[,1])>minLen,]
nrow(seed_I)


mit_loc<-c();mit_pval<-c();mit_ob<-c();mit_permute<-c()
for (i in file_list){
  d<-get(load(i))
  dim(d$test)
  mit_ob<-c(mit_ob,d$observed_test_statistics)
  mit_permute<-rbind(mit_permute,d$permuted_test_statistics)
  mit_pval<-c(mit_pval,d$p_val)
  mit_loc<-c(mit_loc,d$loc)
}
mit_threshold<-quantile(mit_permute[1,],0.95)

#we break the loop, as soon as there is a point doesn't exceed the threshold
mit_depth_max<-c();IND<-c(1:length(mit_ob));I<-c()

for (i in order(mit_ob,decreasing = TRUE)){
  if (i %in% IND){
    #if (mit_pval[i]>0.01){
    if (mit_ob[i]<=mit_threshold){
      break
    }else{mit_depth_max<-c(mit_depth_max,mit_loc[i])}
  }else{next}
  val<-mit_loc[i]
  ind<-which(seed_I[,1]>val|seed_I[,2]<val,)
  IND<-intersect(IND, ind)
  I<-c(I,i)
  
}
sort(mit_depth_max)
#78, 171, 219

X<-read.csv("../real_data/elec_processed.csv")
geo_dist<-function(x,y){
  acos(sum(sqrt(x)*sqrt(y)))
}

distmat<-as.matrix(dist_make(X,geo_dist))
elec_data<-X
df_elec<-as.data.frame(elec_data)
names(df_elec)<-c("coal","petroleum","gas","nuclear",'conventional hydro',"renewable",'solar')
df_elec$time<-seq(as.Date("2001/1/1"), as.Date("2022/12/1"), "months")
df <- melt(df_elec ,  id.vars = 'time', variable.name = 'Resources')
time<-seq(as.Date("2001/1/1"), as.Date("2022/12/1"), "months")
time[77]
time[sort(mit_depth_max)-1]

ggplot(df, aes(time,value)) + 
  geom_line(aes(colour = Resources,lty=Resources))+
  xlab("Month") + ylab("Composition")+geom_vline(xintercept = c(as.numeric(df$time[170]),as.numeric(df$time[77]),as.numeric(df$time[218])), color = "black",lty="dashed")+
  theme_bw()+
  scale_x_date(date_breaks = "40 month", date_labels =  "%b %Y")+
  theme(panel.grid.minor = element_line(linetype="dashed"))+theme(panel.grid.major = element_line(linetype="dashed"))+
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("elec_multiple_cpd.pdf",width = 7, height = 5,units="in")
citation("GreedySBTM")
