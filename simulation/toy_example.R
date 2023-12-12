setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../functions/objTest_fctns.R")
source("../functions/depth_CPD_func.R") 


n1<-100;n2<-200
Data<-c(rnorm(n1,0,1),rnorm(n2,2,1))


distmat<-as.matrix(dist( Data, method = 'euclidean' ))
dSup <- range(distmat)
dSup <- seq( dSup[1], dSup[2], length.out = 1000 )

#look at the depth profile at Y_200 at u=1/3 and u=1/2
r<-getF(distmat,n=n1,m=n2) #list(FXX,FYY,FXY,FYX)
m1_1=50;m2=300-m1_1
r1<-getF(distmat,n=m1_1,m=m2) 
m1_2=150;m2=300-m1_2
r2<-getF(distmat,n=m1_2,m=m2) 
m1_3=225;m2=300-m1_3
r3<-getF(distmat,n=m1_3,m=m2) 

val<-c(r[[2]][,100],r[[4]][,100],r1[[2]][,200-m1_1],r1[[4]][,200-m1_1],r2[[2]][,200-m1_2],r2[[4]][,200-m1_2],r3[[1]][,200],r3[[3]][,200])
lab<-rep(c("scan point at 1/3 (ground truth)",'scan point at 1/6',"scan point at 1/2","scan point at 3/4"),each=2000)
type<-rep(rep(c("b","a"),each=1000),4)
#type<-c(rep(c(TeX("Formula: $Y_1$"),TeX("Formula: $Y_2$")),each=1000),rep(c(TeX("Formula: $Y_1$"),TeX("Formula: $Y_2$")),each=1000))
range<-rep(dSup,8)
df <- data.frame(range,val,lab,type)
names(df) <- c('range',"val","lab","type")
df$lab_f = factor(df$lab, levels=c('scan point at 1/6','scan point at 1/3 (ground truth)','scan point at 1/2','scan point at 3/4'))
p<-ggplot(df,aes(x=range,y=val, group = type,col=type)) + geom_line()+facet_grid(cols = vars(lab_f))
p+xlab("Range") + ylab("Distance Profile")+
  scale_color_manual(labels = c(TeX(r'(distance profile of $Y_{200}$ w.r.t. $Y_1,\ldots,Y_{\[n{u}\]}$)'),TeX(r'(distance profile of $Y_{200}$ w.r.t. $Y_{\[n{u}\]+1}$,\ldots,Y_{n}$)')), values = c("blue", "red"))+
  theme_bw()+theme(legend.text.align = 0)+theme(panel.grid.minor = element_line(linetype="dashed"))+theme(panel.grid.major = element_line(linetype="dashed"))

# TeX(r'($\[km^2\]$)')
# expression(paste("Depth profile w.r.t ",Y[1],ldots,Y[n*tau]))
# expression(paste("Depth profile w.r.t",Y[1],ldots,Y[mu]))
Day <- (1:30)
Type <-  rep(c("A", "B","C"),10)
Value_1 <- runif(30, min=-1, max=2)
Value_2 <- runif(30, min=-1, max=2)
df <- tibble:: tibble(Day, Type, Value_1, Value_2)
df
pivot_longer(df,-c(Day,Type))

library(ggplot2)
library(latex2exp)
library(tidyr)
ggplot(df,aes(x=range,y=val, group = type,col=type)) + geom_line()+facet_grid(cols = vars(lab))


plot(r[[2]][,100],col='red',type='l')+lines(r[[4]][,100])

m1=100;m2=300-m1
r1<-getF(distmat,n=m1,m=m2) 
plot(r1[[2]][,200-m1],col='red',type='l')+lines(r1[[4]][,200-m1])
sum((r1[[2]][,200-m1]-r1[[4]][,200-m1])^2)


S<-c()

for (m1 in seq(30,170)){
  m2=300-m1
  r1<-getF(distmat,n=m1,m=m2) 
  #plot(r1[[2]][,200-m1],col='red',type='l')+lines(r1[[4]][,200-m1])
  S<-c(S,sum((r1[[2]][,200-m1]-r1[[4]][,200-m1])^2))
}
plot(S)
order(S)

