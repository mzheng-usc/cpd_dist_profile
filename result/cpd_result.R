require(gtools)
require(dplyr)
library(ggplot2)
require(reshape2)
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
##### ##### only part need to change, substitute it with your local simulation output path
path_result<-"../simulation_output"
path_result<-"/Users/mzheng/Desktop/Research/Prof_Dubey/cpd_code/MMD_sim"
##### ##### 

result_list<-mixedsort(list.dirs(path = path_result, full.names = TRUE,recursive = F))
result_list


for (path in result_list[grep("toy_example",result_list,invert = TRUE)]){
  file_list = mixedsort(list.files(path = path ,full.names=TRUE,recursive=TRUE))
 
  
  if (grepl("normal_mean",path)){
    set<-1
    save_name<-"normal_mean"
  } else if(grepl("normal_scale",path)) {
    set<-2
    save_name<-"normal_scale"
  }else if(grepl("normal_mix",path)) {
    set<-3
    save_name<-"normal_mix"
  }else if(grepl("normal_tail",path)) {
    set<-4
    save_name<-"normal_tail"
  }else if(grepl("dist_mean",path)) {
    set<-5
    save_name<-"dist_mean"
  }else if(grepl("dist_scale",path)) {
    set<-6
    save_name<-"dist_scale"
  }else if(grepl("network_sim",path)) {
    set<-7
    save_name<-"network_sim"
  }
  
  if (set %in% c(1,3,5)){
    type_seq<-seq(0,1,0.1)
  } else if (set %in% c(2,6)){
    type_seq<-seq(0,0.4,0.04)
  } else if (set %in% c(4)){
    type_seq<-seq(2,22,2)
  } else if (set %in% c(7)){
    type_seq<-seq(0.0000000001,0.5000000001,0.05)
  }
  
  
  result_path<-list()
  for (i in seq_along(1:length(type_seq))){
    value<-c();value_index<-c()
    name<-paste("delta_",type_seq[i],'_',sep="")
    for (file in file_list){
      #cat(file)
      if (grepl(name,file,fixed = TRUE)){
        value<-c(value,file)
      }
    }
    result_path[[i]]<-value
  }
  
  emp_power<-c();emp_loc<-c()
  graph_power<-c();graph_loc<-c()
  energy_power<-c();energy_loc<-c()
  mmd_power<-c();mmd_loc<-c()
  
  for (j in 1:length(type_seq)){
    dec_depth<-c();loc_depth<-c()
    dec_graph<-c();loc_graph<-c()
    dec_ecp<-c();loc_ecp<-c()
    dec_mmd<-c();loc_mmd<-c()
    
    
    for (i in 1:length(result_path[[j]])){
      
      r <- get(load(result_path[[j]][i]))
      if (length(r)==4){
        dec_depth<-rbind(dec_depth,unlist(r$depth$p_val)<=0.05)
        loc_depth<-rbind(loc_depth,unlist(r$depth$loc))
        
        dec_graph<-rbind(dec_graph,unlist(r$graph$pval.appr$generalized)<=0.05)
        loc_graph<-rbind(loc_graph,unlist(r$graph$scanZ$generalized$tauhat))
        
        dec_ecp<-rbind(dec_ecp,unlist(length(r$ecp$estimates))>2)
        ecp_check<-function(x) {
          if (length(x$ecp$estimates)>2)
          {x$ecp$order.found[3]} else{x$ecp$considered.last[1]}
        }
        loc_ecp<-rbind(loc_ecp,unlist(ecp_check(r)))
        
        
        dec_mmd<-rbind(dec_mmd,unlist(r$mmd[1,2]>=quantile(r$mmd[2:dim(r$mmd)[1],2],0.95)))
        loc_mmd<-rbind(loc_mmd,unlist(r$mmd[1,1]))
      }else{
        dec_depth<-rbind(dec_depth,unlist(lapply(r, function(x) x$depth$p_val))<=0.05)
        loc_depth<-rbind(loc_depth,unlist(lapply(r, function(x) x$depth$loc)))
        
        dec_graph<-rbind(dec_graph,unlist(lapply(r, function(x)x$graph$pval.appr$generalized))<=0.05)
        loc_graph<-rbind(loc_graph,unlist(lapply(r, function(x)x$graph$scanZ$generalized$tauhat)))
        
        dec_ecp<-rbind(dec_ecp,unlist(lapply(r, function(x) length(x$ecp$estimates)))>2)
        loc_ecp<-rbind(loc_ecp,unlist(lapply(r, function(x) if (length(x$ecp$estimates)>2){x$ecp$order.found[3]}else{x$ecp$considered.last[1]})))
        
        dec_mmd<-rbind(dec_mmd,unlist(lapply(r, function(x) x$mmd[1,2]>=quantile(x$mmd[2:dim(x$mmd)[1],2],0.95))))
        loc_mmd<-rbind(loc_mmd,unlist(lapply(r, function(x) x$mmd[1,1])))}
    }
    
    #mean(x$mmd[1,2]>=quantile(r[[1]]$mmd[2:501,2],0.95))
    
    emp_power<-rbind(emp_power,colMeans(dec_depth))#each column is one case, each row is one parameter setting, i.e, one delta value
    emp_loc<-rbind(emp_loc,colMeans(abs(loc_depth-100)))
    
    graph_power<-rbind(graph_power,colMeans(dec_graph))
    graph_loc<-rbind(graph_loc,colMeans(abs(loc_graph-100)))
    
    energy_power<-rbind(energy_power,colMeans(dec_ecp))#
    energy_loc<-rbind(energy_loc,colMeans(abs(loc_ecp-100)))
    
    mmd_power<-rbind(mmd_power,colMeans(dec_mmd))#
    mmd_loc<-rbind(mmd_loc,colMeans(abs(loc_mmd-100)))
  }
  
  
  power_to_df<-function(power,type_seq){
    df<-as.data.frame(power) 
    if (ncol(power)==1){
      names(df)<-c("result")
    }else if(set %in% c(1,2,3)){names(df)<-c("p = 30","p = 90","p = 180")
    }else if(set %in% c(4)){names(df)<-c("p = 5","p = 15","p = 60")}
    
    df$delta<-type_seq
    df <- reshape2::melt(df ,  id.vars = 'delta', variable.name = 'case')
    df
  }
  
  #title_pool<-c(bquote("Mean difference" ~ Delta[1]),bquote("Variance difference" ~ Delta[2]),bquote("Mean shift of Gaussian components" ~ Delta[3]),bquote("D.f. of t distribution v" ),bquote("Mean shift of the distribution of mean" ~ delta[1]), bquote("Scale change of the distribution of mean" ~ delta[2]),bquote("Network parameter" ~ gamma))
  title_pool<-c(expression(Delta[1]),expression(Delta[2]),expression(Delta[3]),expression("v" ),expression(delta[1]), expression(delta[2]),expression(gamma))
  title=title_pool[set]
  
  #expression(paste(title, Delta[1]))
  ld<-lapply(list(emp_power,graph_power,energy_power,mmd_power),power_to_df,type_seq=type_seq)
  result_df<-bind_rows(ld, .id = "column_label")
  result_df$method=rep(c("dist-CP","graph-CP","energy-CP","kernel-CP"),each=nrow(emp_power)*ncol(emp_power))
  names(result_df)<-c("column_label","delta","case","value","Test")
  p <- ggplot(result_df,aes(x=delta,y=value, colour = Test, group = Test)) + geom_line(aes(lty=Test),linewidth=1.3)+facet_grid(cols = vars(case))
  #ggplot_build(p)$data
  p<-p+xlab(title) + ylab("Power")+geom_hline(yintercept=0.05,linetype=3)+theme_bw()+theme(panel.grid.minor = element_line(linetype="dashed"))+theme(panel.grid.major = element_line(linetype="dashed"))+
    theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+theme(legend.text=element_text(size=12),legend.title=element_text(size=16),legend.position = "top",legend.key.width = unit(1., 'cm'))
  if (set %in% c(1,2,3,4)){
    p<-p+theme(strip.text=element_text(size=13),legend.position = "right")+scale_color_manual(values = c("dist-CP" = "#F8766D", "energy-CP" = "#7CAE00", "graph-CP" ="#00BFC4","kernel-CP" = "#C77CFF"))+
      scale_linetype_manual(values=c("solid", "22","42","44"))
  }else {p<-p+theme(strip.background = element_blank(), strip.text = element_blank(),legend.position = "right")+scale_color_manual(values = c("dist-CP" = "#F8766D", "energy-CP" = "#7CAE00", "graph-CP" ="#00BFC4","kernel-CP" = "#C77CFF"))+
    scale_linetype_manual(values=c("solid", "22","42","44"))}
  if (set %in% c(1,2,3,4)){
    ggsave(paste(save_name, "power.pdf", sep="_"), plot=p,width = 15, height = 5,units="in")
    
  }else {ggsave(paste(save_name, "power.pdf", sep="_"), plot=p,width = 7, height = 5,units="in")}
  
  
  
  ld<-lapply(list(emp_loc,graph_loc,energy_loc,mmd_loc),power_to_df,type_seq=type_seq)
  result_df<-bind_rows(ld, .id = "column_label")
  result_df$value=result_df$value/300
  result_df$method=rep(c("dist-CP","graph-CP","energy-CP","kernel-CP"),each=nrow(emp_power)*ncol(emp_power))
  names(result_df)<-c("column_label","delta","case","value","Test")
  p <- ggplot(result_df,aes(x=delta,y=value, colour = Test, group = Test)) + geom_line(aes(lty=Test),linewidth=1.3)+facet_grid(cols = vars(case))
  
  p<-p+xlab(title) + ylab("MAE")+theme_bw()+theme(panel.grid.minor = element_line(linetype="dashed"))+theme(panel.grid.major = element_line(linetype="dashed"))+
    theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+theme(legend.text=element_text(size=12),legend.title=element_text(size=16),legend.position = "top",legend.key.width = unit(1., 'cm'))
  
  if (set %in% c(1,2,3,4)){
    p<-p+theme(strip.text=element_text(size=13),legend.position = "right")+
      scale_color_manual(values = c("dist-CP" = "#F8766D", "energy-CP" = "#7CAE00", "graph-CP" ="#00BFC4","kernel-CP" = "#C77CFF"))+
      scale_linetype_manual(values=c("solid", "22","42","44"))
  }else {p<-p+theme(strip.background = element_blank(), strip.text = element_blank(),legend.position = "right")+
    scale_color_manual(values = c("dist-CP" = "#F8766D", "energy-CP" = "#7CAE00", "graph-CP" ="#00BFC4","kernel-CP" = "#C77CFF"))+
    scale_linetype_manual(values=c("solid", "22","42","44"))
  }
  if (set %in% c(1,2,3,4)){
    ggsave(paste(save_name, "loc.pdf", sep="_"), plot=p,width = 15, height = 5,units="in")
    
  }else {ggsave(paste(save_name, "loc.pdf", sep="_"), plot=p,width = 7, height = 5,units="in")}
  
}



