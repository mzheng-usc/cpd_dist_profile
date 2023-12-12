require(gtools)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../functions/binary_segmentation.R") 
path<-"SBM_result"
file_list = mixedsort(list.files(path = path ,full.names=TRUE,recursive=TRUE))

minLen<-20
seed_I<-seeded.intervals(400,decay = sqrt(2))
seed_I<-seed_I[(seed_I[,2]-seed_I[,1])>minLen,]
nrow(seed_I)

depth_loc<-list();ecp_loc<-list();kcp_loc<-list()
for (i in file_list){
  d<-get(load(i))
  result_depth<-d$depth
  
  sbm_ob<-unlist(result_depth[,'observed_test_statistics'])
  sbm_loc<-unlist(result_depth[,'loc'])
  #sbm_pval<-unlist(result_depth[,'p_val'])
  sbm_threshold<-quantile(unlist(result_depth[1,"permuted_test_statistics"]),0.95)
  
  ###### A selection method that based on the quantile of permutated scan statistics on the whole sequence
  #IND contains all the indexes at first, and then remove an index such that  
  #its observed test statistics is larger than the threshold, then it is added to the sbm_depth
  #then, remove all the intervals that contain this index
  sbm_depth<-c();IND<-c(1:length(sbm_ob));I<-c()
  for (j in order(sbm_ob,decreasing = TRUE)){
    if (j %in% IND){
      if (sbm_ob[j]<sbm_threshold){
        break
        }else{sbm_depth<-c(sbm_depth,sbm_loc[j])}
    }else{next}
    ind<-which(seed_I[,1]>sbm_loc[j]|seed_I[,2]<sbm_loc[j],)
    IND<-intersect(IND, ind)
    I<-c(I,j)
  }
  depth_loc[[length(depth_loc) + 1]]<-sort(sbm_depth)
  ecp_loc[[length(ecp_loc) + 1]] <- d$ecp$estimates[-c(1,length(d$ecp$estimates))]
  kcp_loc[[length(kcp_loc) + 1]] <- d$kcp[-c(1,length(d$kcp))]
}

########## check the dimension of the result, so that we know how many cp are detected 
########## we expect the result to be of dimension 500*3 
depth_loc_mat<-t(sapply(depth_loc, "length<-", max(lengths(depth_loc))))# 500 * 3
dim(depth_loc_mat)

ecp_loc_mat<-t(sapply(ecp_loc, "length<-", max(lengths(ecp_loc))))# 500 * 5
dim(ecp_loc_mat)

kcp_loc_mat<-t(sapply(kcp_loc, "length<-", max(lengths(kcp_loc))))# 500 * 3
dim(kcp_loc_mat)

proportion_number_mcpd<-function(est_matrix){
  numer_cpd<-rowSums((!is.na(est_matrix)))
  length(numer_cpd[numer_cpd==3])/500
}
L<-list(depth_loc_mat,ecp_loc_mat,kcp_loc_mat)
#derive the proportion of correctly estimated number of change points.
sapply(L, proportion_number_mcpd)# 1.000 0.958 1.000

max_dist<-function(x,y){
  dist(rbind(x,y),method="manhattan")
}

calculate_MAE<-function(mat){
  Error1<-c()
  for (i in 1:nrow(mat)){
    x<-mat[i,]
    x<-x[!is.na(x)]
    if (length(x)==3){
      e<-dist(rbind(x,c(101,201,301)),method="manhattan")
      Error1<-c(Error1,e)
    }
  }
  mean(Error1)
}

sapply(X=list(depth_loc_mat,ecp_loc_mat,kcp_loc_mat),FUN=calculate_MAE)
#condition on estimating 3 change points,  the MAE=0 for three methods.


