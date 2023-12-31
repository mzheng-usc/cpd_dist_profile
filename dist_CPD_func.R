source("objTest_fctns.R")
#require(MASS)


########## This is the function to run the depth change point detection algorith
#' @param distmat The pairwise distance matrix 
#' @param c the cutoff for the interval, I_c
#' @param num_permut number of permutation
dist_CPD<-function(distmat,c=0.1,num_permut=500){
  n<-nrow(distmat)
  
  max_stat<-c()
  for (j in 0:num_permut){
    #data permutation
    cat(j,"th iteration","\n")
    if (j!=0){
      set.seed(j)
      ind<-sample(nrow(distmat))
      distmat<-distmat[ind,ind]
    }
    ######## calculate the scan statistics
    test<-c()
    for (cp in seq(ceiling(n*c),n-ceiling(n*c),1)){
      testStat <- getT( distmat = distmat, indices = 1:n, n = cp, m = n-cp, cut_off = 0 )
      test<-cbind(test,testStat)
    }
    ######## check the maximum scan statistics and its index
    max_stat<-cbind(max_stat,max(test))
    if (j==0){
      max_stat_index<-which.max(test)+ceiling(n*c)
    }
  }
  p_val<-(1+sum(max_stat[,1]<max_stat[,2:ncol(max_stat)]))/(1+num_permut)
  loc<-max_stat_index
  
  result=list()
  result[['p_val']]=p_val
  result[['loc']]=loc
  result[['observed_test_statistics']]=max_stat[1]
  #result[['permuted_test_statistics']]=max_stat[-1]
  result
}



dist_CPD_permutated<-function(distmat,c=0.1,num_permut=500){
  n<-nrow(distmat)
  
  max_stat<-c()
  for (j in 0:num_permut){
    #data permutation
    cat(j,"th iteration","\n")
    if (j!=0){
      set.seed(j)
      ind<-sample(nrow(distmat))
      distmat<-distmat[ind,ind]
    }
    ######## calculate the scan statistics
    test<-c()
    for (cp in seq(ceiling(n*c),n-ceiling(n*c),1)){
      testStat <- getT( distmat = distmat, indices = 1:n, n = cp, m = n-cp, cut_off = 0 )
      test<-cbind(test,testStat)
    }
    ######## check the maximum scan statistics and its index
    max_stat<-cbind(max_stat,max(test))
    if (j==0){
      max_stat_index<-which.max(test)+ceiling(n*c)
    }
  }
  p_val<-(1+sum(max_stat[,1]<max_stat[,2:ncol(max_stat)]))/(1+num_permut)
  loc<-max_stat_index
  
  result=list()
  result[['p_val']]=p_val
  result[['loc']]=loc
  result[['observed_test_statistics']]=max_stat[1]
  result[['permuted_test_statistics']]=max_stat[-1]
  result
}

