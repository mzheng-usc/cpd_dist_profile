# [Change Point Detection for Random Objects using Distance Profiles](https://arxiv.org/abs/2311.16025)

This repository hosts the code base for the paper

**Change Point Detection for Random Objects using Distance Profiles**

Paromita Dubey, Minxing Zheng

Department of Data Sciences and Operations, Marshall School of Business, University of Southern California

[Link to Preprint](https://arxiv.org/abs/2311.16025)

## Overview

The core function for implementing the proposed distance profile based change point detection method is the `depth_CPD` function in *depth_CPD_func.R*. The `depth_CPD` function takes a pairwise distance matrix as input and outputs the p-value, estimated change point location and test statistics. 

For example, we generate a data sequence of length 100, with first 50 observations generated from $N(0_p,I)$  and last 50 observations generated from $N(1_p,I)$  , where $p$ Is the dimension of the multivariate random vector.

```R
library(MASS)
source("depth_CPD_func.R") 
p <- 10 # dimension of the multivariate random vector
num_permut <- 100 #number of permutation for the permutation test
I<-diag(x = 1, p, p)

#generated data sequence
Data<-rbind(mvrnorm(50,mu=rep(0,p),Sigma=I),mvrnorm(50,mu=c(rep(1,p)),Sigma=I))
#caculate the pariwise distance matrix, with eculidean distance
distmat<-as.matrix(dist( Data, method = 'euclidean' ))

#the depth_CPD function takes pairwise distance matrix as input, and user could also set number of permutation for permutation test and a cut-off parameter c.
depth_result<-depth_CPD(distmat,num_permut =num_permut,c=0.1)
```

The output will look like the following, the p value is less than 0.01, and the estimated change point location is 51.

```R
> depth_result
$p_val
[1] 0.00990099

$loc
[1] 51

$observed_test_statistics
[1] 8.049063
```



