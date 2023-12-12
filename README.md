# [Change Point Detection for Random Objects using Distance Profiles](https://arxiv.org/abs/2311.16025)

This repository hosts the code base for the paper

**Change Point Detection for Random Objects using Distance Profiles**

Paromita Dubey, Minxing Zhen

Department of Data Sciences and Operations, Marshall School of Business, University of Southern California

[Link to Paper](https://arxiv.org/abs/2311.16025)

## Overview

We store all the core functions for implementing the proposed distance profile based change point detection method along as other alternative methods including energy-based and kernel-based methods in the **functions**. To be more specific,

* depth_CPD_func.R is used to implement the proposed distance profile based change point detection method, the *depth_CPD* function takes a pairwise distance matrix as input and outputs the p-value, estimated change point location and test statistics.
* ecp_distmat_input.R is used to implement the energy-based method (modified based on R package "ecp"), the *e.divisive_distmat* function only takes a pairwise distance matrix as require input and outputs information according to the *e.divisive* function in the "ecp" package.
* Kernel_MMD.R is  used to implement the kernel-based single change point detection (modified based on python package "“Chapydette"), the *MMD_test* function only takes a pairwise distance matrix as require input and the outputs the objective value at each location. The estimated change point location occurs at the maximum of objective value.
* Kip_distmat_input.R is used to implement the kernel-based method (modified based on R package "ecp"), the *kcpa_distmatt* function takes a pairwise distance matrix and a parameter indicating the maximum number of change points as require input and outputs the estimated change point locations.

**real_data** stores all the code and data to process the raw data and implement the change point detection methods on two real-world datasets including:  U.S. electricity generation dataset and MIT reality mining dataset.

**simulation** stores all the source code to implement various settings of simulations on the USC cluster.

**result** stores the code to process the output files generated by simulation and real data.
