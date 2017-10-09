# Download `genotype2.mat` from the following link:
#https://uchicago.app.box.com/v/example2

library(RcppEigenH5)
shrink_R <- read_2d_mat_h5("/home/nwknoblauch/Downloads/genotype2.mat","/","shrink_R")
devtools::use_data(shrink_R)
cohort_SNP <- scale(t(read_2d_mat_h5("/home/nwknoblauch/Downloads/genotype2.mat","/","C")),center=T,scale=F)
devtools::use_data(cohort_SNP)


