# context("gds tests")
# 
# test_that("Simulating from gds works the same as simulating from SNP data",{
#   library(SeqArray)
#   library(dplyr)
#   # library(RSSp)
#   gds.fn <- seqExampleFileName("gds")
#   
#   # display
#   gds <- seqOpen(gds.fn)
#   seqSetFilterCond(gds,maf = 0.01)  
#   aSNP <- seqGetData(gds,"$dosage")
#   snpi <- LDshrink::read_SNPinfo_gds(gds)
#   good_SNP <- apply(aSNP,2,function(x)all(!is.na(x)))
#   snpi <- filter(snpi,good_SNP)
# 
#   seqSetFilter(gds,variant.id = snpi$snp_id,action="intersect")
#   mseed <- c(123L)
#   nSNP <- scale(seqGetData(gds,"$dosage"),center=T,scale=F)
#   pve <- 0.1
#   bias <- 0
#   nreps <- 2
#   p <- length(seqGetData(gds,"variant.id"))
#   
#   #First try with chunksize>=p
#   R_list <- gen_sim_list(nSNP, pve, bias=0, nreps,seed=mseed)
#   gds_list <- gen_sim_gds(gds,pve,bias=0,nreps,seed=mseed,chunksize = p)
#   expect_equal(R_list,gds_list)
# 
#   expect_equal(R_list,gds_list)
# })


