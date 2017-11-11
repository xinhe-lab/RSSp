library(SeqArray)
library(SeqSupport)
library(dplyr)
rssp_dir <- "inst/gds"
raw_haplof <- "/media/nwknoblauch/Data/1kg/gds/ALL_seq_hapmap.gds"
raw_genof <- "/media/nwknoblauch/Data/1kg/gds/ALL_seq_hapmap_geno.gds"


raw_haplo_gds <- seqOpen(raw_haplof)
raw_geno_gds <- seqOpen(raw_genof) 

seqSetFilterChrom(raw_haplo_gds,"19")

hap_leg <- read_SNPinfo_gds(raw_haplo_gds,region_id = T)
nreg <- unique(hap_leg$region_id)

tl <- subset_gds(raw_haplo_gds,filter(hap_leg,region_id %in% nreg[1:2]))
ntl <- subset_gds(raw_geno_gds,filter(hap_leg,region_id %in% nreg[1:2]))

rssp_dir <- "inst/gds"
sub_haplof <- file.path(rssp_dir,"sub_19_haplo.gds")
sub_genof <- file.path(rssp_dir,"sub_19_geno.gds")
                                        #genogds is the SeqArray file with the data stored as genotypes (it's used here for GWAS simulation)
genogds <- SeqArray::seqOpen(sub_genof)
                                        #haplogds is the SeqArray file with the data stored as haplotypes (it's used to calculate LD, which is used by RSSp)
haplogds <- SeqArray::seqOpen(sub_haplof)


                                        #pve is bounded between 0 and 1 
pve <- as.numeric(seq(0.7,0.9,length.out = 10))
                                        #bias is also bounded between 0 and 1
bias <- as.numeric(seq(0.0,0.1,length.out = 2))
                                        #each pve-bias combination is simulated 3 times
nreps <- 3

                                        #count the number of SNPs
p <- calc_p(genogds)
                                        #count the number of individuals
n <- calc_N(genogds)

                                        #gen_sim_gds simulates GWAS summary stats.
msim <- gen_sim_gds(genogds,pve=pve,bias=bias,nreps=nreps)
simulation_R <- calc_LD_gds(haplogds)
devtools::use_data(simulation_R)


simulation_tparam_df<- msim[["tparam_df"]]

simulation_tparam_df <- mutate(
    simulation_tparam_df,
    n=msim$n,
    p=msim$p)
simulation_uh_mat <- msim[["bias_uh_mat"]]
devtools::use_data(simulation_tparam_df,overwrite=T)
devtools::use_data(simulation_uh_mat)



