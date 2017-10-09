#Get ldsc example data in a format we can work with

library(tidyverse)
scz_f <- "~/Dropbox/ldsc/examples/pgc.cross.SCZ17.2013-05.txt"
#scz_df <- read_delim(scz_f,delim="\t")
sczf <- "~/Dropbox/ldsc/examples/scz.sumstats.gz"

ld_dir <- "~/Dropbox/ldsc/examples/eur_w_ld_chr"
ld_files <- dir(ld_dir,full.names = T,pattern = "gz$")
all_ld <- bind_rows(map(ld_files,read_delim,delim="\t"))

geno2_info <- read_df_h5("/home/nwknoblauch/Dropbox/LDshrink/inst/snakemake_files/genotype.h5","snp_info")
colnames(geno2_info) <- c("pos","ref","alt")
geno2_info <-   mutate(geno2_info,
                       ref=intToUtf8(ref,T),alt=intToUtf8(alt,T))



# oscz_df <- read_delim(sczf,"\t")
# goscz_df <- filter(oscz_df,!is.na(N))
# ld_scz <- inner_join(all_ld,goscz_df)
# 
# bscz <- inner_join(scz_df,ld_scz,by=c("snpid"="SNP"))
# 
# 
# 
# lmm <- lm(Z~L2,data=bscz)
# plot(lmm)
# 
# ggplot(bscz,aes(x=Z,y=log(or)/se))+geom_point()

library(SeqArray)
gdsf <- "/media/nwknoblauch/Data/1kg/gds/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.gds"
snpgdsf <- "/media/nwknoblauch/Data/1kg/gds/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.snp_gds"

gdsm <-  seqOpen(gdsf,readonly = T)
# sgdsm <- snpgdsOpen(snpgdsf)
seqGDS2SNP(gdsm,snpgdsf,compress.geno = "LZ4_RA.fast")
seqClose(gdsm)
sgdsm <- snpgdsOpen(gdsf,readonly=T)
hap_legend_f <- "/media/nwknoblauch/Data/1kg/haplotypes/EUR.chr1_1kg_geno.legend.gz"
hap_data_f <- "/media/nwknoblauch/Data/1kg/haplotypes/EUR.chr1_1kg_geno.hap.gz"
EUR_ind <- scan("/media/nwknoblauch/Data/1kg/EUR.samples",what=character())

hap_d <- read_delim(hap_legend_f,delim=" ")
t_hap_d <- slice(hap_d,1:10000)
hap_data <- matrix(scan(hap_data_f,what=integer(),sep = " ",nlines = nrow(t_hap_d)),nrow(t_hap_d),byrow = T)





all_samps <- read.gdsn(index.gdsn(gdsm,"/sample.id"))
all_vars <- read.gdsn(index.gdsn(gdsm,"/variant.id"))
vid <- read.gdsn(index.gdsn(gdsm,"/annotation/id"))
t_hap_ind <- match(t_hap_d$id,vid)
t_hap_d <- mutate(t_hap_d,match_ind=t_hap_ind,t_h_ind=1:n())
seqSetFilter(gdsm,variant.id=t_hap_ind,sample.id=EUR_ind)
hapdat <- seqGetData(gdsm,var.name = "genotype")
h_mat <- matrix(c(hapdat[]),nrow(hap_data),ncol(hap_data),byrow = T)





