#Script for ingesting ldsc examples
library(SeqArray)
library(tidyverse)


 eur_url <- "https://data.broadinstitute.org/alkesgroup/LDSCORE/1kg_eur.tar.bz2"
# temp_dir <- tempdir()
 destfile <- file.path("inst/gds","1kg_eur.tar.bz2")
 download.file(eur_url,destfile = destfile)

 ntdir <- untar(tarfile = destfile,compressed="bzip2",exdir = "inst/gds")
 temp_dir <- "inst/gds"
tar_dir <- file.path(temp_dir,"1kg_eur")
bed.fn <- file.path(tar_dir,"22.bed")
fam.fn <- file.path(tar_dir,"22.fam")
bim.fn <- file.path(tar_dir,"22.bim")

# if(!dir.exists("inst/gds")){
#   dir.create("inst/gds",recursive = T)
# }
gds_file <- file.path("inst/gds","ldsc.gds")
# gds_file <- file.path("inst/gds","ldsc_seq.gds")

if(file.exists(gds_file)){
  file.remove(gds_file)
}
mgds <- seqBED2GDS(bed.fn = bed.fn,
                                 fam.fn = fam.fn,
                                 bim.fn = bim.fn,
                                 out.gdsfn = gds_file)

mapdf <- read_delim(file = bim.fn,col_names = c("snp.chromosome",
                                                "snp.id",
                                                "snp.map",
                                                "snp.pos",
                                                "snp.ref","snp.alt"),delim="\t") %>% select(snp.map)
 genofile <- openfn.gds(gds_file,readonly = F)
 add.gdsn(index.gdsn(genofile,"annotation/info/"), "genetic.map", mapdf$snp.map)
 closefn.gds(genofile)

ldsc_df <- read_delim("~/Dropbox/ldsc/examples/1kg_eur/22.l2.ldscore.gz",delim="\t")

osumf <- "~/Dropbox/ldsc/examples/scz.sumstats.gz"
tscz <- read_delim(osumf,delim="\t")
tscz <- filter(tscz,!is.na(Z))

ldsc_scz <- inner_join(ldsc_df,tscz)

genofile <- seqOpen(gds.fn = gds_file,readonly = T)
ldsc_scz <- data_frame(SNP=seqGetData(genofile,"annotation/id"),snp_id=seqGetData(genofile,"variant.id")) %>% inner_join(ldsc_scz)
seqSetFilter(genofile,variant.id=ldsc_scz$snp_id)
hapdat <- matrix(as.numeric(seqGetData(genofile,var.name = "genotype")[]),nrow(ldsc_scz),byrow = T)
map <- seqGetData(genofile,'annotation/info/genetic.map')

big_R <- LDshrink::calcLD(t(hapdat),map)
ldsc_R <- colSums(big_R^2)-1
ldsc_scz <- mutate(ldsc_scz,R_LD=ldsc_R)


# b_df <- inner_join(R_df,tscz)

##### Test that haplotype data is the with IMPUTE2 and SeqArray

hap_legf <- "/media/nwknoblauch/Data/1kg/haplotypes/EUR.chr19_1kg_geno.legend.gz"
hap_leg <- read_delim(hap_legf,delim=" ") %>% mutate(hap_id=1:n(),chr=as.character(19)) %>% rename(SNP=id,pos=position) %>% select(-a0,-a1)
hap_datf <- "/media/nwknoblauch/Data/1kg/haplotypes/EUR.chr19_1kg_geno.hap.gz"
hap_dat <- matrix(scan(hap_datf,what=integer(),sep=" "),nrow = nrow(hap_leg),byrow = T)

gdsf <- "/media/nwknoblauch/Data/1kg/gds/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_seq_hapmap.gds"

snpi <- read_SNPinfo_gds(gds)
snpi <- inner_join(snpi,hap_leg)

t_hap_dat <- hap_dat[head(snpi$hap_id),]
seqSetFilter(gds,variant.id=head(snpi$snp_id))

g_hap_dat <- matrix(c(seqGetData(gds,"genotype")),nrow=nrow(t_hap_dat),byrow=T)
all.equal(g_hap_dat,t_hap_dat)





# nb_df <- inner_join(b_df,tscz)

p <- nrow(ldsc_scz)
n <- unique(ldsc_scz$N)
tpdf <- gen_tparamdf_norm(pve=0.4,bias=0,nreps=2,n = n,p=p) %>% select(-replicate) %>% nest(fgeneid)
evdR <- eigen(big_R)
Q <- evdR$vectors
d <- evdR$values
aqh <- simuh_dir(sigu = tpdf$tsigu,bias = tpdf$tbias,Q = Q,D = d,fgeneids = tpdf$data[[1]],chr = 22,gen_quh = T,ret_d = T)
aqh <- filter(aqh,fgeneid==1)
aqh <- mutate(ldsc_scz,snp_index=1:n()) %>% inner_join(mutate(aqh,snp_index=as.integer(snp_index)))

n_ldsc_df <- select(aqh,CHR,SNP,BP,L2=R_LD)
write.table(n_ldsc_df,file="~/Dropbox/ldsc/examples/eur_w_ld_chr_test/22.l2.ldscore",col.names = T,row.names = F,sep="\t",quote=F)
ldsc_df <- aqh %>% select(SNP,A1,A2,Z=uh,N=N)
RSSp_df <- aqh 
RSSp_res <- RSSp(fgeneid = RSSp_df$fgeneid,
                 D = RSSp_df$rd,
                 quh = RSSp_df$quh,
                 p_n = nrow(RSSp_df)/RSSp_df$N[1])
write.table(ldsc_df,"~/Dropbox/ldsc/examples/test_sim.txt",col.names = T,row.names=F,sep="\t",quote=F)
h2lf <- "~/Dropbox/ldsc/examples/test2_sim_h2.log"
parse_ldsc_h2log <- function(h2lf){
  library(purrr)
  library(tidyr)
  h2_dat <- scan(h2lf,what=character(),sep = "\n")
  h2_rowi <- grep("Total Observed scale",h2_dat)
  h2_row <- h2_dat[h2_rowi]
  h2_data <- h2_dat[h2_rowi:length(h2_dat)]
  h2_data <- h2_data[grep("^[^:]+:[^:]+$",h2_data)]
  h2_datd <- transpose(strsplit(h2_data,split=":"))
  names(h2_datd) <- c("Variable","Value")
  h2_datdf <- unnest(as_data_frame(h2_datd)) %>% mutate(Variable=chartr(" ^","__",Variable),Value=trimws(Value)) %>% separate(Value,c("Est","SD"),sep = "[ s\\(\\)]+",remove=T,fill="right",convert=T)
  return(h2_datdf)
}








