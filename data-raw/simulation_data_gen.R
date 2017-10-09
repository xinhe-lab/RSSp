library(SeqArray)
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

seqExport(raw_haplo_gds,sub_haplof)
seqExport(raw_geno_gds,sub_genof)
