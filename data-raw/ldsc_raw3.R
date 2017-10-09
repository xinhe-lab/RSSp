#ldscore_raw3

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
add.gdsn(index.gdsn(genofile,"annotation"), "genetic.map", mapdf$snp.map)
closefn.gds(genofile)

genofile <- seqOpen(gds.fn = gds_file,readonly = T)
hapdat <- matrix(as.numeric(seqGetData(genofile,var.name = "genotype")[]),nrow(mapdf),byrow = T)
map <- read.gdsn(index.gdsn(genofile,'annotation/genetic.map'))

big_R <- LDshrink::calcLD(t(hapdat),map)
ldsc_R <- 
R_df <- data_frame(R_LD=ldsc_R,SNP=read.gdsn(index.gdsn(genofile,'annotation/id')))

ldsc_df <- read_delim("~/Dropbox/ldsc/examples/1kg_eur/22.l2.ldscore.gz",delim="\t")

osumf <- "~/Dropbox/ldsc/examples/scz.sumstats.gz"