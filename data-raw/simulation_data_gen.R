library(ldshrink)

data("reference_genotype")
data("reference_map")

evdR <- ldshrink::ldshrink_evd(reference_genotype,map=reference_map)
panel_eigenvalues <- evdR$D
usethis::use_data(panel_eigenvalues)
