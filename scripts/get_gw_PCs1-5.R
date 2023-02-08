## Script for performing PCA to obtain top five principal components

# Load rds object (obtained from bigsnpr)
# Dataset source: https://rdrr.io/github/privefl/mypack/src/data-raw/1000G-phase3-common-norel.R 
# Follow script above to obtain the rds/bk/bgen/bed files required
library(bigsnpr)
test <- snp_attach("~/Documents/research/pgs/082821/1000G_phase3_common_norel/1000G_phase3_common_norel.rds")

# Grab individual IDs for performing PCA
chr <- 17 # can be any autosome
genotype_data <- read.csv(paste0("../data/geno_pheno/genotypes/chr", chr, "_allele_dosage.csv")) # need to unzip the file
geuvadis.ids <- genotype_data[,1]

# PCA 
obj.svd <- snp_autoSVD(test$genotypes, 
                       infos.chr = test$map$chromosome,
                       ind.row = which(test$fam$sample.ID %in% geuvadis.ids), 
                       k = 5,
                       ncores = 1)
pcs.matrix <- predict(obj.svd)
norm.pcs.matrix <- apply(pcs.matrix,2,scale) %>% as.data.frame()
colnames(norm.pcs.matrix) <- paste0("PC",1:5)
norm.pcs.matrix$ID <- test$fam$sample.ID[which(test$fam$sample.ID %in% geuvadis.ids)]

# Save data
save(norm.pcs.matrix, file = "../data/geno_pheno/metadata/PCA/gw_PC1-5.RData")
