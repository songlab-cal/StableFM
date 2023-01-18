###################
# genome-wide PCA #
###################

## Option 1
geno_dir <- "/Users/alanaw/Documents/research/pgs/082821/genotypes/"
gw_data <- matrix(nrow=445)
for (chr in 1:22) {
  genotype_data <- read.csv(paste0(geno_dir, "chr", chr, "_allele_dosage.csv"))[,-1]
  gw_data <- cbind(gw_data, genotype_data)
}

## Option 2
library(bigsnpr)
test <- snp_attach("~/Documents/research/pgs/082821/1000G_phase3_common_norel/1000G_phase3_common_norel.rds")
chr<-3
genotype_data <- read.csv(paste0(geno_dir, "chr", chr, "_allele_dosage.csv"))
geuvadis.ids <- genotype_data[,1]
obj.svd <- snp_autoSVD(test$genotypes, 
                       infos.chr = test$map$chromosome,
                       ind.row = which(test$fam$sample.ID %in% geuvadis.ids), 
                       k = 5,
                       ncores = 1)
pcs.matrix <- predict(obj.svd)
norm.pcs.matrix <- apply(pcs.matrix,2,scale) %>% as.data.frame()
colnames(norm.pcs.matrix) <- paste0("PC",1:5)
norm.pcs.matrix$ID <- test$fam$sample.ID[which(test$fam$sample.ID %in% geuvadis.ids)]

save(norm.pcs.matrix, file = "gw_PC1-5.RData")
