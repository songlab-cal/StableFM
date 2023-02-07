## Combine csvs and convert to perturbation scores
library(dplyr)

# Read processed track prediction CSVs
# Convert difference to magnitude
stab.diff.pred.tss.relevant <-
  read.csv('enformer-data/stab_diff_pred_tss_relevant.csv') %>% abs()
top.diff.pred.tss.relevant <-
  read.csv('enformer-data/top_diff_pred_tss_relevant.csv') %>% abs()
stab.diff.pred.ave.relevant <-
  read.csv('enformer-data/stab_diff_pred_ave_relevant.csv') %>% abs()
top.diff.pred.ave.relevant <-
  read.csv('enformer-data/top_diff_pred_ave_relevant.csv') %>% abs()

stab.diff.pred.tss.relevant.new <-
  read.csv('enformer-data/stab_diff_pred_tss_relevant_new.csv') %>% abs()
top.diff.pred.tss.relevant.new <-
  read.csv('enformer-data/top_diff_pred_tss_relevant_new.csv') %>% abs()
stab.diff.pred.ave.relevant.new <-
  read.csv('enformer-data/stab_diff_pred_ave_relevant_new.csv') %>% abs()
top.diff.pred.ave.relevant.new <-
  read.csv('enformer-data/top_diff_pred_ave_relevant_new.csv') %>% abs()

# Read variant ID files
message(date(), ": Match track data with SNP rsIDs")
stab.snps.strand.final.variants <- read.csv('enformer-data/stab_snps_strand_final_variants.csv')
top.snps.strand.final.variants <- read.csv('enformer-data/top_snps_strand_final_variants.csv')

stab.snps.strand.final.variants.new <- read.csv('enformer-data/stab_snps_new_strand_final_variants.csv')
top.snps.strand.final.variants.new <- read.csv('enformer-data/top_snps_new_strand_final_variants.csv')

# Fill in IDs 
stab.diff.pred.tss.relevant$GENE <- stab.snps.strand.final.variants$GENE
stab.diff.pred.tss.relevant$rsID <- stab.snps.strand.final.variants$rsID
stab.diff.pred.tss.relevant$REF <- stab.snps.strand.final.variants$REF
stab.diff.pred.tss.relevant$ALT <- stab.snps.strand.final.variants$ALT
top.diff.pred.tss.relevant$GENE <- top.snps.strand.final.variants$GENE
top.diff.pred.tss.relevant$rsID <- top.snps.strand.final.variants$rsID
top.diff.pred.tss.relevant$REF <- top.snps.strand.final.variants$REF
top.diff.pred.tss.relevant$ALT <- top.snps.strand.final.variants$ALT

stab.diff.pred.tss.relevant.new$GENE <- stab.snps.strand.final.variants.new$GENE
stab.diff.pred.tss.relevant.new$rsID <- stab.snps.strand.final.variants.new$rsID
stab.diff.pred.tss.relevant.new$REF <- stab.snps.strand.final.variants.new$REF
stab.diff.pred.tss.relevant.new$ALT <- stab.snps.strand.final.variants.new$ALT
top.diff.pred.tss.relevant.new$GENE <- top.snps.strand.final.variants.new$GENE
top.diff.pred.tss.relevant.new$rsID <- top.snps.strand.final.variants.new$rsID
top.diff.pred.tss.relevant.new$REF <- top.snps.strand.final.variants.new$REF
top.diff.pred.tss.relevant.new$ALT <- top.snps.strand.final.variants.new$ALT

stab.diff.pred.ave.relevant$GENE <- stab.snps.strand.final.variants$GENE
stab.diff.pred.ave.relevant$rsID <- stab.snps.strand.final.variants$rsID
stab.diff.pred.ave.relevant$REF <- stab.snps.strand.final.variants$REF
stab.diff.pred.ave.relevant$ALT <- stab.snps.strand.final.variants$ALT
top.diff.pred.ave.relevant$GENE <- top.snps.strand.final.variants$GENE
top.diff.pred.ave.relevant$rsID <- top.snps.strand.final.variants$rsID
top.diff.pred.ave.relevant$REF <- top.snps.strand.final.variants$REF
top.diff.pred.ave.relevant$ALT <- top.snps.strand.final.variants$ALT

stab.diff.pred.ave.relevant.new$GENE <- stab.snps.strand.final.variants.new$GENE
stab.diff.pred.ave.relevant.new$rsID <- stab.snps.strand.final.variants.new$rsID
stab.diff.pred.ave.relevant.new$REF <- stab.snps.strand.final.variants.new$REF
stab.diff.pred.ave.relevant.new$ALT <- stab.snps.strand.final.variants.new$ALT
top.diff.pred.ave.relevant.new$GENE <- top.snps.strand.final.variants.new$GENE
top.diff.pred.ave.relevant.new$rsID <- top.snps.strand.final.variants.new$rsID
top.diff.pred.ave.relevant.new$REF <- top.snps.strand.final.variants.new$REF
top.diff.pred.ave.relevant.new$ALT <- top.snps.strand.final.variants.new$ALT

# ombine new and old into a single dataframe
stab.diff.pred.ave.relevant.comb <- rbind(stab.diff.pred.ave.relevant,
                                          stab.diff.pred.ave.relevant.new)
top.diff.pred.ave.relevant.comb <- rbind(top.diff.pred.ave.relevant,
                                         top.diff.pred.ave.relevant.new)

stab.diff.pred.tss.relevant.comb <- rbind(stab.diff.pred.tss.relevant,
                                          stab.diff.pred.tss.relevant.new)
top.diff.pred.tss.relevant.comb <- rbind(top.diff.pred.tss.relevant,
                                         top.diff.pred.tss.relevant.new)

# Save
save(top.diff.pred.tss.relevant.comb, file="top_diff_tss_relevant_comb.rds")
save(stab.diff.pred.tss.relevant.comb, file="stab_diff_tss_relevant_comb.rds")
save(top.diff.pred.ave.relevant.comb, file="top_diff_ave_relevant_comb.rds")
save(stab.diff.pred.ave.relevant.comb, file="stab_diff_ave_relevant_comb.rds")

