############################################
########### Auxiliary functions ############
############################################
## File storage details
# 1. Fine-mapping results 
#    - new results (regularization approach + stability-guided approach)
#      /clusterfs/nilah/alan/PICS2/PC_correction_results/
#    - old results (DO NOT USE!)
#      /clusterfs/nilah/alan/PICS2/peer_results/
#    - combined results with additional info
#      /clusterfs/nilah/alan/PICS2/fine-mapping-results/
# 2. GEUVADIS Phenotypes
#    - peer factor normalized 
#      /clusterfs/nilah/alan/peer_phenotypes/
#    - peer factor normalized + PC regressed out phenotypes are not saved,
#      they are computed in scripts like PICS2_chr22_rm_PCs.R
#    - metadata of phenotypes 
#      /clusterfs/nilah/alan/peer_phenotypes/E-GEUV-1_subset_pop_anno.txt
# 3. Annotations
#    /clusterfs/nilah/alan/all-func-annots/
#    - Distance to TSS
#    - Ensembl
#    - FAVOR
#    - FIRE 
#    - Enformer (x 177)
# 4. Genotypes
#    - tbd 

# agg_res_dir <- '/clusterfs/nilah/alan/PICS2/fine-mapping-results/'
# pheno_dir <- '/clusterfs/nilah/alan/peer_phenotypes/'
# pheno_metadata <- '/clusterfs/nilah/alan/peer_phenotypes/E-GEUV-1_subset_pop_anno.txt'
# annot_dir <- '/clusterfs/nilah/alan/all-func-annots/'

agg_res_dir <- 'PICS2/fine-mapping-results/'
pheno_dir <- 'peer_phenotypes/'
pheno_metadata <- 'peer_phenotypes/E-GEUV-1_subset_pop_anno.txt'
annot_dir <- 'all-func-annots/'

gene_dfs <- lapply(1:22, readFile <- function(x) {
  read.csv(paste0('peer_phenotypes/chr', x, '_gene_names.csv'))
})
names(gene_dfs) <- 1:22

#' summmarizeGeneExp 
#' 
#' This function takes in user-provided chr and gene name,
#' then proceeds with generating a plot summarizing the 
#' distribution of PEER-corrected gene expression across the five slices
#' 
summarizeGeneExp <- function(chr, input.gene.name) {
  # Get slice annotations
  pop_metadata <- data.table::fread(paste0(pheno_dir,
                                           'E-GEUV-1_subset_pop_anno.txt'))
  
  # Get gene expression matrix for all chr's phenotypes
  phenos_df <- data.table::fread(paste0(pheno_dir, 
                                          'E-GEUV-1_GeneQuantRPKM_chr', 
                                          chr,
                                          '_subset_new.txt'))
  # Convert input gene name to Ensembl gene name
  gene.name <- (gene_dfs[[chr]] %>% 
                  subset(INPUT == input.gene.name))[['ensembl_version_id']]
  
  # Filter to specified gene, then combine with slice annotations
  gene_exp_df <- phenos_df %>% 
    select(gene.name) %>%
    cbind(pop_metadata)
  
  # Rename first column 
  colnames(gene_exp_df)[1] <- 'VALUE'
  
  # Extract pop slice names
  slice_names  <- colnames(gene_exp_df)[-c(1,2)]
  
  # Generate factor feature for plotting
  factor_feature_matrix <- gene_exp_df %>% select(slice_names)
  factor_feature_df <- data.frame(SLICE = sapply(apply(factor_feature_matrix,1,function(x) which(x==1)),
                                                 function(x) {slice_names[x]}),
                                  VALUE = gene_exp_df$VALUE)
  
  # Generate plot
  multi_slice_plot <- ggplot(factor_feature_df, aes(x=SLICE,y=VALUE)) +
    #geom_boxplot(aes(fill = SLICE), alpha= 0.5) +
    geom_violin(aes(fill = SLICE), alpha= 0.5) +
    theme_bw() +
    ggtitle('Distribution of Gene Expression Phenotype') +
    ylab('Value') + xlab('') +
    labs(fill = 'Slice') +
    theme(title=element_text(size=rel(1.2)),
          axis.text.x=element_text(size=rel(1.7), angle=90),
          axis.title.y=element_text(size=rel(1.4)),
          legend.title=element_text(size=rel(1.35)),
          legend.text=element_text(size=rel(1.3)))
  
  # Generate slice-specific statistics
  mean_vec <- c(); se_vec <- c()
  for (slice in slice_names) {
    slice_df <- factor_feature_df %>% subset(SLICE == slice)
    mean_vec <- c(mean_vec, mean(slice_df$VALUE))
    se_vec <- c(se_vec, sd(slice_df$VALUE))
  }
  sum_table <- data.frame(MEAN = mean_vec,
                          SE = se_vec,
                          row.names = slice_names)
  colnames(sum_table) <- c('Mean', 'Standard Deviation')
  
  # Return list of plot and sum_table
  return(list(PLOT = multi_slice_plot,
              SUM_TABLE = sum_table,
              ENSEMBL_GENE = gene.name,
              GENE_SYMBOL = (gene_dfs[[chr]] %>% 
                               subset(INPUT == input.gene.name))[['SYMBOL']],
              BIOTYPE = (gene_dfs[[chr]] %>% 
                           subset(INPUT == input.gene.name))[['GENEBIOTYPE']]))
}

#' summmarizePICSResult
#' 
#' This function takes in user-provided chromosome, gene name, 
#' and potential set, then proceeds with generating a table summarizing
#' key moderating characteristics for the fine-mapped variant
#'
#' Example data for testing functionality:
#' chr <- 22
#' gene.name <- gene_names[[22]][1]
#' ps <- 1
summarizePICSResult <- function(chr, input.gene.name, ps) {
  # Get fine-mapping result df
  pics2.summary.df <- readr::read_csv(paste0(agg_res_dir, 
                                             chr, '/all_pp_max.csv')) 
  pics2.summary.df$GENE <- sapply(pics2.summary.df$GENE,
                                  function(x){
                                    stringr::str_split(x, '_')[[1]][2]
                                  })
  # Convert input gene name to Ensembl gene name
  gene.name <- (gene_dfs[[chr]] %>% 
                  subset(INPUT == input.gene.name))[['ensembl_version_id']]
  
  # Subset to relevant gene and ps
  rel.gene.ps <- pics2.summary.df %>% 
    subset(GENE == gene.name) %>% 
    select(c(paste0('TOP_',ps),
             paste0('TOP_',ps,'_VAR_HG19'),
             paste0('TOP_',ps,'_VAR_HG38'),
             paste0('TOP_',ps,'_SUPPORT'),
             paste0('POSTPROB_',ps),
             paste0('STAB_',ps),
             paste0('STAB_',ps,'_VAR_HG19'),
             paste0('STAB_',ps,'_VAR_HG38'),
             paste0('STAB_', ps,'_SUPPORT'),
             paste0('STAB_POSTPROB_', ps),
             paste0('FST_',ps),
             paste0('AF_',ps),
             paste0('NSETS_',ps),
             paste0('YRI_TOP_',ps),
             paste0('YRI_STAB_',ps))) 
  
  # Create summary dataframe to return
  # Need to include REF/ALT later
  top.vs.stable.df <- data.frame(Top = c(rel.gene.ps[[paste0('TOP_',ps)]],
                                         rel.gene.ps[[paste0('TOP_',ps,'_VAR_HG19')]],
                                         rel.gene.ps[[paste0('TOP_',ps,'_VAR_HG38')]],
                                         round(rel.gene.ps[[paste0('POSTPROB_',ps)]],3),
                                         rel.gene.ps[[paste0('TOP_',ps,'_SUPPORT')]],
                                         rel.gene.ps[[paste0('YRI_TOP_',ps)]]),
                                 Stable = c(rel.gene.ps[[paste0('STAB_',ps)]],
                                            rel.gene.ps[[paste0('STAB_',ps,'_VAR_HG19')]],
                                            rel.gene.ps[[paste0('STAB_',ps,'_VAR_HG38')]],
                                            round(rel.gene.ps[[paste0('STAB_POSTPROB_',ps)]],3),
                                            rel.gene.ps[[paste0('STAB_',ps,'_SUPPORT')]],
                                            rel.gene.ps[[paste0('YRI_STAB_',ps)]]),
                                 row.names = c('rsID',
                                               'POS:REF:ALT (GRCh37 build)',
                                               'POS:REF:ALT (GRCh38 build)',
                                               'Posterior Probability',
                                               'Positive Posterior Probability Support',
                                               'Variant has Positive Probability in YRI Slice?'))
  stable.stats.df <- data.frame(FST = rel.gene.ps[[paste0('FST_',ps)]],
                                AF = rel.gene.ps[[paste0('AF_',ps)]],
                                NSLICES = as.integer((rel.gene.ps[[paste0('NSETS_',ps)]] -1)))
  colnames(stable.stats.df) <- c('Max Slice-Slice F_ST', 
                                 'Max Slice-Slice AF Difference',
                                 'No. Slices Reporting Variant with Positive Probability')
  
  # Create URL to dbSNP
  top.url <- paste0('https://www.ncbi.nlm.nih.gov/snp/', rel.gene.ps[[paste0('TOP_',ps)]])
  stable.url <- paste0('https://www.ncbi.nlm.nih.gov/snp/', rel.gene.ps[[paste0('STAB_',ps)]])
  return(list(TOP_VS_STABLE = top.vs.stable.df,
              STABLE_STATS = stable.stats.df,
              TOP_DBSNP = top.url,
              STAB_DBSNP = stable.url))
}

#' getFuncAnnotTable
#' 
#' This function takes in user-provided chromosome, gene name, potential set, 
#' and functional annotation choice, then proceeds with generating 
#' a table summarizing the score for the stable vs top variant.
#'
#' Example data for testing functionality:
#' chr <- 22
#' gene.name <- gene_names[[22]][1]
#' ps <- 1
#' func.class <- 'FAVOR'
#' func.annot <- 'CADD.RawScore'
getFuncAnnotTable <- function(chr, input.gene.name, ps, func.class, func.annot) {
  # Get functional annotations df for both stable and top variants
  top.annots.df <- readr::read_csv(paste0(annot_dir, 
                                             chr, '/top_annots.csv')) 
  stab.annots.df <- readr::read_csv(paste0(annot_dir, 
                                          chr, '/stab_annots.csv')) 
  
  # Convert input gene name to Ensembl gene name
  gene.name <- (gene_dfs[[chr]] %>% 
                  subset(INPUT == input.gene.name))[['ensembl_version_id']]
  
  # Subset to relevant gene, ps and func.annot
  if (func.class == 'Enformer') {
    top.value.ave <- (top.annots.df %>% 
                        subset(GENE == gene.name))[[paste0(gsub('.{1}$','',stringr::str_split(func.annot, '\\(')[[1]][2]),'_AVE_',ps)]]
    stab.value.ave <- (stab.annots.df %>% 
                         subset(GENE == gene.name))[[paste0(gsub('.{1}$','',stringr::str_split(func.annot, '\\(')[[1]][2]),'_AVE_',ps)]]
    top.value.tss <- (top.annots.df %>% 
                        subset(GENE == gene.name))[[paste0(gsub('.{1}$','',stringr::str_split(func.annot, '\\(')[[1]][2]),'_TSS_',ps)]]
    stab.value.tss <- (stab.annots.df %>% 
                         subset(GENE == gene.name))[[paste0(gsub('.{1}$','',stringr::str_split(func.annot, '\\(')[[1]][2]),'_TSS_',ps)]]
    
    to.return <- data.frame(TOP = abs(c(top.value.tss,top.value.ave)),
                            STABLE = abs(c(stab.value.tss,stab.value.ave)),
                            row.names = c(paste0(func.annot,' , TSS-centered'), paste0(func.annot,' , Averaged')))
    colnames(to.return) <- c('Top Variant', 'Stable Variant')
  } else {
    top.value <- (top.annots.df %>% subset(GENE == gene.name))[[paste0(func.annot,'_',ps)]]
    stab.value <- (stab.annots.df %>% subset(GENE == gene.name))[[paste0(func.annot,'_',ps)]]
    
    to.return <- data.frame(TOP = top.value,
                            STABLE = stab.value,
                            row.names = func.annot)
    colnames(to.return) <- c('Top Variant', 'Stable Variant')
  }
  
  # Return
  return(to.return)
}

#' getDensityPlot
#' 
#' This function takes in user-provided functional annotation class,
#' functional annotation, and potential set, and generates 
#' the density plots for the matching vs non-matching (stable / top) variants.
#'
getDensityPlot <- function(func.class, func.annot, ps) {
  func.class.plot.list <- readRDS(file = paste0(annot_dir,func.class,'_density_plots.rds'))
  # print('hello')
  # print(func.class)
  # print(func.annot)
  # print(ps)
  # print(class(func.class.plot.list[[func.annot]][[ps]]))
  if (func.class == 'Enformer') {
    to.return <- func.class.plot.list[[gsub('.{1}$','',stringr::str_split(func.annot, '\\(')[[1]][2])]][[as.numeric(ps)]]
  } else {
    to.return <- func.class.plot.list[[func.annot]][[as.numeric(ps)]]
  }
  return(plot(to.return))
}

#' getPairwisePlot
#' 
#' This function takes in user-provided functional annotation class,
#' functional annotation, and potential set, and generates 
#' the pairwise plot for the stable vs top variant.
#'
getPairwisePlot <- function(func.class, func.annot, ps) {
  func.class.plot.list <- readRDS(file = paste0(annot_dir,func.class,'_plots.rds'))
  # print('hello')
  # print(func.class)
  # print(func.annot)
  # print(ps)
  # print(class(func.class.plot.list[[func.annot]][[ps]]))
  if (func.class == 'Enformer') {
    to.return <- func.class.plot.list[[gsub('.{1}$','',stringr::str_split(func.annot, '\\(')[[1]][2])]][[as.numeric(ps)]]
  } else {
    to.return <- func.class.plot.list[[func.annot]][[as.numeric(ps)]]
  }
  # Handle ggplot2 version incompatibility by rebuilding plot from data
  tryCatch({
    return(plot(to.return))
  }, error = function(e) {
    if (grepl("continuous_range", e$message)) {
      # Extract data from the old ggplot object and rebuild with original layers:
      # 1. geom_point(alpha = 0.02)
      # 2. geom_bin2d(bins = 100)
      # 3. geom_abline(intercept = 0, slope = 1, linetype = "dashed")
      # 4. coord_fixed()
      plot_data <- to.return$data
      plot_labels <- to.return$labels

      p <- ggplot(plot_data, aes(x = TOP, y = STABLE)) +
        geom_point(alpha = 0.02) +
        geom_bin2d(bins = 100) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        coord_fixed() +
        labs(x = plot_labels$x, y = plot_labels$y) +
        theme_bw()
      return(plot(p))
    } else {
      stop(e)
    }
  })
}
