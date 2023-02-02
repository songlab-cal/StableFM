################################################
#### Stand-alone script for executing PICS2 ####
################################################

######################
##                  ##
## Helper functions ##
##                  ##
######################

#' getLeadSNPs
#' 
#' Computes the leading column of a \eqn{N\times P} covariate matrix
#' when regressing against the \eqn{N\times 1} response variable, \eqn{y}.
#' Returns a list with three elements: (1) All the marginal p-values,
#' (2) Leading SNP, (3) All SNPs with LD (\eqn{r^2}) > ld_thres to leading SNP
#' 
#' @param X The matrix on which column-wise regression is performed
#' @param y The response variable used in the regression
#' @param ld_thres How much LD to the lead SNP should we use 
#' as criteria to include other SNPs? 
#' 
#' @return A list containing three elements: log10(p-values), 
#' lead indices, and the set of SNPs with LD > ld_thres to the lead SNP.
#' 
getLeadSNPs <- function(X, 
                        y,
                        ld_thres) {
  # If column names are NULL, give names to columns(X)
  if (is.null(colnames(X))) {
    message(date(), ": Column of X unnamed, proceeding to give names rs1, rs2, ...")
    colnames(X) <- paste0("rs", 1:ncol(X))
  }
  
  # Define some local variables used in univariate regression
  S_y <- sum(y)
  S_yy = sum(y^2)
  
  # Local function to compute p-value for a column
  getMarginalP <- function(x) {
    S_x <- sum(x)
    S_xx <- sum(x^2)
    S_xy <- sum(x * y)
    n <- length(x)
    beta_hat <- (n * S_xy - S_x * S_y) / (n * S_xx - S_x^2)
    
    s_epsilon2 <- (n * S_yy - S_y^2 - beta_hat^2 * (n * S_xx - S_x^2)) / 
      (n * (n - 2))
    s_beta_hat2 <- n * s_epsilon2 / 
      (n * S_xx - S_x^2)
    
    t_score <- beta_hat / sqrt(s_beta_hat2)
    
    if (n > 50) {
      return(-log10(pchisq(t_score^2, df = 1, lower.tail = FALSE)))
    } else {
      return(-log10(pt(abs(t_score), df = n - 2, lower.tail = FALSE)))
    }
  }
  
  # Compute -log10(p) for each column
  pvals <- apply(X, 2, getMarginalP)
  
  # Find the lead SNPs
  lead_indices <- colnames(X)[which(pvals == max(pvals, na.rm = TRUE))]
  lead_indices <- max(lead_indices) # pick the maximum 
  
  # Find all SNPs with r^2 > ld_thres to the lead SNP 
  neighbor_set <- c()
  for (j in colnames(X)) {
    if (var(X[,j]) == 0) {
      message(date(), ": No polymorphism detected at ", j)
    } else {
      if (any(cor(X[,j], X[,lead_indices])^2 > ld_thres^2)) {
        neighbor_set <- c(neighbor_set, j)
      }
    }
  }
  
  # Return
  return(list(log10P = pvals,
              lead_indices = lead_indices,
              lead_set = neighbor_set))
}

#' restrictShuffle
#' 
#' For a given triplet of matrix \eqn{X}, vector of responses \eqn{y}, and 
#' integer \eqn{i}, returns a resampled vector of responses, \eqn{y^{\ast}}. 
#' This is done by fixing a focal column i of the original matrix \eqn{X}. 
#' Assuming the focal column is discrete, a random permutation of the indices is
#' performed such that only indices in the focal column sharing the same level 
#' can be permuted. For example, for the focal column \eqn{(1,1,1,2,2,0,0,1)}, 
#' any random permutation can only permute the zeros amongst themselves, 
#' the ones amongst themselves, and so on. The resulting random permutation is 
#' then applied to the vector \eqn{y} to obtain \eqn{y^{\ast}}.
#' 
#' @param X The matrix a resampled version of which is to be produced
#' @param y The vector of responses
#' @param i The column to fix 
#' 
#' @return The resampled vector, \eqn{y^{\ast}}
#' 
restrictShuffle <- function(X,
                            y,
                            i) {
  # Define permutation
  sigma <- rep(0,dim(X)[1])
  
  # Pull out focal column
  x_i <- as.factor(X[,i])
  # (Optional) If x_i has k levels, with sizes (n_1,...,n_k)
  # for the levels and n_1!...n_k! not too large, consider
  # generating all possible permutations
  
  # Generate a random permutation
  for (level in levels(x_i)) {
    level_inds <- which(x_i == level)
    if (length(level_inds) == 1) {
      sigma[level_inds] <- level_inds
    } else {
      sigma[level_inds] <- sample(x = level_inds, 
                                  size = length(level_inds), 
                                  replace = FALSE)
    }
  }
  
  # Return
  return(y[sigma])
}

#' PICS Base Function
#' 
#' Computes potential sets by running PICS, a Bayesian non-parametric
#' approach to inferring causal variants for fine-mapping. See the 
#' manuscript for details on this method. Note that any necessary pre-processing 
#' of response variable, y, should be performed. Returns a list of length L,
#' with each element itself a list containing the posterior probabilities of
#' variants and the variant with highest posterior probability.  
#' 
#' Dependencies
#' @param X The allele dosage matrix (\eqn{N \times P})
#' @param y The phenotype vector (case-control or quantitative) (\eqn{N\times 1})
#' @param pr_prob A vector of prior probabilities of each variant being causal. 
#' Default is \eqn{1/P}
#' @param resample_n The resampling number used to estimate \eqn{P(A\text{ lead}|B\text{ causal})}
#' @param ld_thres The LD threshold used to include SNPs in LD with the lead SNP 
#' @param L The number of causal variants / credible sets to detect 
#' 
#' @return A list of lists. The higher list has length L. Within each element is 
#' a lower list of (1) posterior probabilities of each SNP being causal; and 
#' (2) index of SNPs with maximum posterior probability of being causal
#' 
picsBase <- function(X, 
                     y, 
                     pr_prob = rep(1/dim(X)[2], dim(X)[2]), 
                     resample_n = 500,
                     ld_thres = 0.5,
                     L = 1) {
  # If X has no column names, add them
  if (is.null(colnames(X))) {
    message(date(), ": Column of X unnamed, proceeding to give standard names rs1, rs2, ...")
    colnames(X) <- paste0("rs", 1:ncol(X))
  }
  
  # If prior probabilities have no names, 
  # allot them names from X
  if (is.null(names(pr_prob))) {
    message(date(), ": Prior probabilities are unnamed, giving matching names to SNPs")
    names(pr_prob) <- colnames(X)
  }
  
  # Find the lead SNP in original data
  set_list <- list()
  ell <- 0
  message(date(), paste0(": Checking dimension of genotype matrix -- ", dim(X)[1], " by ", dim(X)[2]))
  X_ <- X # make local copy to modify 
  while (ell < L) {
    message(date(), ": Computing lead SNPs and lead set ", ell + 1, "...")
    comp_lead_snp <- getLeadSNPs(X_, y, ld_thres)
    lead_snp <- max(comp_lead_snp$lead_indices, na.rm = TRUE) 
    message(date(), ": Lead SNP found is ", lead_snp, ".")
    lead_set <- comp_lead_snp$lead_set
    message(date(), ": Found ", length(lead_set) - 1, " SNPs with |r| > ", ld_thres, " to the lead SNP.")
    
    # Save lead_snp and lead_set 
    set_list[[ell + 1]] <- list(LEAD_SNP = lead_snp,
                                LEAD_SET = lead_set)
    # Reset X 
    X_ <- X_[,!(colnames(X_) %in% lead_set)]
    
    # Increment ell
    ell <- ell + 1
  }
  
  # Make post_prob_list
  post_prob_list <- list()
  
  # Construct P(A_lead | B) vector 
  for (ell in 1:L) {
    A_given_B <- c()
    
    # For each variant
    if (length(set_list[[ell]][[2]]) == 1) {
      post_prob_list[[ell]] <- list(post_prob = 1,
                                    max_index = set_list[[ell]][[2]])
      # Name of single causal variant is max index
      names(post_prob_list[[ell]]$post_prob) <- post_prob_list[[ell]]$max_index
      
    } else {
      for (i in set_list[[ell]][[2]]) {
        # Compute P(A_lead | B_i causal)
        repeat_leads <- foreach(j=1:resample_n, .combine = c) %dopar% 
          (set_list[[ell]][[1]] %in% getLeadSNPs(X[,set_list[[ell]][[2]]], 
                                                 restrictShuffle(X[,set_list[[ell]][[2]]],y,i), 
                                                 ld_thres = 0)$lead_indices) 
        A_given_B <- c(A_given_B,mean(repeat_leads))
      }
      
      # Give names to conditional probabilities 
      names(A_given_B) <- set_list[[ell]][[2]]
      
      # Compute posterior probability under conditionally exchangeable null
      post_prob <- A_given_B * pr_prob[set_list[[ell]][[2]]] / 
        sum(A_given_B * pr_prob[set_list[[ell]][[2]]])
      
      post_prob_list[[ell]] <- list(post_prob = post_prob,
                                    max_index = names(post_prob)[which(post_prob == max(post_prob))])
    } 
  }
  
  # Return
  return(post_prob_list)
}

#' oneHotEncode
#' 
#' Given an array containing external categorical slice labels S, 
#' returns a new array that one-hot encodes each slice. An additional column is
#' generated to represent the full sample.
#' 
#' Dependencies: one_hot
#' @param S The \eqn{K} external slice labels (\eqn{N \times K}). Columns need 
#' not be binary.
#' 
#' @return A \eqn{N \times K'} array, which is a one-hot encoded version of the 
#' original slice label array. The first column captures inclusion of all 
#' individuals  
#' 
oneHotEncode <- function(S) {
  # One-hot encode if column not already binary
  S_col <- 0
  for (i in 1:ncol(S)){
    n_unique_vals <- length(unique(S[,i]))
    if (n_unique_vals == 2){ # one col needed for binary variable
      S_col <- S_col + 1
    } else {
      S_col <- S_col + n_unique_vals # unique_val cols needed under one-hot encoding
    }
  }
  
  message(date(), ": Found ", S_col, 
          " one-hot-encoded external annotations from ", 
          ncol(S), " original annotations.")
  
  # Create updated matrix with one-hot encoding
  S_onehot <- matrix(NA, nrow = nrow(S), ncol = S_col) 
  S_colnames <- colnames(S)
  S_onehot_colnames <- c()
  
  current_index <- 1
  for (i in 1:ncol(S)){
    current_col <- S[,i]
    unique_vals <- sort(unique(current_col))
    n_unique_vals <- length(unique_vals)
    if (n_unique_vals == 2){ # force values to be 0 or 1
      current_col[current_col == unique_vals[1]] <- 0
      current_col[current_col == unique_vals[2]] <- 1
      
      S_onehot[,current_index] <- current_col
      S_onehot_colnames <- c(S_onehot_colnames, S_colnames[i])
      current_index <- current_index + 1
    } else{
      new_cols <- as.matrix(
        mltools::one_hot(as.data.table(factor(sapply(current_col, as.character)))))
      S_onehot[,seq(current_index, current_index + n_unique_vals - 1)] <- new_cols
      S_onehot_colnames <- c(S_onehot_colnames, paste0(S_colnames[i], '_', 0:(n_unique_vals-1)))
      current_index <- current_index + n_unique_vals
    }
  }
  colnames(S_onehot) <- S_onehot_colnames
  
  S <- S_onehot # Use one-hot encoded S for downstream steps
  to_return <- cbind(rep(1,nrow(S)), S_onehot)
  colnames(to_return) <- c('ALL', S_onehot_colnames)
  
  # Return
  return(to_return)
}

#' summarizeStabResults
#' 
#' Given a list containing posterior probabilities of running PICS on multiple
#' slices and for L potential sets, returns a dataframe that summarizes the 
#' stable variants for each potential set. The posterior probability of the
#' stable variant in the ALL slice is also provided, as is the number of variants
#' in the ALL slice appearing with positive probability after running PICS2. 
#' 
#' @param pics_result An output list generated by picsStable with summarize=FALSE
#' 
#' @return A \eqn{L\times 3} dataframe with each row reporting the stable variant
#' for the corresponding potential set, alongside some fine-mapping characteristics. 
#' 
summarizeStabResults <- function(pics_result) {
  # Get list of all variants appearing in at least one slice
  n_annots <- length(pics_result)
  var_list <- c()
  for (s in 1:n_annots) {
    # edge case: NA (because only one SNP is available for fine-mapping)
    if (all(is.na(pics_result[[s]]))) {
      # Return 0
      return(-99)
    }
    L <- length(pics_result[[s]]) # number of potential sets
    for (l in 1:L) {
      var_list <- c(var_list, names(pics_result[[s]][[l]]$post_prob)) 
      var_list <- c(var_list, pics_result[[s]][[l]]$max_index)
    }
  }
  var_list <- unique(var_list) 
  
  # Get list of posterior probabilities for each slice,
  # as well as the potential sets containing that variant
  # for a slice
  summary_df <- data.frame(marker.ID = var_list)
  set_labels <- data.frame(marker.ID = var_list)
  for (s in 1:n_annots) {
    post_probs <- rep(0, length(var_list))
    set_label <- rep(0, length(var_list))
    names(post_probs) <- summary_df$marker.ID
    names(set_label) <- summary_df$marker.ID
    L <- length(pics_result[[s]])
    for (l in 1:L) {
      if (all(is.na(pics_result[[s]][[l]]))) {
        post_probs[pics_result[[s]][[l]]$max_index] <- NA
        set_label[pics_result[[s]][[l]]$max_index] <- l
      } else if (length(pics_result[[s]][[l]]$post_prob) == 1) {
        post_probs[pics_result[[s]][[l]]$max_index] <- 1
        set_label[pics_result[[s]][[l]]$max_index] <- l
      } else {
        post_probs[names(pics_result[[s]][[l]]$post_prob)] <- pics_result[[s]][[l]]$post_prob
        set_label[names(pics_result[[s]][[l]]$post_prob)] <- l
      }
    }
    summary_df <- cbind(summary_df, post_probs)
    set_labels <- cbind(set_labels, set_label)
  }
  
  names(summary_df)[-c(1)] <- names(pics_result)
  names(set_labels)[-c(1)] <- paste0(names(pics_result), "_set")
  
  # Generate table of all SNPs with positive posterior probability and slice data
  # [!] OPTIONAL - return gene_summary_df can be useful
  gene_summary_df <- summary_df %>% 
    mutate(TOTAL = apply(summary_df[,-1], 1, function(x) {sum(x != 0)})) %>%
    mutate(TOTAL_SET = apply(summary_df[,-1], 1, function(x) { paste(names(which(x != 0)), collapse = ", ")}))
  
  gene_summary_df <- merge(gene_summary_df, set_labels, by = "marker.ID")
  
  # Transform table into summary dataframe for returning
  L <- length(pics_result[['ALL']])
  stab_variant_vec <- c()
  post_prob_vec <- c()
  support_vec <- c()
  for (l in 1:L) {
    stab_variant_vec <- c(stab_variant_vec,
                          (gene_summary_df %>% 
                             filter(ALL_set == l) %>% 
                             filter(TOTAL == max(TOTAL)) %>% 
                             filter(ALL == max(ALL)))$marker.ID[1])
    post_prob_vec <- c(post_prob_vec,
                       (gene_summary_df %>% 
                          filter(ALL_set == l) %>% 
                          filter(TOTAL == max(TOTAL)) %>% 
                          filter(ALL == max(ALL)))$ALL[1])
    support_vec <- c(support_vec,
                     length(pics_result[['ALL']][[l]]$post_prob))
  }
  to_return <- data.frame(VARIANT = stab_variant_vec,
                          POSTERIOR_PROB = post_prob_vec,
                          SUPPORT_SIZE = support_vec)
  rownames(to_return) <- paste0('Potential_Set_', 1:L)
  
  # Return
  return(to_return)
}

#' PICS Stability-Guided Approach
#'
#' Computes potential sets by running PICS on multiple slices of the same
#' training sample, where each slice is determined by an exogenous label.
#' For example, an exogenous label might be a subpopulation membership,
#' a superpopulation membership, or any other label apart from the case-control
#' label.   
#'
#' Takes in covariates, phenotypes, and slice labels and then 
#' returns a list of potential sets generated by running PICS on each slice. 
#' 
#' Example: picsStable(X=genotype_data_,y=y,S=S_test,pr_prob=rep(1/dim(X)[2], dim(X)[2]),resample_n=1000,ld_thres=0.5,L=3)
#' 
#' Dependencies: picsBase
#' @param X The binary or real matrix (\eqn{N \times P})
#' @param y The phenotype (case-control or quantitative) on which to run (\eqn{N \times 1})
#' @param pr_prob A vector of prior probabilities of each variant being causal. 
#' Default is \eqn{1/P} where \eqn{P} is the number of variants in \eqn{X}
#' @param S The \eqn{K} external slice labels (\eqn{N \times K}).
#' Non-binary categorical slice labels must be binarized
#' @param resample_n The resampling number used to estimate \eqn{P(A\text{ lead}|B\text{ causal})}
#' @param ld_thres The LD threshold used to include SNPs in LD with the lead SNP 
#' @param L The number of causal variants / credible sets to detect 
#' @param summarize Whether to convert the list to just a table of SNPs and posterior probabilities
#' 
#' @return A list of potential sets and posterior probabilities (if summarize = FALSE),
#' or a dataframe summarizing the top variants (if summarize = TRUE)
#' @export
#' @importFrom doParallel %dopar%
#'
picsStable <- function(X, 
                       y, 
                       pr_prob=rep(1/dim(X)[2], dim(X)[2]),
                       S, 
                       resample_n,
                       ld_thres,
                       L,
                       summarize) {
  # Give S column names if they don't exist
  if (is.null(colnames(S))) {
    colnames(S) <- paste0("S", 1:ncol(S))
  }
  
  one_hot_enc_slices <- oneHotEncode(S)
  
  all_post_probs <- list()
  
  # Proceed to perform slice-specific fine-mapping
  for (annot in colnames(one_hot_enc_slices)) {
    # Include subsample belonging to the slice
    X_slice <- X[which(one_hot_enc_slices[,annot] == 1),]
    y_slice <- y[which(one_hot_enc_slices[,annot] == 1)]
    
    # If phenotype has 0 variation for the slice, output NA for that annotation
    if (var(y_slice) == 0) {
      message(date(), ": Phenotype has zero variation for annotation = ", annot, "...")
      all_post_probs[[annot]] <- NA 
    } else if (is.null(dim(X_slice))) {
      # If genotype matrix consists of only 1 variant, return that variant
      message(date(), ": Only one variant is present")
      all_post_probs[[annot]] <- list(post_prob_list = NA,
                                      boot_df = NA,
                                      jackknife_est = NA) 
    } else {
      # Discard non-polymorphic variants
      message(date(), ": Removing non-polymorphic SNPs from smaller genotype matrix for annotation = ", annot, "...")
      poly_snvs <- which(apply(X_slice, 2, var) != 0)
      
      genotype_data_1 <- X_slice[,poly_snvs]
      snp_metadata_1 <- snp_metadata[poly_snvs,]
      #colnames(genotype_data_1) <- X_slice[,poly_snvs] # 10/14/21: ERROR HERE. Column names ALREADY given, no need to rename!
      
      genotype_data_1 <- scale(genotype_data_1)  
      y_slice <- scale(y_slice)
      
      # Run PICS2
      slice_fitted <- picsBase(X = genotype_data_1, 
                               y = y_slice,
                               pr_prob = pr_prob, 
                               resample_n = resample_n,
                               ld_thres = ld_thres,
                               L = L)
      all_post_probs[[annot]] <- slice_fitted
    }
  }
  
  # Return
  if (summarize) {
    to_return <- summarizeStabResults(all_post_probs)
    return(to_return)
  } else {
    return(all_post_probs)
  }
}


#' PICS Residualization Approach
#' 
#' Computes potential sets by residualizing phenotype by covariates that
#' capture potential confounders \eqn{Z} (e.g., genetic PCs that capture structure).
#' First residualizes the phenotype \eqn{y} against \eqn{Z} to obtain \eqn{y^{\ast}}, 
#' then runs PICS on the genotype and residualized phenotypes. 
#' 
#' Example: picsTop(X=genotype_data_,y=y,pr_prob=rep(1/dim(X)[2], dim(X)[2]),Z=Z,resample_n=1000,ld_thres=0.5,L=3)
#' 
#' Dependencies: picsBase
#' @param X The binary or real matrix (\eqn{N \times P})
#' @param y The phenotype (case-control or quantitative) on which to run (\eqn{N \times 1})
#' @param pr_prob A vector of prior probabilities of each variant being causal. 
#' Default is \eqn{1/P} where \eqn{P} is the number of variants in \eqn{X}
#' @param Z The \eqn{N\times K} matrix, where each column encodes a potential confounder
#' @param resample_n The resampling number used to estimate \eqn{P(A\text{ lead}|B\text{ causal})}
#' @param ld_thres The LD threshold used to include SNPs in LD with the lead SNP 
#' @param L The number of causal variants / credible sets to detect 
#' @param summarize Whether to convert the list to just a table of SNPs and posterior probabilities
#' 
#' @return A list of potential sets and posterior probabilities (if summarize = FALSE),
#' or a dataframe summarizing the top variants (if summarize = TRUE)
#' @export
#' @importFrom doParallel %dopar%
#' 
picsTop <- function(X, 
                    y, 
                    pr_prob=rep(1/dim(X)[2], dim(X)[2]),
                    Z,
                    resample_n,
                    ld_thres,
                    L,
                    summarize) {
  # Return error message if Z is null
  assertthat::assert_that(!is.null(Z),
                          msg = 'Please specify a vector or a matrix of potential confounders')
  
  # Regress against potential confounders
  num.confounders <- ifelse(is.null(dim(Z)), 1, dim(Z)[2])
  message(date(), paste0(": Residualize against ", 
                         num.confounders, 
                         " potential confounders"))
  
  lm.input.df <- cbind(data.frame(y=y),
                       Z)
  colnames(lm.input.df)
  y_residuals <- lm(y~.-1, data = lm.input.df)$residuals
  
  poly_snvs <- which(apply(X, 2, var) != 0)
  genotype_data_1 <- X[,poly_snvs]
  snp_metadata_1 <- snp_metadata[poly_snvs,]
  genotype_data_1 <- scale(genotype_data_1)
  y_residuals <- scale(y_residuals)

  # Run PICS2
  pics_output <- picsBase(X = genotype_data_1,
                          y = y_residuals,
                          pr_prob = pr_prob,
                          resample_n = resample_n,
                          ld_thres = ld_thres,
                          L = L)
  
  # Return
  if (summarize) {
    to_return <- data.frame(VARIANT = unlist(lapply(pics_output, "[[", 2)),
                            POSTERIOR_PROB = unlist(lapply(pics_output, function(x){max(x[['post_prob']])})),
                            SUPPORT_SIZE = unlist(lapply(pics_output, function(x){length(x[['post_prob']])})))
    rownames(to_return) <- paste0('Potential_Set_', 1:L)
    return(to_return)
  } else {
    return(pics_output)
  }
  
}

#' runPICS2
#' 
#' Either runs the sliced version,
#' the regular version, or both.
#' 
#' Example: runPICS2(X=genotype_data_,y=y,S=S_test,Z=Z,resample_n=1000,ld_thres=0.5,L=3)
#' #' Dependencies: picsStable, picsTop
#' @param X The binary or real matrix (\eqn{N \times P})
#' @param y The phenotype (case-control or quantitative) on which to run (\eqn{N \times 1})
#' @param pr_prob A vector of prior probabilities of each variant being causal. 
#' Default is \eqn{1/P} where \eqn{P} is the number of variants in \eqn{X}
#' @param S The \eqn{K} external slice labels (\eqn{N \times K}).
#' Non-binary categorical slice labels must be binarized
#' @param Z The \eqn{N\times K} matrix, where each column encodes a potential confounder
#' @param resample_n The resampling number used to estimate \eqn{P(A\text{ lead}|B\text{ causal})}
#' @param ld_thres The LD threshold used to include SNPs in LD with the lead SNP 
#' @param L The number of causal variants / credible sets to detect 
#' @param option The approach(es) on which to run PICS --- either 'stable', 'top',
#' or 'both'
#' @param summarize Whether to convert the list to just a table of SNPs and posterior probabilities
#' 
#' @return Fine-mapping variants --- either stable, top or both.
#' 
runPICS2 <- function(X, 
                     y, 
                     pr_prob=NULL,
                     S=data.frame('ALL' = rep(1,nrow(X))),
                     Z=NULL,
                     resample_n,
                     ld_thres,
                     L,
                     option = "both",
                     summarize = TRUE) {
  # Check for mismatches in input data dimensions
  assertthat::assert_that(nrow(X) == length(y),
                          msg = 'The number of individuals in phenotype and allele dosage matrix must match.')
  if (!is.null(Z)) {
    assertthat::assert_that((nrow(X) == nrow(Z) | nrow(X) == length(Z)),
                            msg = 'The number of individuals in confounder vector/matrix and allele dosage matrix must match.')
  }
  
  if (is.null(pr_prob)) {
    pr_prob<-rep(1/dim(X)[2], dim(X)[2])
  } else {
    assertthat::assert_that(ncol(X) == length(pr_prob),
                            msg = 'The number of variants in allele dosage matrix and prior probability vector must match.') 
  }
  
  # Check for correct types
  assertthat::assert_that(round(L) == L,
                          msg = 'L must be a whole number.')
  assertthat::assert_that(ld_thres < 1 & ld_thres >= 0,
                          msg = 'The ld_thres parameter must be less than 1 and at least 0.')
  assertthat::assert_that(resample_n >= 500,
                          msg = 'Please use resample_n at least 500 for reliable estimates of posterior probabilities.')
  assertthat::assert_that((option == 'both' | option == 'stable' | option == 'top'),
                          msg = 'Please only choose top, stable or both as the option.')
  
  # Perform fine-mapping
  if (option == 'stable') {
    to_return <- picsStable(X=X,y=y,pr_prob=pr_prob,S=S,
                         resample_n=resample_n,ld_thres=ld_thres,
                         L=L,summarize=summarize) 
  } else if (option == 'top') {
    to_return <- picsTop(X=X,y=y,pr_prob=pr_prob,Z=Z,
                         resample_n=resample_n,ld_thres=ld_thres,
                         L=L,summarize=summarize) 
  } else {
    message(date(), ': Running stability-guided approach')
    stable_to_return <- picsStable(X=X,y=y,pr_prob=pr_prob,S=S,
                                   resample_n=resample_n,ld_thres=ld_thres,
                                   L=L,summarize=summarize) 
    message(date(), ': Running residualization approach')
    top_to_return <- picsTop(X=X,y=y,pr_prob=pr_prob,Z=Z,
                         resample_n=resample_n,ld_thres=ld_thres,
                         L=L,summarize=summarize) 
    to_return <- list(STABLE = stable_to_return,
                      TOP = top_to_return)
  }
  
  # Return
  return(to_return)
}

######################
##                  ##
##       Main       ##
##                  ##
######################

library(dplyr)
library(doParallel)

## Inputs provided by user
# commandArgs picks up the variables you pass from the command line
## Ex: Rscript runPICS2.R 500 3 toy-data/genotype_array.csv toy-data/phenotype_vector.txt toy-data/snp_metadata.txt toy-data/sample_metadata.txt toy-data/gene_tss_metadata.txt toy-data/covars_array.csv toy-results/
args <- commandArgs(trailingOnly = TRUE)
n_perms <- args[1] # Ex: 500
user_L <- args[2] # Ex: 3
geno_file <- args[3] # Ex: "toy-data/genotype_array.csv"
pheno_file <- args[4] # Ex: "toy-data/phenotype_vector.txt"
snp_metadata_file <- args[5] # Ex: "toy-data/snp_metadata.txt"
sample_metadata_file <- args[6] # Ex: "toy-data/sample_metadata.txt"
tss_metadata_file <- args[7] # Ex: "toy-data/gene_tss_metadata.txt"
covar_file <- args[8] # Ex: "toy-data/covars_array.csv"
out_dir <- args[9] # Ex: "toy-results/"

file.cxn <- file("logs/runPICS2_LOGFILE.txt", open = "wt")
sink(file = file.cxn, type = "message")

# Perform fine-mapping 
geno.df <- readr::read_csv(geno_file) # genotypes
pheno.df <- readr::read_delim(pheno_file, delim = '\t') # phenotypes
snp.metadata.df <- readr::read_delim(snp_metadata_file, delim = '\t') # SNP IDs
samp.metadata.df <- readr::read_delim(sample_metadata_file, delim = '\t') # indiv IDs
tss.metadata.df <- readr::read_delim(tss_metadata_file, delim = '\t') # TSS metadata
covars.df <- readr::read_csv(covar_file) # genetic PCs

# Check that IDs match across all individual-level data
assertthat::assert_that(identical(geno.df$ID,pheno.df$ID),
                        msg = "Individual IDs are not matched between genotype array and phenotype vector")
assertthat::assert_that(identical(geno.df$ID,samp.metadata.df$ID),
                        msg = "Individual IDs are not matched between genotype array and sample metadata")
assertthat::assert_that(identical(geno.df$ID,covars.df$ID),
                        msg = "Individual IDs are not matched between genotype array and covariates array")

# Check that genotype matrix has same number of variants as variants in SNP metadata file
assertthat::assert_that((ncol(geno.df)-1) == nrow(snp.metadata.df),
                        msg = "No. variants provided in genotype array does not match number of variants in metadata file")

# Extract relevant objects
ids <- geno.df[,1]
genos <- geno.df[,-1]
poly.snvs <- which(apply(genos, 2, var) != 0)
genotype.data <- genos[,poly.snvs]
snp_metadata <- snp.metadata.df[poly.snvs,]
colnames(genotype.data) <- snp_metadata$rsid  

# Filter variants to include only those lying in 1Mb of gene TSS
gene.name <- colnames(pheno.df)[2]
gene.tss <- tss.metadata.df %>% subset(GEUVADIS == gene.name)

message(date(), paste0(": GEUVADIS ID = ", gene.tss$GEUVADIS))
message(date(), paste0(": TSS = ", gene.tss$TSS))
message(date(), paste0(": LEFT ENDPOINT = ", gene.tss$LEFT))
message(date(), paste0(": RIGHT ENDPOINT = ", gene.tss$RIGHT))

included.snps <- snp_metadata %>% subset(pos >= gene.tss$LEFT & pos <= gene.tss$RIGHT)

X.input <- as.data.frame(genotype.data) %>% select(included.snps$rsid)
S.input <- as.data.frame(samp.metadata.df[,-1])
Z.input <- covars.df[,-1]
y.input <- pheno.df[[gene.name]]

user_result <- runPICS2(X=X.input, 
                        y=y.input, 
                        pr_prob=NULL,
                        S=S.input,
                        Z=Z.input,
                        resample_n=as.numeric(n_perms),
                        ld_thres=0.5,
                        L=as.numeric(user_L),
                        option = "both",
                        summarize = TRUE)

save(user_result, file = paste0(out_dir, 'test_result.rds'))
sink()

