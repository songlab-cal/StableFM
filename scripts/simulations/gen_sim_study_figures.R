################################################################################
########################## Simulation Study Figures ############################
################################################################################
#' 
#' Created: 2/14/2025
#' 
#' This script generates the following figures for the simulation study. 
#' 
#' Main Text:
#' A. Comparison of Plain PICS and Stability-guided PICS on environmental heterogeneity simulations.
#' B. SNR vs Recovery plot (Combined PICS vs Stable PICS vs Top PICS)
#' C. SNR vs Recovery plot (Top SuSiE vs Stable SuSiE)
#' D. Matching versus Non-matching causal variant recovery rate (PICS and SuSiE).
#' 
#' Supplement:
#' 1. Comparison of Plain PICS and Stability-guided PICS (PS2 and PS3)
#' 2. PICS SNR vs Recovery plot (PS2 and PS3)
#' 3. SuSiE SNR vs Recovery plot (PS2 and PS3)
#' 4. Matching vs Non-matching causal variant recovery rate (PS2 and PS3)
#' 5. Impact of including more potential sets on matching frequency
#' 6. Difference in posterior probability between non-matching variants
#' 7. SuSiE vs PICS
#' 8. Split results by no. causal variants and SNR-stratified analysis

#### Load directories and packages ---------------------------------------------
## Define directories (large simulation study)
pics2_res_dir <- "~/Documents/het_and_strat/082924/"
pics2_top_res_dir <- "~/Documents/het_and_strat/092324/"
pics2_plain_res_dir <- "~/Documents/het_and_strat/092624/"
pics2_combined_res_dir <- "~/Documents/het_and_strat/092724/"
susie_res_dir <- "~/Documents/het_and_strat/092424/"
stable_susie_res_dir <- "~/Documents/het_and_strat/092924/"

## Define directories (small env het simulation study)
t8_het_env_res_dir <- "~/Documents/het_and_strat/012525/"
t16_het_env_res_dir <- "~/Documents/het_and_strat/012725/"
t128_het_env_res_dir <- "~/Documents/het_and_strat/012825/"
t256_het_env_res_dir <- "~/Documents/het_and_strat/012925/"
im3_het_env_res_dir <- "~/Documents/het_and_strat/013025/"
ie3_het_env_res_dir <- "~/Documents/het_and_strat/013125/"

## Define directory for saving plots
plot_dir <- "~/Documents/het_and_strat/sim_study_figures/"

## Load libraries
library(dplyr)
library(ggplot2)

#### Helper functions ----------------------------------------------------------
## doesSetContainCausalSNP =====================================================
#' Compares a set of inferred causal SNPS with the true causal SNP to check
#' if the latter is contained in the former
doesSetContainCausalSNP <- function(causal_snp,snp_set) {
  return(grepl(causal_snp,snp_set)+0)
}

## getIntersect ================================================================
#' Takes two strings, splits each by comma-delimiter, then checks if the two 
#' sets of substrings contain common elements
getIntersect <- function(string_1, string_2) {
  v1 <- unlist(strsplit(string_1, split = ","))
  v2 <- unlist(strsplit(string_2, split = ","))
  common_elements <- intersect(v1, v2)
  
  # Return TRUE if there's at least one common element, otherwise FALSE
  return(length(common_elements) > 0)
}

## getIntersect2 ===============================================================
#' Takes two strings, splits each by comma-delimiter, then checks if the two 
#' sets of substrings contain common elements. If true, returns the common 
#' elements as a single string with elements separated by comma. If false, 
#' returns NA. This function outputs strings to be passed into getIntersect 
#' by design. 
getIntersect2 <- function(string_1,string_2) {
  v1 <- unlist(strsplit(string_1, split = ","))
  v2 <- unlist(strsplit(string_2, split = ","))
  common_elements <- intersect(v1, v2)
  
  if (length(common_elements) > 0) {
    return(paste(common_elements,sep=","))
  } else {
    return(NA)
  }
}

## getIntersect3 ===============================================================
#' Takes three strings, splits each by comma-delimiter, then checks if the three 
#' sets of substrings contain common elements. If true, returns the common 
#' elements as a single string with elements separated by comma. If false, 
#' returns NA. This function outputs strings to be passed into getIntersect 
#' by design. 
getIntersect3 <- function(string_1,string_2,string_3) {
  v1 <- unlist(strsplit(string_1, split = ","))
  v2 <- unlist(strsplit(string_2, split = ","))
  v3 <- unlist(strsplit(string_3, split = ","))
  common_elements <- intersect(intersect(v1, v2),v3)
  
  if (length(common_elements) > 0) {
    return(paste(common_elements,sep=","))
  } else {
    return(NA)
  }
}

## getSNRIntersect =============================================================
#' Takes two dataframe of fine-mapping results, stratifies by simulation SNR,
#' and tabulates agreement fraction between the two approaches. Uses getIntersect.
#' getSNRIntersect(top_pics2_one_causal_results,stable_pics2_one_causal_results
#'                 snp=1)
getSNRIntersect <- function(df_1,df_2,snp) {
  df_to_return <- NULL
  for (phi in c(0.05,0.1,0.2,0.4)) {
    snr_df_1 <- df_1 %>% subset(Phi==phi)
    snr_df_2 <- df_2 %>% subset(Phi==phi)
    true_false_vec <- mapply(getIntersect,
                             snr_df_1[[paste0("SNP",snp)]], 
                             snr_df_2[[paste0("SNP",snp)]]) 
    falses <- length(true_false_vec) - sum(true_false_vec)
    trues <- sum(true_false_vec)
    message("SNR = ", phi/(1-phi))
    message("Falses = ", falses)
    message("Trues = ", trues)
    df_to_return <- rbind(df_to_return, 
                          data.frame(SNR=phi/(1-phi),
                                     Match=trues,
                                     NonMatch=falses))
    
  }
  df_to_return$MatchPct <- 100*df_to_return$Match/
    (df_to_return$Match+df_to_return$NonMatch)
  return(df_to_return)
}

## countIntersect ==============================================================
#' Takes two strings, splits each by comma-delimiter, then counts the number of
#' elements in common.
countIntersect <- function(string_1, string_2) {
  v1 <- unlist(strsplit(string_1, split = ","))
  v2 <- unlist(strsplit(string_2, split = ","))
  common_elements <- intersect(v1, v2)
  
  # Return the number of common elements
  return(length(common_elements))
}
## countCausalSNPs =============================================================
#' Compares a set of true causal SNPS with the inferred causal SNPs across three
#' potential sets, and tabulates the number of matching SNPs. Uses getIntersect2
#' and countIntersect.
#' 
#' Test case 1 
#' countCausalSNPs(all_stable_susie_results$CausalSNPs[42],all_stable_susie_results$SNP1[42],all_stable_susie_results$SNP2[42],all_stable_susie_results$SNP3[42])
#' countCausalSNPs(all_stable_susie_results$CausalSNPs[2],all_stable_susie_results$SNP1[2],all_stable_susie_results$SNP2[2],all_stable_susie_results$SNP3[2])
countCausalSNPs <- function(causal_string,
                            snp1,
                            snp2,
                            snp3) {
  
  # Get comma-delimited string of causal SNPs intersecting each set's SNPs
  if (identical(snp1,NA)) {
    snp1_intersect <- NA
  } else {
    snp1_intersect <- getIntersect2(snp1,causal_string)
  }
  if (identical(snp2,NA)) {
    snp2_intersect <- NA
  } else {
    snp2_intersect <- getIntersect2(snp2,causal_string)
  }
  if (identical(snp3,NA)) {
    snp3_intersect <- NA
  } else {
    snp3_intersect <- getIntersect2(snp3,causal_string)
  }
  all_intersect <- union(union(snp1_intersect,snp2_intersect),snp3_intersect)
  
  # Count the number of SNPs recovered
  # NULL means snp1, snp2 and snp3 all have NAs
  # NA means at least one of snp1,snp2 or snp3 is not NA, but no matches found
  causal_snp_vec <- unlist(strsplit(causal_string, split = ","))
  all_count <- length(intersect(all_intersect,causal_snp_vec))
  # if (identical(all_intersect,NA)) {
  #   all_count <- 0
  # } else {
  #   all_count <- length(intersect(all_intersect,causal_snp_vec))
  # }
  return(all_count)
}

## getMatchRecProbDF ===========================================================
#' [!] DO NOT USE THIS FUNCTION [!]
#' Takes two dataframe of fine-mapping results (top vs stable), and stratifies
#' by whether the top and stable variants MATCH. Computes the recovery probability
#' for three cases: matching variants, non-matching top variant and non-matching
#' stable variant. Uses getIntersect. By default, will split the analysis by 
#' SNR used in simulations. If SNR is set to `FALSE`, it will simply compute
#' the statistics across all SNR scenarios in the provided datasets.
#' getMatchRecProbDF(top_pics2_one_causal_results,
#'                   stable_pics2_one_causal_results,
#'                   SNR=FALSE)
getMatchRecProbDF <- function(top_df,stable_df,SNR=TRUE) {
  df_intermediate <- NULL
  if (SNR) {
    for (phi in c(0.05,0.1,0.2,0.4)) {
      snr_df_1 <- top_df %>% subset(Phi==phi)
      snr_df_2 <- stable_df %>% subset(Phi==phi)
      true_false_vec <- mapply(getIntersect,
                               snr_df_1[[paste0("SNP1")]], 
                               snr_df_2[[paste0("SNP1")]]) 
      
      snr_df_1$Match <- true_false_vec
      snr_df_2$Match <- true_false_vec
      match_snr_df_2 <- snr_df_2 %>% subset(Match) # matching top and stable
      nonmatch_snr_df_1 <- snr_df_1 %>% subset(!Match) # non-matching top 
      nonmatch_snr_df_2 <- snr_df_2 %>% subset(!Match) # non-matching stable
      
      match_snr_df_2_contains_causal <- mapply(
        getIntersect,
        match_snr_df_2$CausalSNPs,match_snr_df_2$SNP1
      )  
      nonmatch_snr_df_2_contains_causal <- mapply(
        getIntersect,
        nonmatch_snr_df_2$CausalSNPs,nonmatch_snr_df_2$SNP1
      )  
      nonmatch_snr_df_1_contains_causal <- mapply(
        getIntersect,
        nonmatch_snr_df_1$CausalSNPs,nonmatch_snr_df_1$SNP1
      )  
      
      df_intermediate <- rbind(
        df_intermediate, 
        data.frame(
          Phi=phi,
          Matching=mean(match_snr_df_2_contains_causal),
          NonMatching_Top=mean(nonmatch_snr_df_1_contains_causal),
          NonMatching_Stable=mean(nonmatch_snr_df_2_contains_causal),
          N_matching=length(match_snr_df_2_contains_causal),
          N_nonMatching_top=length(nonmatch_snr_df_1_contains_causal),
          N_nonMatching_stable=length(nonmatch_snr_df_2_contains_causal)))
      
    }
    
    # Gather intermediate results into desired dataframe
    df_to_return <- NULL
    for (phi in c(0.05,0.1,0.2,0.4)) {
      rel_row <- df_intermediate %>% subset(Phi==phi)
      df_to_return <- 
        rbind(df_to_return, 
              data.frame(SNR=phi/(1-phi),
                         RecProb=rel_row[["Matching"]],
                         BinomSE=sqrt(rel_row[["Matching"]]*
                                        (1-rel_row[["Matching"]])/
                                        rel_row[["N_matching"]]),
                         Type="Matching"))
      df_to_return <- 
        rbind(df_to_return, 
              data.frame(SNR=phi/(1-phi),
                         RecProb=rel_row[["NonMatching_Top"]],
                         BinomSE=sqrt(rel_row[["NonMatching_Top"]]*
                                        (1-rel_row[["NonMatching_Top"]])/
                                        rel_row[["N_nonMatching_top"]]),
                         Type="Non-Matching Top"))
      df_to_return <- 
        rbind(df_to_return, 
              data.frame(SNR=phi/(1-phi),
                         RecProb=rel_row[["NonMatching_Stable"]],
                         BinomSE=sqrt(rel_row[["NonMatching_Stable"]]*
                                        (1-rel_row[["NonMatching_Stable"]])/
                                        rel_row[["N_nonMatching_stable"]]),
                         Type="Non-Matching Stable"))
    }
  } else {
    # Aggregate results across all SNRs
    true_false_vec <- mapply(getIntersect,
                             top_df[[paste0("SNP1")]], 
                             stable_df[[paste0("SNP1")]])
    
    top_df$Match <- true_false_vec
    stable_df$Match <- true_false_vec
    match_snr_df_2 <- stable_df %>% subset(Match) # matching top and stable
    nonmatch_snr_df_1 <- top_df %>% subset(!Match) # non-matching top 
    nonmatch_snr_df_2 <- stable_df %>% subset(!Match) # non-matching stable
    
    match_snr_df_2_contains_causal <- mapply(
      getIntersect,
      match_snr_df_2$CausalSNPs,match_snr_df_2$SNP1
    )  
    nonmatch_snr_df_2_contains_causal <- mapply(
      getIntersect,
      nonmatch_snr_df_2$CausalSNPs,nonmatch_snr_df_2$SNP1
    )  
    nonmatch_snr_df_1_contains_causal <- mapply(
      getIntersect,
      nonmatch_snr_df_1$CausalSNPs,nonmatch_snr_df_1$SNP1
    )  
    
    df_intermediate <- rbind(
      df_intermediate, 
      data.frame(Matching=mean(match_snr_df_2_contains_causal),
                 NonMatching_Top=mean(nonmatch_snr_df_1_contains_causal),
                 NonMatching_Stable=mean(nonmatch_snr_df_2_contains_causal),
                 N_matching=length(match_snr_df_2_contains_causal),
                 N_nonMatching_top=length(nonmatch_snr_df_1_contains_causal),
                 N_nonMatching_stable=length(nonmatch_snr_df_2_contains_causal)))
    # Gather intermediate results into desired dataframe
    df_to_return <- NULL
    df_to_return <- 
      rbind(df_to_return,
            data.frame(
              RecProb=df_intermediate[["Matching"]],
              BinomSE=sqrt(df_intermediate[["Matching"]]*
                             (1-df_intermediate[["Matching"]])/
                             df_intermediate[["N_matching"]]),
              Type="Matching"))
    df_to_return <- 
      rbind(df_to_return, 
            data.frame(
              RecProb=df_intermediate[["NonMatching_Top"]],
              BinomSE=sqrt(df_intermediate[["NonMatching_Top"]]*
                             (1-df_intermediate[["NonMatching_Top"]])/
                             df_intermediate[["N_nonMatching_top"]]),
              Type="Non-Matching Top"))
    df_to_return <- 
      rbind(df_to_return, 
            data.frame(
              RecProb=df_intermediate[["NonMatching_Stable"]],
              BinomSE=sqrt(df_intermediate[["NonMatching_Stable"]]*
                             (1-df_intermediate[["NonMatching_Stable"]])/
                             df_intermediate[["N_nonMatching_stable"]]),
              Type="Non-Matching Stable"))
    
  }
  # Return
  return(df_to_return)
}

## getMatchRecProbDF2 ==========================================================
#' Takes two dataframe of fine-mapping results (top vs stable), and stratifies
#' by whether the top and stable variants match. Computes the recovery probability
#' for three cases: matching variants, non-matching top variant and non-matching
#' stable variant. It uses getIntersect2 instead of getIntersect, unlike 
#' getMatchRecProbDF. This allows the matching variants to be explicitly identified,
#' since it is technically possible for stable and top results to have non-empty
#' intersection, but that intersection does not contain the causal variant while
#' one of the two results does (e.g., {A1,A3,A4} and {A2,A3,A4} intersect but this
#' intersection does not contain A1, which is the true causal).
#' By default, will split the analysis by SNR used in 
#' simulations. If SNR is set to `FALSE`, it will simply compute the statistics
#' across all SNR scenarios in the provided datasets. 
#' Uses getIntersect2 and getIntersect. 
#' getMatchRecProbDF2(top_pics2_one_causal_results,
#'                    stable_pics2_one_causal_results,
#'                    SNR=FALSE)
getMatchRecProbDF2 <- function(top_df,stable_df,SNR=TRUE) {
  df_intermediate <- NULL
  # SNR-stratified calculations
  if (SNR) {
    for (snp in 1:3) {
      for (phi in c(0.05,0.1,0.2,0.4)) {
        snr_df_1 <- top_df %>% subset(Phi==phi)
        snr_df_2 <- stable_df %>% subset(Phi==phi)
        string_vec <- mapply(getIntersect2,
                             snr_df_1[[paste0("SNP",snp)]], 
                             snr_df_2[[paste0("SNP",snp)]]) 
        
        true_false_vec <- !is.na(string_vec)
        
        snr_df_1$Match <- true_false_vec
        snr_df_2$Match <- true_false_vec
        snr_df_2$MatchString <- string_vec
        
        match_snr_df_2 <- snr_df_2 %>% subset(Match) # matching top and stable
        nonmatch_snr_df_1 <- snr_df_1 %>% subset(!Match) # non-matching top 
        nonmatch_snr_df_2 <- snr_df_2 %>% subset(!Match) # non-matching stable
        
        match_snr_df_2_contains_causal <- mapply(
          getIntersect,
          match_snr_df_2$CausalSNPs,match_snr_df_2$MatchString
        )  
        nonmatch_snr_df_2_contains_causal <- mapply(
          getIntersect,
          nonmatch_snr_df_2$CausalSNPs,nonmatch_snr_df_2[[paste0("SNP",snp)]]
        )  
        nonmatch_snr_df_1_contains_causal <- mapply(
          getIntersect,
          nonmatch_snr_df_1$CausalSNPs,nonmatch_snr_df_1[[paste0("SNP",snp)]]
        )  
        
        df_intermediate <- rbind(
          df_intermediate, 
          data.frame(Phi=phi,
                     Potential_Set=paste0("Potential Set ", snp),
                     Matching=mean(match_snr_df_2_contains_causal),
                     NonMatching_Top=mean(nonmatch_snr_df_1_contains_causal),
                     NonMatching_Stable=mean(nonmatch_snr_df_2_contains_causal),
                     N_matching=length(match_snr_df_2_contains_causal),
                     N_nonMatching_top=length(nonmatch_snr_df_1_contains_causal),
                     N_nonMatching_stable=length(nonmatch_snr_df_2_contains_causal)))
        
      }
    }
    
    # Gather intermediate results into desired dataframe
    df_to_return <- NULL
    for (ps in 1:3) {
      for (phi in c(0.05,0.1,0.2,0.4)) {
        rel_row <- df_intermediate %>% 
          subset(Phi==phi & Potential_Set==paste0("Potential Set ", ps))
        df_to_return <- 
          rbind(df_to_return, 
                data.frame(PS=rel_row[["Potential_Set"]],
                           SNR=phi/(1-phi),
                           RecProb=rel_row[["Matching"]],
                           BinomSE=sqrt(rel_row[["Matching"]]*
                                          (1-rel_row[["Matching"]])/
                                          rel_row[["N_matching"]]),
                           Type="Matching"))
        df_to_return <- 
          rbind(df_to_return, 
                data.frame(PS=rel_row[["Potential_Set"]],
                           SNR=phi/(1-phi),
                           RecProb=rel_row[["NonMatching_Top"]],
                           BinomSE=sqrt(rel_row[["NonMatching_Top"]]*
                                          (1-rel_row[["NonMatching_Top"]])/
                                          rel_row[["N_nonMatching_top"]]),
                           Type="Non-Matching Top"))
        df_to_return <- 
          rbind(df_to_return, 
                data.frame(PS=rel_row[["Potential_Set"]],
                           SNR=phi/(1-phi),
                           RecProb=rel_row[["NonMatching_Stable"]],
                           BinomSE=sqrt(rel_row[["NonMatching_Stable"]]*
                                          (1-rel_row[["NonMatching_Stable"]])/
                                          rel_row[["N_nonMatching_stable"]]),
                           Type="Non-Matching Stable"))
      }
    }
  } else {
    # SNR-agnostic calculations
    for (snp in 1:3) {
      string_vec <- mapply(getIntersect2,
                           top_df[[paste0("SNP",snp)]], 
                           stable_df[[paste0("SNP",snp)]]) 
      
      true_false_vec <- !is.na(string_vec)
      
      top_df$Match <- true_false_vec
      stable_df$Match <- true_false_vec
      stable_df$MatchString <- string_vec
      
      match_snr_df_2 <- stable_df %>% subset(Match) # matching top and stable
      nonmatch_snr_df_1 <- top_df %>% subset(!Match) # non-matching top 
      nonmatch_snr_df_2 <- stable_df %>% subset(!Match) # non-matching stable
      
      match_snr_df_2_contains_causal <- mapply(
        getIntersect,
        match_snr_df_2$CausalSNPs,match_snr_df_2$MatchString
      )  
      nonmatch_snr_df_2_contains_causal <- mapply(
        getIntersect,
        nonmatch_snr_df_2$CausalSNPs,nonmatch_snr_df_2[[paste0("SNP",snp)]]
      )  
      nonmatch_snr_df_1_contains_causal <- mapply(
        getIntersect,
        nonmatch_snr_df_1$CausalSNPs,nonmatch_snr_df_1[[paste0("SNP",snp)]]
      )  
      
      df_intermediate <- rbind(
        df_intermediate, 
        data.frame(Potential_Set=paste0("Potential Set ", snp),
                   Matching=mean(match_snr_df_2_contains_causal),
                   NonMatching_Top=mean(nonmatch_snr_df_1_contains_causal),
                   NonMatching_Stable=mean(nonmatch_snr_df_2_contains_causal),
                   N_matching=length(match_snr_df_2_contains_causal),
                   N_nonMatching_top=length(nonmatch_snr_df_1_contains_causal),
                   N_nonMatching_stable=length(nonmatch_snr_df_2_contains_causal)))
    }
    
    # Gather intermediate results into desired dataframe
    df_to_return <- NULL
    for (ps in 1:3) {
      rel_row <- df_intermediate %>% 
        subset(Potential_Set==paste0("Potential Set ", ps))
      df_to_return <- 
        rbind(df_to_return, 
              data.frame(PS=rel_row[["Potential_Set"]],
                         RecProb=rel_row[["Matching"]],
                         BinomSE=sqrt(rel_row[["Matching"]]*
                                        (1-rel_row[["Matching"]])/
                                        rel_row[["N_matching"]]),
                         Type="Matching"))
      df_to_return <- 
        rbind(df_to_return, 
              data.frame(PS=rel_row[["Potential_Set"]],
                         RecProb=rel_row[["NonMatching_Top"]],
                         BinomSE=sqrt(rel_row[["NonMatching_Top"]]*
                                        (1-rel_row[["NonMatching_Top"]])/
                                        rel_row[["N_nonMatching_top"]]),
                         Type="Non-Matching Top"))
      df_to_return <- 
        rbind(df_to_return, 
              data.frame(PS=rel_row[["Potential_Set"]],
                         RecProb=rel_row[["NonMatching_Stable"]],
                         BinomSE=sqrt(rel_row[["NonMatching_Stable"]]*
                                        (1-rel_row[["NonMatching_Stable"]])/
                                        rel_row[["N_nonMatching_stable"]]),
                         Type="Non-Matching Stable"))
      
    }
  }
  
  # Return
  return(df_to_return)
}

## getCausalVarCountsDF ========================================================
#' Takes two dataframe of fine-mapping results (top vs stable), and stratifies
#' by whether the top and stable variants match. Computes the number of causal 
#' variants recovered for three cases: matching variants, non-matching 
#' top variant and non-matching stable variant. Uses countIntersect. 
#' Reports Potential Set-specific count distributions. By default, 
#' will split the analysis by SNR used in simulations. If SNR is set to `FALSE`, 
#' it will simply compute the statistics across all SNR scenarios in the 
#' provided datasets. It counts the number of causal variants recovered 
#' from 0 to 2. Relies on countIntersect, getIntersect2
#' getCausalVarCountsDF(top_pics2_one_causal_results,
#'                      stable_pics2_one_causal_results,
#'                      SNR=FALSE)
getCausalVarCountsDF <- function(top_df,
                                 stable_df,
                                 SNR=TRUE,
                                 algo="PICS") {
  # Define keyword for appending to string to return in the Set column
  KEYWORD <- ifelse(algo=="PICS","Potential", "Credible")
  df_to_return <- NULL
  # SNR-stratified calculations
  if (SNR) {
    for (snp in 1:3) {
      for (phi in c(0.05,0.1,0.2,0.4)) {
        snr_df_1 <- top_df %>% subset(Phi==phi)
        snr_df_2 <- stable_df %>% subset(Phi==phi)
        string_vec <- mapply(getIntersect2,
                             snr_df_1[[paste0("SNP",snp)]], 
                             snr_df_2[[paste0("SNP",snp)]]) 
        
        true_false_vec <- !is.na(string_vec)
        
        snr_df_1$Match <- true_false_vec
        snr_df_2$Match <- true_false_vec
        snr_df_2$MatchString <- string_vec
        
        match_snr_df_2 <- snr_df_2 %>% subset(Match) # matching top and stable
        nonmatch_snr_df_1 <- snr_df_1 %>% subset(!Match) # non-matching top 
        nonmatch_snr_df_2 <- snr_df_2 %>% subset(!Match) # non-matching stable
        
        # [!] CHANGE TO COUNT CAUSAL SNPS RECOVERED
        match_snr_df_2_num_causal <- mapply(
          countIntersect,
          match_snr_df_2$CausalSNPs,match_snr_df_2$MatchString
        )  
        nonmatch_snr_df_2_num_causal <- mapply(
          countIntersect,
          nonmatch_snr_df_2$CausalSNPs,nonmatch_snr_df_2[[paste0("SNP",snp)]]
        )  
        nonmatch_snr_df_1_num_causal <- mapply(
          countIntersect,
          nonmatch_snr_df_1$CausalSNPs,nonmatch_snr_df_1[[paste0("SNP",snp)]]
        )  
        
        match_df <- 
          DescTools::MultinomCI(table(match_snr_df_2_num_causal)) %>% 
          as.data.frame() # matching
        nonmatch_2_df <- 
          DescTools::MultinomCI(table(nonmatch_snr_df_2_num_causal)) %>% 
          as.data.frame() # nonmatching stable
        nonmatch_1_df <- 
          DescTools::MultinomCI(table(nonmatch_snr_df_1_num_causal)) %>% 
          as.data.frame() # nonmatching top
        
        # add matching 
        for (i in 1:3) {
          # add matching
          df_to_return <- 
            rbind(df_to_return, 
                  data.frame(Phi=phi,
                             Potential_Set=paste0(KEYWORD," Set ", snp),
                             Type="Matching",
                             NoCausalRec=(i-1),
                             est=ifelse(!identical(match_df[rownames(match_df)==
                                                              (i-1),]$est,numeric(0)),
                                        match_df[rownames(match_df)==(i-1),]$est,0),
                             lwr.ci=ifelse(!identical(match_df[rownames(match_df)==
                                                                 (i-1),]$lwr.ci,numeric(0)),
                                           match_df[rownames(match_df)==(i-1),]$lwr.ci,0),
                             upr.ci=ifelse(!identical(match_df[rownames(match_df)==
                                                                 (i-1),]$upr.ci,numeric(0)),
                                           match_df[rownames(match_df)==(i-1),]$upr.ci,0)))
          # add non-matching stable
          df_to_return <- 
            rbind(df_to_return, 
                  data.frame(Phi=phi,
                             Potential_Set=paste0(KEYWORD," Set ", snp),
                             Type="Non-Matching Stable",
                             NoCausalRec=(i-1),
                             est=ifelse(!identical(nonmatch_2_df[rownames(nonmatch_2_df)==
                                                                   (i-1),]$est,numeric(0)),
                                        nonmatch_2_df[rownames(nonmatch_2_df)==
                                                        (i-1),]$est,0),
                             lwr.ci=ifelse(!identical(nonmatch_2_df[rownames(nonmatch_2_df)==
                                                                      (i-1),]$lwr.ci,numeric(0)),
                                           nonmatch_2_df[rownames(nonmatch_2_df)==
                                                           (i-1),]$lwr.ci,0),
                             upr.ci=ifelse(!identical(nonmatch_2_df[rownames(nonmatch_2_df)==
                                                                      (i-1),]$upr.ci,numeric(0)),
                                           nonmatch_2_df[rownames(nonmatch_2_df)==
                                                           (i-1),]$upr.ci,0)))
          # add non-matching top
          df_to_return <- 
            rbind(df_to_return, 
                  data.frame(Phi=phi,
                             Potential_Set=paste0(KEYWORD," Set ", snp),
                             Type="Non-Matching Top",
                             NoCausalRec=(i-1),
                             est=ifelse(!identical(nonmatch_1_df[rownames(nonmatch_1_df)==
                                                                   (i-1),]$est,numeric(0)),
                                        nonmatch_1_df[rownames(nonmatch_1_df)==
                                                        (i-1),]$est,0),
                             lwr.ci=ifelse(!identical(nonmatch_1_df[rownames(nonmatch_1_df)==
                                                                      (i-1),]$lwr.ci,numeric(0)),
                                           nonmatch_1_df[rownames(nonmatch_1_df)==
                                                           (i-1),]$lwr.ci,0),
                             upr.ci=ifelse(!identical(nonmatch_1_df[rownames(nonmatch_1_df)==
                                                                      (i-1),]$upr.ci,numeric(0)),
                                           nonmatch_1_df[rownames(nonmatch_1_df)==
                                                           (i-1),]$upr.ci,0)))
        }
      }
    }
  } else {
    # SNR-agnostic calculations
    for (snp in 1:3) {
      string_vec <- mapply(getIntersect2,
                           top_df[[paste0("SNP",snp)]], 
                           stable_df[[paste0("SNP",snp)]]) 
      
      true_false_vec <- !is.na(string_vec)
      
      top_df$Match <- true_false_vec
      stable_df$Match <- true_false_vec
      stable_df$MatchString <- string_vec
      
      match_snr_df_2 <- stable_df %>% subset(Match) # matching top and stable
      nonmatch_snr_df_1 <- top_df %>% subset(!Match) # non-matching top 
      nonmatch_snr_df_2 <- stable_df %>% subset(!Match) # non-matching stable
      
      match_snr_df_2_num_causal <- mapply(
        countIntersect,
        match_snr_df_2$CausalSNPs,match_snr_df_2$MatchString
      )  
      nonmatch_snr_df_2_num_causal <- mapply(
        countIntersect,
        nonmatch_snr_df_2$CausalSNPs,nonmatch_snr_df_2[[paste0("SNP",snp)]]
      )  
      nonmatch_snr_df_1_num_causal <- mapply(
        countIntersect,
        nonmatch_snr_df_1$CausalSNPs,nonmatch_snr_df_1[[paste0("SNP",snp)]]
      )  
      
      match_df <- 
        DescTools::MultinomCI(table(match_snr_df_2_num_causal)) %>% 
        as.data.frame() # matching
      nonmatch_2_df <- 
        DescTools::MultinomCI(table(nonmatch_snr_df_2_num_causal)) %>% 
        as.data.frame() # nonmatching stable
      nonmatch_1_df <- 
        DescTools::MultinomCI(table(nonmatch_snr_df_1_num_causal)) %>% 
        as.data.frame() # nonmatching top
      
      # add matching 
      for (i in 1:3) {
        # add matching
        df_to_return <- 
          rbind(df_to_return, 
                data.frame(Phi=phi,
                           Potential_Set=paste0(KEYWORD," Set ", snp),
                           Type="Matching",
                           NoCausalRec=(i-1),
                           est=ifelse(!identical(match_df[rownames(match_df)==
                                                            (i-1),]$est,numeric(0)),
                                      match_df[rownames(match_df)==(i-1),]$est,0),
                           lwr.ci=ifelse(!identical(match_df[rownames(match_df)==
                                                               (i-1),]$lwr.ci,numeric(0)),
                                         match_df[rownames(match_df)==(i-1),]$lwr.ci,0),
                           upr.ci=ifelse(!identical(match_df[rownames(match_df)==
                                                               (i-1),]$upr.ci,numeric(0)),
                                         match_df[rownames(match_df)==(i-1),]$upr.ci,0)))
        # add non-matching stable
        df_to_return <- 
          rbind(df_to_return, 
                data.frame(Phi=phi,
                           Potential_Set=paste0(KEYWORD," Set ", snp),
                           Type="Non-Matching Stable",
                           NoCausalRec=(i-1),
                           est=ifelse(!identical(nonmatch_2_df[rownames(nonmatch_2_df)==
                                                                 (i-1),]$est,numeric(0)),
                                      nonmatch_2_df[rownames(nonmatch_2_df)==
                                                      (i-1),]$est,0),
                           lwr.ci=ifelse(!identical(nonmatch_2_df[rownames(nonmatch_2_df)==
                                                                    (i-1),]$lwr.ci,numeric(0)),
                                         nonmatch_2_df[rownames(nonmatch_2_df)==
                                                         (i-1),]$lwr.ci,0),
                           upr.ci=ifelse(!identical(nonmatch_2_df[rownames(nonmatch_2_df)==
                                                                    (i-1),]$upr.ci,numeric(0)),
                                         nonmatch_2_df[rownames(nonmatch_2_df)==
                                                         (i-1),]$upr.ci,0)))
        # add non-matching top
        df_to_return <- 
          rbind(df_to_return, 
                data.frame(Phi=phi,
                           Potential_Set=paste0(KEYWORD," Set ", snp),
                           Type="Non-Matching Top",
                           NoCausalRec=(i-1),
                           est=ifelse(!identical(nonmatch_1_df[rownames(nonmatch_1_df)==
                                                                 (i-1),]$est,numeric(0)),
                                      nonmatch_1_df[rownames(nonmatch_1_df)==
                                                      (i-1),]$est,0),
                           lwr.ci=ifelse(!identical(nonmatch_1_df[rownames(nonmatch_1_df)==
                                                                    (i-1),]$lwr.ci,numeric(0)),
                                         nonmatch_1_df[rownames(nonmatch_1_df)==
                                                         (i-1),]$lwr.ci,0),
                           upr.ci=ifelse(!identical(nonmatch_1_df[rownames(nonmatch_1_df)==
                                                                    (i-1),]$upr.ci,numeric(0)),
                                         nonmatch_1_df[rownames(nonmatch_1_df)==
                                                         (i-1),]$upr.ci,0)))
      }
    }
  }
  # Return
  return(df_to_return)
}

## getCausalVarCountsDF2 =======================================================
#' Takes two dataframe of fine-mapping results (top vs stable), and stratifies
#' by whether the top and stable variants match. Computes the number of causal 
#' variants recovered for three cases: matching variants, non-matching 
#' top variant and non-matching stable variant. Uses countIntersect. Aggregates 
#' counts across Potential Sets, unlike getCausalVarCountsDF. It counts the 
#' number of causal variants recovered from 0 to 2. Also stratifies
#' results by the number of matching variants between Top and Stable PICS. This
#' is to answer the question, "How frequently are causal variants recovered
#' CONDITIONED ON the number of matching variants?"
#' getCausalVarCountsDF2(top_pics2_one_causal_results,
#'                       stable_pics2_one_causal_results)
getCausalVarCountsDF2 <- function(top_df,stable_df) {
  top_df_to_return <- NULL
  stable_df_to_return <- NULL
  for (phi in c(0.05,0.1,0.2,0.4)) {
    snr_df_1 <- top_df %>% subset(Phi==phi)
    snr_df_2 <- stable_df %>% subset(Phi==phi)
    
    # Get matching SNP count between top and stable 
    for (snp in 1:3) {
      string_vec <- mapply(getIntersect2,
                           snr_df_1[[paste0("SNP",snp)]], 
                           snr_df_2[[paste0("SNP",snp)]]) 
      
      true_false_vec <- !is.na(string_vec)
      
      snr_df_1[[paste0("Match",snp)]] <- true_false_vec
      snr_df_2[[paste0("Match",snp)]] <- true_false_vec
      snr_df_2[[paste0("Match",snp,"_String")]] <- string_vec
    }
    
    # Count the number of potential sets with matching SNPs
    snr_df_2[["TotalMatchCount"]] <- snr_df_2[[paste0("Match1")]]+
      snr_df_2[[paste0("Match2")]]+snr_df_2[[paste0("Match3")]]
    snr_df_1[["TotalMatchCount"]] <- snr_df_1[[paste0("Match1")]]+
      snr_df_1[[paste0("Match2")]]+snr_df_1[[paste0("Match3")]]
    
    # Collect all SNPs for stable and for top PICS
    # If interested in evaluating the recovery of matching SNPs rather
    # than the SNPs returned by each individual algorithm, use Match1_String
    # in place of SNP1, Match2_String in place of SNP2, etc.
    stable_ps1_ps2 <- mapply(function(x,y){
      x_new <- paste(x,collapse=",")
      y_new <- paste(y,collapse=",")
      paste(x_new,y_new,sep=",")}, 
      snr_df_2[["SNP1"]],
      snr_df_2[["SNP2"]])
    stable_ps1_ps2_ps3 <- mapply(function(x,y){
      x_new <- paste(x,collapse=",")
      y_new <- paste(y,collapse=",")
      paste(x_new,y_new,sep=",")},
      stable_ps1_ps2,
      snr_df_2[["SNP3"]])
    snr_df_2[["AllMatchingSNPs"]] <- stable_ps1_ps2_ps3 # this is all SNPs across 3 potential sets
    
    top_ps1_ps2 <- mapply(function(x,y){
      x_new <- paste(x,collapse=",")
      y_new <- paste(y,collapse=",")
      paste(x_new,y_new,sep=",")}, 
      snr_df_1[["SNP1"]],
      snr_df_1[["SNP2"]])
    top_ps1_ps2_ps3 <- mapply(function(x,y){
      x_new <- paste(x,collapse=",")
      y_new <- paste(y,collapse=",")
      paste(x_new,y_new,sep=",")},
      top_ps1_ps2,
      snr_df_1[["SNP3"]]) 
    snr_df_1[["AllMatchingSNPs"]] <- top_ps1_ps2_ps3 
    
    # For (i-1) matches, (ranges from 0 to #[Potential Sets])
    for (i in 1:4) {
      # Look at stable SNP dataframes
      match_i_snr_df_2 <- snr_df_2 %>% 
        subset(TotalMatchCount==(i-1)) # (i-1) matches
      
      # Look at top SNP dataframes
      match_i_snr_df_1 <- snr_df_1 %>% 
        subset(TotalMatchCount==(i-1)) # (i-1) matches
      
      # ... then tabulate number of causal SNPs recovered for each
      match_i_snr_df_2_num_causal <- mapply(
        countIntersect,
        match_i_snr_df_2$CausalSNPs,match_i_snr_df_2$AllMatchingSNPs
      )  
      match_i_snr_df_1_num_causal <- mapply(
        countIntersect,
        match_i_snr_df_1$CausalSNPs,match_i_snr_df_1$AllMatchingSNPs
      )
      
      # Get causal SNP count distribution
      match_i_snr_df_1_dist_df <- 
        DescTools::MultinomCI(table(match_i_snr_df_1_num_causal)) %>% as.data.frame()
      match_i_snr_df_2_dist_df <- 
        DescTools::MultinomCI(table(match_i_snr_df_2_num_causal)) %>% as.data.frame()
      
      # Add to dataframe to return
      for (no_causal_var in 0:2) {
        if (!identical(match_i_snr_df_1_dist_df[rownames(match_i_snr_df_1_dist_df)==no_causal_var,]$est,
                       numeric(0))) {
          top_df_to_return <- 
            rbind(top_df_to_return,
                  data.frame(Phi=phi,
                             Algorithm="Top PICS",
                             NumMatchingVar=(i-1),
                             NoCausalRec=no_causal_var,
                             est=match_i_snr_df_1_dist_df[rownames(match_i_snr_df_1_dist_df)==no_causal_var,]$est,
                             lwr.ci=match_i_snr_df_1_dist_df[rownames(match_i_snr_df_1_dist_df)==no_causal_var,]$lwr.ci,
                             upr.ci=match_i_snr_df_1_dist_df[rownames(match_i_snr_df_1_dist_df)==no_causal_var,]$upr.ci))
        } else {
          top_df_to_return <- rbind(top_df_to_return,
                                    data.frame(Phi=phi,
                                               Algorithm="Top PICS",
                                               NumMatchingVar=(i-1),
                                               NoCausalRec=no_causal_var,
                                               est=0,
                                               lwr.ci=0,
                                               upr.ci=0))
          
        }
        
        if (!identical(match_i_snr_df_2_dist_df[rownames(match_i_snr_df_2_dist_df)==no_causal_var,]$est,
                       numeric(0))) {
          stable_df_to_return <- 
            rbind(stable_df_to_return,
                  data.frame(Phi=phi,
                             Algorithm="Stable PICS",
                             NumMatchingVar=(i-1),
                             NoCausalRec=no_causal_var,
                             est=match_i_snr_df_2_dist_df[rownames(match_i_snr_df_2_dist_df)==no_causal_var,]$est,
                             lwr.ci=match_i_snr_df_2_dist_df[rownames(match_i_snr_df_2_dist_df)==no_causal_var,]$lwr.ci,
                             upr.ci=match_i_snr_df_2_dist_df[rownames(match_i_snr_df_2_dist_df)==no_causal_var,]$upr.ci))
        } else {
          stable_df_to_return <- rbind(stable_df_to_return,
                                       data.frame(Phi=phi,
                                                  Algorithm="Stable PICS",
                                                  NumMatchingVar=(i-1),
                                                  NoCausalRec=no_causal_var,
                                                  est=0,
                                                  lwr.ci=0,
                                                  upr.ci=0))
        }
      }
    }
  }
  
  list_to_return <- list(Stable=stable_df_to_return,
                         Top=top_df_to_return)
  # Return
  return(list_to_return)
}

## getCausalVarCountsDF3 =======================================================
#' Takes two dataframe of fine-mapping results (top vs stable), and stratifies
#' by whether the top and stable variants match. Computes the number of causal 
#' variants recovered for three cases: matching variants, non-matching 
#' top variant and non-matching stable variant. Uses countIntersect. 
#' Reports Potential Set-specific count distributions. By default, 
#' will split the analysis by SNR used in simulations. If SNR is set to `FALSE`, 
#' it will simply compute the statistics across all SNR scenarios in the 
#' provided datasets. It counts the number of causal variants recovered 
#' from 0 to 3. Relies on countIntersect, getIntersect2
#' getCausalVarCountsDF3(top_pics2_one_causal_results,
#'                      stable_pics2_one_causal_results,
#'                      SNR=FALSE)
getCausalVarCountsDF3 <- function(top_df,
                                  stable_df,
                                  SNR=TRUE,
                                  algo="PICS") {
  # Define keyword for appending to string to return in the Set column
  KEYWORD <- ifelse(algo=="PICS","Potential", "Credible")
  df_to_return <- NULL
  if (SNR) {
    for (snp in 1:3) {
      for (phi in c(0.05,0.1,0.2,0.4)) {
        snr_df_1 <- top_df %>% subset(Phi==phi)
        snr_df_2 <- stable_df %>% subset(Phi==phi)
        string_vec <- mapply(getIntersect2,
                             snr_df_1[[paste0("SNP",snp)]], 
                             snr_df_2[[paste0("SNP",snp)]]) 
        
        true_false_vec <- !is.na(string_vec)
        
        snr_df_1$Match <- true_false_vec
        snr_df_2$Match <- true_false_vec
        snr_df_2$MatchString <- string_vec
        
        match_snr_df_2 <- snr_df_2 %>% subset(Match) # matching top and stable
        nonmatch_snr_df_1 <- snr_df_1 %>% subset(!Match) # non-matching top 
        nonmatch_snr_df_2 <- snr_df_2 %>% subset(!Match) # non-matching stable
        
        match_snr_df_2_num_causal <- mapply(
          countIntersect,
          match_snr_df_2$CausalSNPs,match_snr_df_2$MatchString
        )  
        nonmatch_snr_df_2_num_causal <- mapply(
          countIntersect,
          nonmatch_snr_df_2$CausalSNPs,nonmatch_snr_df_2[[paste0("SNP",snp)]]
        )  
        nonmatch_snr_df_1_num_causal <- mapply(
          countIntersect,
          nonmatch_snr_df_1$CausalSNPs,nonmatch_snr_df_1[[paste0("SNP",snp)]]
        )  
        
        match_df <- 
          DescTools::MultinomCI(table(match_snr_df_2_num_causal)) %>% 
          as.data.frame() # matching
        nonmatch_2_df <- 
          DescTools::MultinomCI(table(nonmatch_snr_df_2_num_causal)) %>% 
          as.data.frame() # nonmatching stable
        nonmatch_1_df <- 
          DescTools::MultinomCI(table(nonmatch_snr_df_1_num_causal)) %>% 
          as.data.frame() # nonmatching top
        
        # for each number of causal variants recovered
        for (i in 1:4) {
          # add matching
          df_to_return <- 
            rbind(df_to_return, 
                  data.frame(Phi=phi,
                             Potential_Set=paste0(KEYWORD," Set ", snp),
                             Type="Matching",
                             NoCausalRec=(i-1),
                             est=ifelse(!identical(match_df[rownames(match_df)==
                                                              (i-1),]$est,numeric(0)),
                                        match_df[rownames(match_df)==(i-1),]$est,0),
                             lwr.ci=ifelse(!identical(match_df[rownames(match_df)==
                                                                 (i-1),]$lwr.ci,numeric(0)),
                                           match_df[rownames(match_df)==(i-1),]$lwr.ci,0),
                             upr.ci=ifelse(!identical(match_df[rownames(match_df)==
                                                                 (i-1),]$upr.ci,numeric(0)),
                                           match_df[rownames(match_df)==
                                                      (i-1),]$upr.ci,0)))
          # add non-matching stable
          df_to_return <- 
            rbind(df_to_return, 
                  data.frame(Phi=phi,
                             Potential_Set=paste0(KEYWORD," Set ", snp),
                             Type="Non-Matching Stable",
                             NoCausalRec=(i-1),
                             est=ifelse(!identical(nonmatch_2_df[rownames(nonmatch_2_df)==
                                                                   (i-1),]$est,numeric(0)),
                                        nonmatch_2_df[rownames(nonmatch_2_df)==
                                                        (i-1),]$est,0),
                             lwr.ci=ifelse(!identical(nonmatch_2_df[rownames(nonmatch_2_df)==
                                                                      (i-1),]$lwr.ci,numeric(0)),
                                           nonmatch_2_df[rownames(nonmatch_2_df)==
                                                           (i-1),]$lwr.ci,0),
                             upr.ci=ifelse(!identical(nonmatch_2_df[rownames(nonmatch_2_df)==
                                                                      (i-1),]$upr.ci,numeric(0)),
                                           nonmatch_2_df[rownames(nonmatch_2_df)==
                                                           (i-1),]$upr.ci,0)))
          # add non-matching top
          df_to_return <- 
            rbind(df_to_return, 
                  data.frame(Phi=phi,
                             Potential_Set=paste0(KEYWORD," Set ", snp),
                             Type="Non-Matching Top",
                             NoCausalRec=(i-1),
                             est=ifelse(!identical(nonmatch_1_df[rownames(nonmatch_1_df)==
                                                                   (i-1),]$est,numeric(0)),
                                        nonmatch_1_df[rownames(nonmatch_1_df)==
                                                        (i-1),]$est,0),
                             lwr.ci=ifelse(!identical(nonmatch_1_df[rownames(nonmatch_1_df)==
                                                                      (i-1),]$lwr.ci,numeric(0)),
                                           nonmatch_1_df[rownames(nonmatch_1_df)==
                                                           (i-1),]$lwr.ci,0),
                             upr.ci=ifelse(!identical(nonmatch_1_df[rownames(nonmatch_1_df)==
                                                                      (i-1),]$upr.ci,numeric(0)),
                                           nonmatch_1_df[rownames(nonmatch_1_df)==
                                                           (i-1),]$upr.ci,0)))
        }
      }
    }
  } else {
    for (snp in 1:3) {
      string_vec <- mapply(getIntersect2,
                           top_df[[paste0("SNP",snp)]], 
                           stable_df[[paste0("SNP",snp)]]) 
      
      true_false_vec <- !is.na(string_vec)
      
      top_df$Match <- true_false_vec
      stable_df$Match <- true_false_vec
      stable_df$MatchString <- string_vec
      
      match_snr_df_2 <- stable_df %>% subset(Match) # matching top and stable
      nonmatch_snr_df_1 <- top_df %>% subset(!Match) # non-matching top 
      nonmatch_snr_df_2 <- stable_df %>% subset(!Match) # non-matching stable
      
      match_snr_df_2_num_causal <- mapply(
        countIntersect,
        match_snr_df_2$CausalSNPs,match_snr_df_2$MatchString
      )  
      nonmatch_snr_df_2_num_causal <- mapply(
        countIntersect,
        nonmatch_snr_df_2$CausalSNPs,nonmatch_snr_df_2[[paste0("SNP",snp)]]
      )  
      nonmatch_snr_df_1_num_causal <- mapply(
        countIntersect,
        nonmatch_snr_df_1$CausalSNPs,nonmatch_snr_df_1[[paste0("SNP",snp)]]
      )  
      
      match_df <- 
        DescTools::MultinomCI(table(match_snr_df_2_num_causal)) %>% 
        as.data.frame() # matching
      nonmatch_2_df <- 
        DescTools::MultinomCI(table(nonmatch_snr_df_2_num_causal)) %>% 
        as.data.frame() # nonmatching stable
      nonmatch_1_df <- 
        DescTools::MultinomCI(table(nonmatch_snr_df_1_num_causal)) %>% 
        as.data.frame() # nonmatching top
      
      # for each number of causal variants recovered
      for (i in 1:4) {
        # add matching
        df_to_return <- 
          rbind(df_to_return, 
                data.frame(Phi=phi,
                           Potential_Set=paste0(KEYWORD," Set ", snp),
                           Type="Matching",
                           NoCausalRec=(i-1),
                           est=ifelse(!identical(match_df[rownames(match_df)==
                                                            (i-1),]$est,numeric(0)),
                                      match_df[rownames(match_df)==(i-1),]$est,0),
                           lwr.ci=ifelse(!identical(match_df[rownames(match_df)==
                                                               (i-1),]$lwr.ci,numeric(0)),
                                         match_df[rownames(match_df)==(i-1),]$lwr.ci,0),
                           upr.ci=ifelse(!identical(match_df[rownames(match_df)==
                                                               (i-1),]$upr.ci,numeric(0)),
                                         match_df[rownames(match_df)==
                                                    (i-1),]$upr.ci,0)))
        # add non-matching stable
        df_to_return <- 
          rbind(df_to_return, 
                data.frame(Phi=phi,
                           Potential_Set=paste0(KEYWORD," Set ", snp),
                           Type="Non-Matching Stable",
                           NoCausalRec=(i-1),
                           est=ifelse(!identical(nonmatch_2_df[rownames(nonmatch_2_df)==
                                                                 (i-1),]$est,numeric(0)),
                                      nonmatch_2_df[rownames(nonmatch_2_df)==
                                                      (i-1),]$est,0),
                           lwr.ci=ifelse(!identical(nonmatch_2_df[rownames(nonmatch_2_df)==
                                                                    (i-1),]$lwr.ci,numeric(0)),
                                         nonmatch_2_df[rownames(nonmatch_2_df)==
                                                         (i-1),]$lwr.ci,0),
                           upr.ci=ifelse(!identical(nonmatch_2_df[rownames(nonmatch_2_df)==
                                                                    (i-1),]$upr.ci,numeric(0)),
                                         nonmatch_2_df[rownames(nonmatch_2_df)==
                                                         (i-1),]$upr.ci,0)))
        # add non-matching top
        df_to_return <- 
          rbind(df_to_return, 
                data.frame(Phi=phi,
                           Potential_Set=paste0(KEYWORD," Set ", snp),
                           Type="Non-Matching Top",
                           NoCausalRec=(i-1),
                           est=ifelse(!identical(nonmatch_1_df[rownames(nonmatch_1_df)==
                                                                 (i-1),]$est,numeric(0)),
                                      nonmatch_1_df[rownames(nonmatch_1_df)==
                                                      (i-1),]$est,0),
                           lwr.ci=ifelse(!identical(nonmatch_1_df[rownames(nonmatch_1_df)==
                                                                    (i-1),]$lwr.ci,numeric(0)),
                                         nonmatch_1_df[rownames(nonmatch_1_df)==
                                                         (i-1),]$lwr.ci,0),
                           upr.ci=ifelse(!identical(nonmatch_1_df[rownames(nonmatch_1_df)==
                                                                    (i-1),]$upr.ci,numeric(0)),
                                         nonmatch_1_df[rownames(nonmatch_1_df)==
                                                         (i-1),]$upr.ci,0)))
      }
    }
  }
  
  # Return
  return(df_to_return)
}

## getCausalVarCountsDF4 =======================================================
# Takes two dataframe of fine-mapping results (top vs stable), and stratifies
# by whether the top and stable variants match. Computes the number of causal 
# variants recovered for three cases: matching variants, non-matching 
# top variant and non-matching stable variant. Uses countIntersect. Aggregates 
#' counts across Potential Sets, unlike getCausalVarCountsDF3. It counts the 
#' number of causal variants recovered from 0 to 3. Also stratifies
#' results by the number of matching variants between Top and Stable PICS. This
#' is to answer the question, "How frequently are causal variants recovered
#' CONDITIONED ON the number of matching variants?"
# Aggregates counts across Potential Sets.
# getCausalVarCountsDF4(top_pics2_one_causal_results,
#                       stable_pics2_one_causal_results)
getCausalVarCountsDF4 <- function(top_df,stable_df) {
  top_df_to_return <- NULL
  stable_df_to_return <- NULL
  for (phi in c(0.05,0.1,0.2,0.4)) {
    snr_df_1 <- top_df %>% subset(Phi==phi)
    snr_df_2 <- stable_df %>% subset(Phi==phi)
    
    # Get matching SNPs between top and stable 
    for (snp in 1:3) {
      string_vec <- mapply(getIntersect2,
                           snr_df_1[[paste0("SNP",snp)]], 
                           snr_df_2[[paste0("SNP",snp)]]) 
      
      true_false_vec <- !is.na(string_vec)
      
      snr_df_1[[paste0("Match",snp)]] <- true_false_vec
      snr_df_2[[paste0("Match",snp)]] <- true_false_vec
      snr_df_2[[paste0("Match",snp,"_String")]] <- string_vec
    }
    
    snr_df_2[["TotalMatchCount"]] <- snr_df_2[[paste0("Match1")]]+
      snr_df_2[[paste0("Match2")]]+snr_df_2[[paste0("Match3")]]
    snr_df_1[["TotalMatchCount"]] <- snr_df_1[[paste0("Match1")]]+
      snr_df_1[[paste0("Match2")]]+snr_df_1[[paste0("Match3")]]
    
    stable_ps1_ps2 <- mapply(function(x,y){
      x_new <- paste(x,collapse=",")
      y_new <- paste(y,collapse=",")
      paste(x_new,y_new,sep=",")}, 
      snr_df_2[["SNP1"]],
      snr_df_2[["SNP2"]])
    stable_ps1_ps2_ps3 <- mapply(function(x,y){
      x_new <- paste(x,collapse=",")
      y_new <- paste(y,collapse=",")
      paste(x_new,y_new,sep=",")},
      stable_ps1_ps2,
      snr_df_2[["SNP3"]])
    snr_df_2[["AllMatchingSNPs"]] <- stable_ps1_ps2_ps3 
    
    top_ps1_ps2 <- mapply(function(x,y){
      x_new <- paste(x,collapse=",")
      y_new <- paste(y,collapse=",")
      paste(x_new,y_new,sep=",")}, 
      snr_df_1[["SNP1"]],
      snr_df_1[["SNP2"]])
    top_ps1_ps2_ps3 <- mapply(function(x,y){
      x_new <- paste(x,collapse=",")
      y_new <- paste(y,collapse=",")
      paste(x_new,y_new,sep=",")},
      top_ps1_ps2,
      snr_df_1[["SNP3"]])
    snr_df_1[["AllMatchingSNPs"]] <- top_ps1_ps2_ps3 
    
    # For (i-1) matches, (ranges from 0 to #[Potential Sets])
    for (i in 1:4) {
      # Look at stable SNP dataframes
      match_i_snr_df_2 <- snr_df_2 %>% 
        subset(TotalMatchCount==(i-1)) # (i-1) matches
      
      # Look at top SNP dataframes
      match_i_snr_df_1 <- snr_df_1 %>% 
        subset(TotalMatchCount==(i-1)) # (i-1) matches
      
      # ... then tabulate number of causal SNPs recovered for each
      match_i_snr_df_2_num_causal <- mapply(
        countIntersect,
        match_i_snr_df_2$CausalSNPs,match_i_snr_df_2$AllMatchingSNPs
      )  
      match_i_snr_df_1_num_causal <- mapply(
        countIntersect,
        match_i_snr_df_1$CausalSNPs,match_i_snr_df_1$AllMatchingSNPs
      )
      
      # Get causal SNP count distribution
      match_i_snr_df_1_dist_df <- 
        DescTools::MultinomCI(table(match_i_snr_df_1_num_causal)) %>% as.data.frame()
      match_i_snr_df_2_dist_df <- 
        DescTools::MultinomCI(table(match_i_snr_df_2_num_causal)) %>% as.data.frame()
      
      # Add to dataframe to return
      for (no_causal_var in 0:3) {
        if (!identical(match_i_snr_df_1_dist_df[rownames(match_i_snr_df_1_dist_df)==no_causal_var,]$est,
                       numeric(0))) {
          top_df_to_return <- 
            rbind(top_df_to_return,
                  data.frame(Phi=phi,
                             Algorithm="Top PICS",
                             NumMatchingVar=(i-1),
                             NoCausalRec=no_causal_var,
                             est=match_i_snr_df_1_dist_df[rownames(match_i_snr_df_1_dist_df)==no_causal_var,]$est,
                             lwr.ci=match_i_snr_df_1_dist_df[rownames(match_i_snr_df_1_dist_df)==no_causal_var,]$lwr.ci,
                             upr.ci=match_i_snr_df_1_dist_df[rownames(match_i_snr_df_1_dist_df)==no_causal_var,]$upr.ci))
        } else {
          top_df_to_return <- rbind(top_df_to_return,
                                    data.frame(Phi=phi,
                                               Algorithm="Top PICS",
                                               NumMatchingVar=(i-1),
                                               NoCausalRec=no_causal_var,
                                               est=0,
                                               lwr.ci=0,
                                               upr.ci=0))
          
        }
        
        if (!identical(match_i_snr_df_2_dist_df[rownames(match_i_snr_df_2_dist_df)==no_causal_var,]$est,
                       numeric(0))) {
          stable_df_to_return <- 
            rbind(stable_df_to_return,
                  data.frame(Phi=phi,
                             Algorithm="Stable PICS",
                             NumMatchingVar=(i-1),
                             NoCausalRec=no_causal_var,
                             est=match_i_snr_df_2_dist_df[rownames(match_i_snr_df_2_dist_df)==no_causal_var,]$est,
                             lwr.ci=match_i_snr_df_2_dist_df[rownames(match_i_snr_df_2_dist_df)==no_causal_var,]$lwr.ci,
                             upr.ci=match_i_snr_df_2_dist_df[rownames(match_i_snr_df_2_dist_df)==no_causal_var,]$upr.ci))
        } else {
          stable_df_to_return <- rbind(stable_df_to_return,
                                       data.frame(Phi=phi,
                                                  Algorithm="Stable PICS",
                                                  NumMatchingVar=(i-1),
                                                  NoCausalRec=no_causal_var,
                                                  est=0,
                                                  lwr.ci=0,
                                                  upr.ci=0))
        }
      }
    }
  }
  
  list_to_return <- list(Stable=stable_df_to_return,
                         Top=top_df_to_return)
  # Return
  return(list_to_return)
}

## getMatchStats ===============================================================
#' Takes two dataframe of fine-mapping results, one focal and another, and 
#' computes matching variant statistics for each potential set in the focal result 
#' and all potential sets in the other result file. 
#' getMatchStats(top_pics2_one_causal_results,stable_pics2_one_causal_results)
getMatchStats <- function(focal_df,other_df) {
  # For each potential set,
  for (s in 1:3) {
    # Count the number of intersecting variants with focal df
    focal1_other <- sapply(mapply(getIntersect2,focal_df[["SNP1"]],
                                  other_df[[paste0("SNP",s)]]), 
                           function(x) {paste(x,collapse=",")})
    focal2_other <- sapply(mapply(getIntersect2,focal_df[["SNP2"]],
                                  other_df[[paste0("SNP",s)]]),
                           function(x) {paste(x,collapse=",")})
    focal3_other <- sapply(mapply(getIntersect2,focal_df[["SNP3"]],
                                  other_df[[paste0("SNP",s)]]),
                           function(x) {paste(x,collapse=",")})
    focal_df[[paste0("int_focal_SNP1_other_SNP",s)]] <- focal1_other
    focal_df[[paste0("int_focal_SNP2_other_SNP",s)]] <- focal2_other
    focal_df[[paste0("int_focal_SNP3_other_SNP",s)]] <- focal3_other
    
    focal_df[[paste0("focal_SNP1_match_other_SNP",s)]] <- sapply(
      focal1_other,function(x){!identical(x,"NA")})
    focal_df[[paste0("focal_SNP2_match_other_SNP",s)]] <- sapply(
      focal2_other,function(x){!identical(x,"NA")})
    focal_df[[paste0("focal_SNP3_match_other_SNP",s)]] <- sapply(
      focal3_other,function(x){!identical(x,"NA")})
  }
  
  # Return
  return(focal_df)
}


#### Main Text Figures ---------------------------------------------------------
## Fig A. ======================================================================
plain_pics_res_string_vec <- c("PICS2_Chr1_plain_results",
                               "PICS2_Chr20_plain_results",
                               "PICS2_Chr21_plain_results",
                               "PICS2_Chr22_plain_results")
stable_pics_res_string_vec <- c("PICS2_Chr1_results",
                                "PICS2_Chr20_results",
                                "PICS2_Chr21_results",
                                "PICS2_Chr22_results")

## Initialize empty dataframes
big_plain_pics_df <- NULL
big_stable_pics_df <- NULL

## t=8
plain_pics_df <- do.call(
  rbind,
  lapply(paste0(t8_het_env_res_dir,plain_pics_res_string_vec,".csv"),
         function(x)readr::read_csv(x)))
plain_pics_df$MODEL <- rep("t=8",nrow(plain_pics_df))
stable_pics_df <- do.call(
  rbind,
  lapply(paste0(t8_het_env_res_dir,stable_pics_res_string_vec,".csv"),
         function(x)readr::read_csv(x)))
stable_pics_df$MODEL <- rep("t=8",nrow(stable_pics_df))
big_plain_pics_df <- rbind(big_plain_pics_df,plain_pics_df)
big_stable_pics_df <- rbind(big_stable_pics_df,stable_pics_df)

## t=16
plain_pics_df <- do.call(
  rbind,
  lapply(paste0(t16_het_env_res_dir,plain_pics_res_string_vec,".csv"),
         function(x)readr::read_csv(x)))
plain_pics_df$MODEL <- rep("t=16",nrow(plain_pics_df))
stable_pics_df <- do.call(
  rbind,
  lapply(paste0(t16_het_env_res_dir,stable_pics_res_string_vec,".csv"),
         function(x)readr::read_csv(x)))
stable_pics_df$MODEL <- rep("t=16",nrow(stable_pics_df))
big_plain_pics_df <- rbind(big_plain_pics_df,plain_pics_df)
big_stable_pics_df <- rbind(big_stable_pics_df,stable_pics_df)

## t=128
plain_pics_df <- do.call(
  rbind,
  lapply(paste0(t128_het_env_res_dir,plain_pics_res_string_vec,".csv"),
         function(x)readr::read_csv(x)))
plain_pics_df$MODEL <- rep("t=128",nrow(plain_pics_df))
stable_pics_df <- do.call(
  rbind,
  lapply(paste0(t128_het_env_res_dir,stable_pics_res_string_vec,".csv"),
         function(x)readr::read_csv(x)))
stable_pics_df$MODEL <- rep("t=128",nrow(stable_pics_df))
big_plain_pics_df <- rbind(big_plain_pics_df,plain_pics_df)
big_stable_pics_df <- rbind(big_stable_pics_df,stable_pics_df)

## t=256
plain_pics_df <- do.call(
  rbind,
  lapply(paste0(t256_het_env_res_dir,plain_pics_res_string_vec,".csv"),
         function(x)readr::read_csv(x)))
plain_pics_df$MODEL <- rep("t=256",nrow(plain_pics_df))
stable_pics_df <- do.call(
  rbind,
  lapply(paste0(t256_het_env_res_dir,stable_pics_res_string_vec,".csv"),
         function(x)readr::read_csv(x)))
stable_pics_df$MODEL <- rep("t=256",nrow(stable_pics_df))
big_plain_pics_df <- rbind(big_plain_pics_df,plain_pics_df)
big_stable_pics_df <- rbind(big_stable_pics_df,stable_pics_df)

## |i-3|
plain_pics_df <- do.call(
  rbind,
  lapply(paste0(im3_het_env_res_dir,plain_pics_res_string_vec,".csv"),
         function(x)readr::read_csv(x)))
plain_pics_df$MODEL <- rep("|i-3|",nrow(plain_pics_df))
stable_pics_df <- do.call(
  rbind,
  lapply(paste0(im3_het_env_res_dir,stable_pics_res_string_vec,".csv"),
         function(x)readr::read_csv(x)))
stable_pics_df$MODEL <- rep("|i-3|",nrow(stable_pics_df))
big_plain_pics_df <- rbind(big_plain_pics_df,plain_pics_df)
big_stable_pics_df <- rbind(big_stable_pics_df,stable_pics_df)

## i==3
plain_pics_df <- do.call(
  rbind,
  lapply(paste0(ie3_het_env_res_dir,plain_pics_res_string_vec,".csv"),
         function(x)readr::read_csv(x)))
plain_pics_df$MODEL <- rep("i=3",nrow(plain_pics_df))
stable_pics_df <- do.call(
  rbind,
  lapply(paste0(ie3_het_env_res_dir,stable_pics_res_string_vec,".csv"),
         function(x)readr::read_csv(x)))
stable_pics_df$MODEL <- rep("i=3",nrow(stable_pics_df))
big_plain_pics_df <- rbind(big_plain_pics_df,plain_pics_df)
big_stable_pics_df <- rbind(big_stable_pics_df,stable_pics_df)

## One Causal SNP results 
big_env_het_plain_pics_one_causal_results <- big_plain_pics_df %>% 
  subset(S==1)
big_env_het_stable_pics_one_causal_results <- big_stable_pics_df %>% 
  subset(S==1)

## Two Causal SNP results 
big_env_het_plain_pics_two_causal_results <- big_plain_pics_df %>% 
  subset(S==2)
big_env_het_stable_pics_two_causal_results <- big_stable_pics_df %>% 
  subset(S==2)

## Three Causal SNP results 
big_env_het_plain_pics_three_causal_results <- big_plain_pics_df %>% 
  subset(S==3)
big_env_het_stable_pics_three_causal_results <- big_stable_pics_df %>% 
  subset(S==3)

## Create dataframe for plotting
# Compute whether PS1 Contains at least one causal variant (Stable PICS)
big_env_het_stable_pics_one_causal_results$PS1ContainsCausal <- mapply(
  getIntersect,
  big_env_het_stable_pics_one_causal_results$CausalSNPs,
  big_env_het_stable_pics_one_causal_results$SNP1
) 

big_env_het_stable_pics_two_causal_results$PS1ContainsCausal <- mapply(
  getIntersect,
  big_env_het_stable_pics_two_causal_results$CausalSNPs,
  big_env_het_stable_pics_two_causal_results$SNP1
)

big_env_het_stable_pics_three_causal_results$PS1ContainsCausal <- mapply(
  getIntersect,
  big_env_het_stable_pics_three_causal_results$CausalSNPs,
  big_env_het_stable_pics_three_causal_results$SNP1
)

rec_count_vec <- c(as.numeric(table(big_env_het_stable_pics_one_causal_results$PS1ContainsCausal))[2],
                  as.numeric(table(big_env_het_stable_pics_two_causal_results$PS1ContainsCausal))[2],
                  as.numeric(table(big_env_het_stable_pics_three_causal_results$PS1ContainsCausal))[2])
rec_prob_vec <- rec_count_vec/480
rec_prob_se <- sqrt(rec_prob_vec*(1-rec_prob_vec)/480)
stable_pics_plotting_df <- data.frame(Algorithm=rep("Stable PICS",3),
                                     No_Causal=c(1,2,3),
                                     Rec_Prob=rec_prob_vec,
                                     BinomSE=rec_prob_se)

# Compute whether PS1 Contains at least one causal variant (Plain PICS)
big_env_het_plain_pics_one_causal_results$PS1ContainsCausal <- mapply(
  getIntersect,
  big_env_het_plain_pics_one_causal_results$CausalSNPs,
  big_env_het_plain_pics_one_causal_results$SNP1
) 

big_env_het_plain_pics_two_causal_results$PS1ContainsCausal <- mapply(
  getIntersect,
  big_env_het_plain_pics_two_causal_results$CausalSNPs,
  big_env_het_plain_pics_two_causal_results$SNP1
)

big_env_het_plain_pics_three_causal_results$PS1ContainsCausal <- mapply(
  getIntersect,
  big_env_het_plain_pics_three_causal_results$CausalSNPs,
  big_env_het_plain_pics_three_causal_results$SNP1
)

rec_count_vec <- c(as.numeric(table(big_env_het_plain_pics_one_causal_results$PS1ContainsCausal))[2],
                   as.numeric(table(big_env_het_plain_pics_two_causal_results$PS1ContainsCausal))[2],
                   as.numeric(table(big_env_het_plain_pics_three_causal_results$PS1ContainsCausal))[2])
rec_prob_vec <- rec_count_vec/480
rec_prob_se <- sqrt(rec_prob_vec*(1-rec_prob_vec)/480)
plain_pics_plotting_df <- data.frame(Algorithm=rep("Plain PICS",3),
                                      No_Causal=c(1,2,3),
                                      Rec_Prob=rec_prob_vec,
                                      BinomSE=rec_prob_se)

# Merge tabulated arrays
figA_plotting_df <- rbind(stable_pics_plotting_df,
                              plain_pics_plotting_df)

## Construct plot
figA_plot <- ggplot(figA_plotting_df,aes(x=factor(No_Causal),y=Rec_Prob)) +
  geom_point(aes(colour=Algorithm),
             alpha=0.8)+
  geom_errorbar(aes(ymin=Rec_Prob-1.96*BinomSE,
                    ymax=Rec_Prob+1.96*BinomSE,
                    colour=Algorithm),
                alpha=0.8,
                width=.15)+
  theme_bw() +
  ylim(c(0.5,0.8))+
  ylab("Recovery Frequency") +
  xlab("Number of Causal Variants") +
  ggtitle("A. Plain PICS vs Stable PICS") +
  scale_color_manual(values=c("#f781bf","#2171b5")) +
  labs(colour="Approach") + 
  theme(plot.title=element_text(face="bold"),
        legend.position ='bottom',
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_text(size=13),
        legend.text=element_text(size=11))

## Check intersecting counts
print(getSNRIntersect(big_plain_pics_df,
                      big_stable_pics_df,
                      snp=1))

## Fig B. ======================================================================
## Load all simulation results for stable, top and combined PICS
# Stable PICS2 
all_stable_pics2_results <- do.call(rbind, lapply(22:1, function(x) {
  readr::read_csv(paste0(pics2_res_dir, "PICS2_Chr",x,"_results.csv"))
}))

# Top PICS2
all_top_pics2_results <- do.call(rbind, lapply(22:1, function(x) {
  readr::read_csv(paste0(pics2_top_res_dir, "PICS2_Chr",x,"_top_results.csv"))
}))

# Combined PICS2
all_combined_pics2_results <- do.call(rbind, lapply(c(22:1), function(x) {
  readr::read_csv(paste0(pics2_combined_res_dir, 
                         "PICS2_Chr",x,"_combined_results.csv"))
}))

## Gather recovery frequencies for Potential Set 1
# Top PICS2
all_top_pics2_results$PS1ContainsCausal <- mapply(
  getIntersect,
  all_top_pics2_results$CausalSNPs,
  all_top_pics2_results$SNP1
) 

top_pics2_df_1 <- data.frame(
  SNR=sapply(c(0.05,0.1,0.2,0.4),
             function(x) {return(x/(1-x))}),
  RecProb=as.numeric(
    table(all_top_pics2_results$Phi,
          all_top_pics2_results$PS1ContainsCausal)[,2])/600
)
top_pics2_df_1$BinomSE <- sqrt(top_pics2_df_1$RecProb*
                                 (1-top_pics2_df_1$RecProb)/600)
top_pics2_df_1$Algorithm <- 
  rep("Top PICS", nrow(top_pics2_df_1))

# Stable PICS2
all_stable_pics2_results$PS1ContainsCausal <- mapply(
  getIntersect,
  all_stable_pics2_results$CausalSNPs,
  all_stable_pics2_results$SNP1
) 

stable_pics2_df_1 <- data.frame(
  SNR=sapply(c(0.05,0.1,0.2,0.4),
             function(x) {return(x/(1-x))}),
  RecProb=as.numeric(
    table(all_stable_pics2_results$Phi,
          all_stable_pics2_results$PS1ContainsCausal)[,2])/600
)
stable_pics2_df_1$BinomSE <- sqrt(stable_pics2_df_1$RecProb*
                                    (1-stable_pics2_df_1$RecProb)/600)
stable_pics2_df_1$Algorithm <- 
  rep("Stable PICS", nrow(stable_pics2_df_1))

# Combined PICS2
all_combined_pics2_results$PS1ContainsCausal <- mapply(
  getIntersect,
  all_combined_pics2_results$CausalSNPs,
  all_combined_pics2_results$SNP1
) 

combined_pics2_df_1 <- data.frame(
  SNR=sapply(c(0.05,0.1,0.2,0.4),
             function(x) {return(x/(1-x))}),
  RecProb=as.numeric(
    table(all_combined_pics2_results$Phi,
          all_combined_pics2_results$PS1ContainsCausal)[,2])/600
)
combined_pics2_df_1$BinomSE <- sqrt(combined_pics2_df_1$RecProb*
                                      (1-combined_pics2_df_1$RecProb)/600)
combined_pics2_df_1$Algorithm <- 
  rep("Combined PICS", nrow(combined_pics2_df_1))

## Merge tabulated arrays
merged_df_1 <- rbind(stable_pics2_df_1,
                     top_pics2_df_1,
                     combined_pics2_df_1)

## Construct plot
figB_plot <- ggplot(merged_df_1,aes(x=SNR,y=RecProb)) +
  geom_point(aes(colour=Algorithm),
             alpha=0.8)+
  geom_errorbar(aes(ymin=RecProb-1.96*BinomSE,
                    ymax=RecProb+1.96*BinomSE,
                    colour=Algorithm),
                alpha=0.8,
                width=.05)+
  geom_line(lty="dashed",
            aes(colour=Algorithm),
            alpha=0.8) +
  theme_bw() +
  xlim(c(0.05,0.67)) +
  ylim(c(0,1))+
  scale_x_continuous(breaks=c(0.05,0.1,0.25,0.67)) +
  ylab("Recovery Frequency") +
  xlab("Signal-to-Noise Ratio (SNR)") +
  ggtitle("B. Performance of PICS Algorithms") +
  labs(colour="Approach") +
  scale_color_manual(values=c("#f16913","#2171b5","#238b45")) +
  theme(plot.title=element_text(face="bold"),
        legend.position ='bottom',
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=13),
        legend.text=element_text(size=11))

## Check intersecting counts
print(getSNRIntersect(all_stable_pics2_results,
                      all_top_pics2_results,
                      snp=1))

print(getSNRIntersect(all_stable_pics2_results,
                      all_combined_pics2_results,
                      snp=1))

print(getSNRIntersect(all_combined_pics2_results,
                      all_top_pics2_results,
                      snp=1))

sum(!is.na(mapply(getIntersect3,
                  all_stable_pics2_results$SNP1,
                  all_top_pics2_results$SNP1,
                  all_combined_pics2_results$SNP1)))

## Fig C. ======================================================================
## Load all simulation results for stable and top SuSiE
# SuSiE
all_top_susie_results <- do.call(rbind, lapply(22:1, function(x) {
  readr::read_csv(paste0(susie_res_dir, 
                         "PICS2_Chr",x,"_susieR_L3_results.csv"))
}))

# Stability-guided SuSiE
all_stable_susie_results <- readr::read_csv(
  paste0(stable_susie_res_dir,
         "stable_susie_all_results.csv"))

## Gather recovery frequencies for Potential Set 1
# Top SuSiE
all_top_susie_results$PS1ContainsCausal <- mapply(
  getIntersect,
  all_top_susie_results$CausalSNPs,
  all_top_susie_results$SNP1
) 

top_susie_df_1 <- data.frame(
  SNR=sapply(c(0.05,0.1,0.2,0.4),
             function(x) {return(x/(1-x))}),
  RecProb=as.numeric(
    table(all_top_susie_results$Phi,
          all_top_susie_results$PS1ContainsCausal)[,2])/600
)
top_susie_df_1$BinomSE <- sqrt(top_susie_df_1$RecProb*
                                 (1-top_susie_df_1$RecProb)/600)
top_susie_df_1$Algorithm <- rep("Top SuSiE", nrow(top_susie_df_1))

# Stable SuSiE
all_stable_susie_results$PS1ContainsCausal <- mapply(
  getIntersect,
  all_stable_susie_results$CausalSNPs,
  all_stable_susie_results$SNP1
) 

stable_susie_df_1 <- data.frame(
  SNR=sapply(c(0.05,0.1,0.2,0.4),
             function(x) {return(x/(1-x))}),
  RecProb=as.numeric(
    table(all_stable_susie_results$Phi,
          all_stable_susie_results$PS1ContainsCausal)[,2])/600
)
stable_susie_df_1$BinomSE <- sqrt(stable_susie_df_1$RecProb*
                                    (1-stable_susie_df_1$RecProb)/600)
stable_susie_df_1$Algorithm <- rep("Stable SuSiE", 
                                   nrow(stable_susie_df_1))

## Merge tabulated arrays
merged_susie_df_1 <- rbind(stable_susie_df_1,
                           top_susie_df_1)

## Construct plots
figC_plot <- ggplot(merged_susie_df_1,aes(x=SNR,y=RecProb)) +
  geom_point(aes(colour=Algorithm),
             alpha=0.8)+
  geom_errorbar(aes(ymin=RecProb-1.96*BinomSE,
                    ymax=RecProb+1.96*BinomSE,
                    colour=Algorithm),
                alpha=0.8,
                width=.05)+
  geom_line(lty="dashed",
            aes(colour=Algorithm),
            alpha=0.8) +
  theme_bw() +
  xlim(c(0.05,0.67)) +
  ylim(c(0,1))+
  scale_x_continuous(breaks=c(0.05,0.1,0.25,0.67)) +
  ylab("Recovery Frequency") +
  xlab("Signal-to-Noise Ratio (SNR)") +
  ggtitle("C. Performance of SuSiE Algorithms") +
  scale_color_manual(values=c("#2171b5","#238b45")) +
  labs(colour="Approach") + 
  theme(plot.title=element_text(face="bold"),
        legend.position ='bottom',
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=13),
        legend.text=element_text(size=11))

## Check intersecting counts
print(getSNRIntersect(all_stable_susie_results,
                      all_top_susie_results,
                      snp=1))

## Fig D. ======================================================================
## Gather PICS matching vs non-matching
pics_match_vs_nonmatch <- getMatchRecProbDF2(all_top_pics2_results,
                                             all_stable_pics2_results) %>%
  subset(PS == "Potential Set 1")
pics_match_vs_nonmatch$Algorithm <- rep("PICS", nrow(pics_match_vs_nonmatch))

# To solve jitter problem (improve visualization), manually perturb SNR by +0.02
pics_match_vs_nonmatch$SNR <- pics_match_vs_nonmatch$SNR+0.004

## Gather SuSiE matching vs non-matching
susie_match_vs_nonmatch <- getMatchRecProbDF2(all_top_susie_results,
                                             all_stable_susie_results) %>%
  subset(PS == "Potential Set 1")

susie_match_vs_nonmatch$Algorithm <- rep("SuSiE", nrow(susie_match_vs_nonmatch))

# To solve jitter problem (improve visualization), manually perturb SNR by -0.02
susie_match_vs_nonmatch$SNR <- susie_match_vs_nonmatch$SNR-0.004

## Combine PICS and SuSiE
combined_match_vs_nonmatch <- rbind(pics_match_vs_nonmatch,
                                    susie_match_vs_nonmatch)
## Construct plot
figD_plot <- ggplot(combined_match_vs_nonmatch,
       aes(x=SNR,y=RecProb)) +
  geom_point(aes(colour=Type,
                 shape=Algorithm),
             alpha=0.8)+
  geom_errorbar(aes(ymin=RecProb-1.96*BinomSE,
                    ymax=RecProb+1.96*BinomSE,
                    colour=Type),
                alpha=0.8,
                width=.05)+
  geom_line(aes(colour=Type,
                lty=Algorithm,
                group=interaction(Algorithm,Type)),
            alpha=0.8) +
  theme_bw() +
  xlim(c(0.05,0.67)) +
  ylim(c(0,1))+
  scale_x_continuous(breaks=c(0.05,0.1,0.25,0.67)) +
  ylab("Recovery Frequency") +
  xlab("Signal-to-Noise Ratio (SNR)") +
  ggtitle("D. Matching vs Non-matching Variants") +
  scale_color_manual(values=c("#d7301f","#2171b5","#238b45")) +
  scale_shape_manual(values=c("PICS" = 1, "SuSiE" = 16)) +
  #guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme(plot.title=element_text(face="bold"),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10))

## Combined Figure =============================================================
combined_maintext_fig <- gridExtra::grid.arrange(figA_plot,figB_plot,
                                                 figC_plot,figD_plot,
                                                 nrow=2,ncol=2)
ggsave(combined_maintext_fig,
       file=paste0(plot_dir,"Sim_MainText_Fig.pdf"),
       height=11,
       width=10.5)

#### Supplement Figures --------------------------------------------------------
## Fig 1. Comparison of Plain PICS and Stability-guided PICS (PS2 and PS3) =====
## Fig 1A (Potential Set 2) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Compute whether PS2 Contains at least one causal variant (Stable PICS)
big_env_het_stable_pics_one_causal_results$PS2ContainsCausal <- mapply(
  getIntersect,
  big_env_het_stable_pics_one_causal_results$CausalSNPs,
  big_env_het_stable_pics_one_causal_results$SNP2
) 

big_env_het_stable_pics_two_causal_results$PS2ContainsCausal <- mapply(
  getIntersect,
  big_env_het_stable_pics_two_causal_results$CausalSNPs,
  big_env_het_stable_pics_two_causal_results$SNP2
)

big_env_het_stable_pics_three_causal_results$PS2ContainsCausal <- mapply(
  getIntersect,
  big_env_het_stable_pics_three_causal_results$CausalSNPs,
  big_env_het_stable_pics_three_causal_results$SNP2
)

rec_count_vec <- c(as.numeric(table(big_env_het_stable_pics_one_causal_results$PS2ContainsCausal))[2],
                   as.numeric(table(big_env_het_stable_pics_two_causal_results$PS2ContainsCausal))[2],
                   as.numeric(table(big_env_het_stable_pics_three_causal_results$PS2ContainsCausal))[2])
rec_prob_vec <- rec_count_vec/480
rec_prob_se <- sqrt(rec_prob_vec*(1-rec_prob_vec)/480)
stable_pics_plotting_df <- data.frame(Algorithm=rep("Stable PICS",3),
                                      No_Causal=c(1,2,3),
                                      Rec_Prob=rec_prob_vec,
                                      BinomSE=rec_prob_se)

# Compute whether PS2 Contains at least one causal variant (Plain PICS)
big_env_het_plain_pics_one_causal_results$PS2ContainsCausal <- mapply(
  getIntersect,
  big_env_het_plain_pics_one_causal_results$CausalSNPs,
  big_env_het_plain_pics_one_causal_results$SNP2
) 

big_env_het_plain_pics_two_causal_results$PS2ContainsCausal <- mapply(
  getIntersect,
  big_env_het_plain_pics_two_causal_results$CausalSNPs,
  big_env_het_plain_pics_two_causal_results$SNP2
)

big_env_het_plain_pics_three_causal_results$PS2ContainsCausal <- mapply(
  getIntersect,
  big_env_het_plain_pics_three_causal_results$CausalSNPs,
  big_env_het_plain_pics_three_causal_results$SNP2
)

rec_count_vec <- c(as.numeric(table(big_env_het_plain_pics_one_causal_results$PS2ContainsCausal))[2],
                   as.numeric(table(big_env_het_plain_pics_two_causal_results$PS2ContainsCausal))[2],
                   as.numeric(table(big_env_het_plain_pics_three_causal_results$PS2ContainsCausal))[2])
rec_prob_vec <- rec_count_vec/480
rec_prob_se <- sqrt(rec_prob_vec*(1-rec_prob_vec)/480)
plain_pics_plotting_df <- data.frame(Algorithm=rep("Plain PICS",3),
                                     No_Causal=c(1,2,3),
                                     Rec_Prob=rec_prob_vec,
                                     BinomSE=rec_prob_se)

# Merge tabulated arrays
supp_fig1A_plotting_df <- rbind(stable_pics_plotting_df,
                          plain_pics_plotting_df)

## Construct plot
supp_fig1A_plot <- ggplot(supp_fig1A_plotting_df,
                          aes(x=factor(No_Causal),y=Rec_Prob)) +
  geom_point(aes(colour=Algorithm),
             alpha=0.8,
             size=2.5)+
  geom_errorbar(aes(ymin=Rec_Prob-1.96*BinomSE,
                    ymax=Rec_Prob+1.96*BinomSE,
                    colour=Algorithm),
                alpha=0.8,
                width=.2)+
  theme_bw() +
  ylim(c(0,0.3))+
  ylab("Recovery Frequency") +
  xlab("Number of Causal Variants") +
  ggtitle("A. Plain PICS vs Stable PICS (Potential Set 2)") +
  scale_color_manual(values=c("#f781bf","#2171b5")) +
  labs(colour="Approach") + 
  theme(plot.title=element_text(face="bold",size=16),
        legend.position ='None',
        axis.text=element_text(size=13),
        axis.title=element_text(size=15),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))

## Check intersecting counts
print(getSNRIntersect(big_plain_pics_df,
                      big_stable_pics_df,
                      snp=2))
# 1088 vs 352

## Fig 1B (Potential Set 3) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Compute whether PS3 Contains at least one causal variant (Stable PICS)
big_env_het_stable_pics_one_causal_results$PS3ContainsCausal <- mapply(
  getIntersect,
  big_env_het_stable_pics_one_causal_results$CausalSNPs,
  big_env_het_stable_pics_one_causal_results$SNP3
) 

big_env_het_stable_pics_two_causal_results$PS3ContainsCausal <- mapply(
  getIntersect,
  big_env_het_stable_pics_two_causal_results$CausalSNPs,
  big_env_het_stable_pics_two_causal_results$SNP3
)

big_env_het_stable_pics_three_causal_results$PS3ContainsCausal <- mapply(
  getIntersect,
  big_env_het_stable_pics_three_causal_results$CausalSNPs,
  big_env_het_stable_pics_three_causal_results$SNP3
)

rec_count_vec <- c(as.numeric(table(big_env_het_stable_pics_one_causal_results$PS3ContainsCausal))[2],
                   as.numeric(table(big_env_het_stable_pics_two_causal_results$PS3ContainsCausal))[2],
                   as.numeric(table(big_env_het_stable_pics_three_causal_results$PS3ContainsCausal))[2])
rec_prob_vec <- rec_count_vec/480
rec_prob_se <- sqrt(rec_prob_vec*(1-rec_prob_vec)/480)
stable_pics_plotting_df <- data.frame(Algorithm=rep("Stable PICS",3),
                                      No_Causal=c(1,2,3),
                                      Rec_Prob=rec_prob_vec,
                                      BinomSE=rec_prob_se)

# Compute whether PS3 Contains at least one causal variant (Plain PICS)
big_env_het_plain_pics_one_causal_results$PS3ContainsCausal <- mapply(
  getIntersect,
  big_env_het_plain_pics_one_causal_results$CausalSNPs,
  big_env_het_plain_pics_one_causal_results$SNP3
) 

big_env_het_plain_pics_two_causal_results$PS3ContainsCausal <- mapply(
  getIntersect,
  big_env_het_plain_pics_two_causal_results$CausalSNPs,
  big_env_het_plain_pics_two_causal_results$SNP3
)

big_env_het_plain_pics_three_causal_results$PS3ContainsCausal <- mapply(
  getIntersect,
  big_env_het_plain_pics_three_causal_results$CausalSNPs,
  big_env_het_plain_pics_three_causal_results$SNP3
)

rec_count_vec <- c(as.numeric(table(big_env_het_plain_pics_one_causal_results$PS3ContainsCausal))[2],
                   as.numeric(table(big_env_het_plain_pics_two_causal_results$PS3ContainsCausal))[2],
                   as.numeric(table(big_env_het_plain_pics_three_causal_results$PS3ContainsCausal))[2])
rec_prob_vec <- rec_count_vec/480
rec_prob_se <- sqrt(rec_prob_vec*(1-rec_prob_vec)/480)
plain_pics_plotting_df <- data.frame(Algorithm=rep("Plain PICS",3),
                                     No_Causal=c(1,2,3),
                                     Rec_Prob=rec_prob_vec,
                                     BinomSE=rec_prob_se)

# Merge tabulated arrays
supp_fig1B_plotting_df <- rbind(stable_pics_plotting_df,
                                plain_pics_plotting_df)

## Construct plot
supp_fig1B_plot <- ggplot(supp_fig1B_plotting_df,
                          aes(x=factor(No_Causal),y=Rec_Prob)) +
  geom_point(aes(colour=Algorithm),
             alpha=0.8,
             size=2.5)+
  geom_errorbar(aes(ymin=Rec_Prob-1.96*BinomSE,
                    ymax=Rec_Prob+1.96*BinomSE,
                    colour=Algorithm),
                alpha=0.8,
                width=.2)+
  theme_bw() +
  #ylim(c(0.5,0.8))+
  ylab("Recovery Frequency") +
  xlab("Number of Causal Variants") +
  ggtitle("B. Plain PICS vs Stable PICS (Potential Set 3)") +
  scale_color_manual(values=c("#f781bf","#2171b5")) +
  labs(colour="Approach") + 
  theme(plot.title=element_text(face="bold",size=16),
        legend.position ='bottom',
        axis.text=element_text(size=13),
        axis.title=element_text(size=15),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))

## Check intersecting counts
print(getSNRIntersect(big_plain_pics_df,
                      big_stable_pics_df,
                      snp=3))
# 1090 vs 350 

## Combine subfigures A and B <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
supp_fig1_plot <- gridExtra::grid.arrange(supp_fig1A_plot,
                                          supp_fig1B_plot,
                                          nrow=2,
                                          heights = c(2.4, 2.8))
ggsave(supp_fig1_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig1.pdf"),
       height=11.8,
       width=5.9)

## Fig 2. PICS SNR vs Recovery plot (PS2 and PS3) ==============================
## Fig 2A (Potential Set 2) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
## Gather recovery frequencies for Potential Set 2
# Top PICS2
all_top_pics2_results$PS2ContainsCausal <- mapply(
  getIntersect,
  all_top_pics2_results$CausalSNPs,
  all_top_pics2_results$SNP2
) 

top_pics2_df_1 <- data.frame(
  SNR=sapply(c(0.05,0.1,0.2,0.4),
             function(x) {return(x/(1-x))}),
  RecProb=as.numeric(
    table(all_top_pics2_results$Phi,
          all_top_pics2_results$PS2ContainsCausal)[,2])/600
)
top_pics2_df_1$BinomSE <- sqrt(top_pics2_df_1$RecProb*
                                 (1-top_pics2_df_1$RecProb)/600)
top_pics2_df_1$Algorithm <- 
  rep("Top PICS", nrow(top_pics2_df_1))

# Stable PICS2
all_stable_pics2_results$PS2ContainsCausal <- mapply(
  getIntersect,
  all_stable_pics2_results$CausalSNPs,
  all_stable_pics2_results$SNP2
) 

stable_pics2_df_1 <- data.frame(
  SNR=sapply(c(0.05,0.1,0.2,0.4),
             function(x) {return(x/(1-x))}),
  RecProb=as.numeric(
    table(all_stable_pics2_results$Phi,
          all_stable_pics2_results$PS2ContainsCausal)[,2])/600
)
stable_pics2_df_1$BinomSE <- sqrt(stable_pics2_df_1$RecProb*
                                    (1-stable_pics2_df_1$RecProb)/600)
stable_pics2_df_1$Algorithm <- 
  rep("Stable PICS", nrow(stable_pics2_df_1))

# Combined PICS2
all_combined_pics2_results$PS2ContainsCausal <- mapply(
  getIntersect,
  all_combined_pics2_results$CausalSNPs,
  all_combined_pics2_results$SNP2
) 

combined_pics2_df_1 <- data.frame(
  SNR=sapply(c(0.05,0.1,0.2,0.4),
             function(x) {return(x/(1-x))}),
  RecProb=as.numeric(
    table(all_combined_pics2_results$Phi,
          all_combined_pics2_results$PS2ContainsCausal)[,2])/600
)
combined_pics2_df_1$BinomSE <- sqrt(combined_pics2_df_1$RecProb*
                                      (1-combined_pics2_df_1$RecProb)/600)
combined_pics2_df_1$Algorithm <- 
  rep("Combined PICS", nrow(combined_pics2_df_1))

## Merge tabulated arrays
merged_df_ps2<- rbind(stable_pics2_df_1,
                     top_pics2_df_1,
                     combined_pics2_df_1)

## Construct plot
supp_fig2A_plot <- ggplot(merged_df_ps2,aes(x=SNR,y=RecProb)) +
  geom_point(aes(colour=Algorithm),
             alpha=0.8)+
  geom_errorbar(aes(ymin=RecProb-1.96*BinomSE,
                    ymax=RecProb+1.96*BinomSE,
                    colour=Algorithm),
                alpha=0.8,
                width=.05)+
  geom_line(lty="dashed",
            aes(colour=Algorithm),
            alpha=0.8) +
  theme_bw() +
  xlim(c(0.05,0.67)) +
  ylim(c(0,0.25))+
  scale_x_continuous(breaks=c(0.05,0.1,0.25,0.67)) +
  ylab("Recovery Frequency") +
  xlab("Signal-to-Noise Ratio (SNR)") +
  ggtitle("A. Performance of PICS Algorithms\n     (Potential Set 2)") +
  labs(colour="Approach") +
  scale_color_manual(values=c("#f16913","#2171b5","#238b45")) +
  theme(plot.title=element_text(face="bold",size=16),
        legend.position ='None',
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=15),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))

## Check intersecting counts
print(getSNRIntersect(all_stable_pics2_results,
                      all_top_pics2_results,
                      snp=2))

print(getSNRIntersect(all_stable_pics2_results,
                      all_combined_pics2_results,
                      snp=2))

print(getSNRIntersect(all_combined_pics2_results,
                      all_top_pics2_results,
                      snp=2))

sum(!is.na(mapply(getIntersect3,
                  all_stable_pics2_results$SNP2,
                  all_top_pics2_results$SNP2,
                  all_combined_pics2_results$SNP2)))

## Fig 2B (Potential Set 3) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
## Gather recovery frequencies for Potential Set 2
# Top PICS2
all_top_pics2_results$PS3ContainsCausal <- mapply(
  getIntersect,
  all_top_pics2_results$CausalSNPs,
  all_top_pics2_results$SNP3
) 

top_pics2_df_1 <- data.frame(
  SNR=sapply(c(0.05,0.1,0.2,0.4),
             function(x) {return(x/(1-x))}),
  RecProb=as.numeric(
    table(all_top_pics2_results$Phi,
          all_top_pics2_results$PS3ContainsCausal)[,2])/600
)
top_pics2_df_1$BinomSE <- sqrt(top_pics2_df_1$RecProb*
                                 (1-top_pics2_df_1$RecProb)/600)
top_pics2_df_1$Algorithm <- 
  rep("Top PICS", nrow(top_pics2_df_1))

# Stable PICS2
all_stable_pics2_results$PS3ContainsCausal <- mapply(
  getIntersect,
  all_stable_pics2_results$CausalSNPs,
  all_stable_pics2_results$SNP3
) 

stable_pics2_df_1 <- data.frame(
  SNR=sapply(c(0.05,0.1,0.2,0.4),
             function(x) {return(x/(1-x))}),
  RecProb=as.numeric(
    table(all_stable_pics2_results$Phi,
          all_stable_pics2_results$PS3ContainsCausal)[,2])/600
)
stable_pics2_df_1$BinomSE <- sqrt(stable_pics2_df_1$RecProb*
                                    (1-stable_pics2_df_1$RecProb)/600)
stable_pics2_df_1$Algorithm <- 
  rep("Stable PICS", nrow(stable_pics2_df_1))

# Combined PICS2
all_combined_pics2_results$PS3ContainsCausal <- mapply(
  getIntersect,
  all_combined_pics2_results$CausalSNPs,
  all_combined_pics2_results$SNP3
) 

combined_pics2_df_1 <- data.frame(
  SNR=sapply(c(0.05,0.1,0.2,0.4),
             function(x) {return(x/(1-x))}),
  RecProb=as.numeric(
    table(all_combined_pics2_results$Phi,
          all_combined_pics2_results$PS3ContainsCausal)[,2])/600
)
combined_pics2_df_1$BinomSE <- sqrt(combined_pics2_df_1$RecProb*
                                      (1-combined_pics2_df_1$RecProb)/600)
combined_pics2_df_1$Algorithm <- 
  rep("Combined PICS", nrow(combined_pics2_df_1))

## Merge tabulated arrays
merged_df_ps3<- rbind(stable_pics2_df_1,
                      top_pics2_df_1,
                      combined_pics2_df_1)

## Construct plot
supp_fig2B_plot <- ggplot(merged_df_ps3,aes(x=SNR,y=RecProb)) +
  geom_point(aes(colour=Algorithm),
             alpha=0.8)+
  geom_errorbar(aes(ymin=RecProb-1.96*BinomSE,
                    ymax=RecProb+1.96*BinomSE,
                    colour=Algorithm),
                alpha=0.8,
                width=.05)+
  geom_line(lty="dashed",
            aes(colour=Algorithm),
            alpha=0.8) +
  theme_bw() +
  xlim(c(0.05,0.67)) +
  ylim(c(0,0.12))+
  scale_x_continuous(breaks=c(0.05,0.1,0.25,0.67)) +
  ylab("Recovery Frequency") +
  xlab("Signal-to-Noise Ratio (SNR)") +
  ggtitle("B. Performance of PICS Algorithms\n     (Potential Set 3)") +
  labs(colour="Approach") +
  scale_color_manual(values=c("#f16913","#2171b5","#238b45")) +
  theme(plot.title=element_text(face="bold",size=16),
        legend.position ='bottom',
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=15),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))

## Check intersecting counts
print(getSNRIntersect(all_stable_pics2_results,
                      all_top_pics2_results,
                      snp=3))

print(getSNRIntersect(all_stable_pics2_results,
                      all_combined_pics2_results,
                      snp=3))

print(getSNRIntersect(all_combined_pics2_results,
                      all_top_pics2_results,
                      snp=3))

sum(!is.na(mapply(getIntersect3,
                  all_stable_pics2_results$SNP3,
                  all_top_pics2_results$SNP3,
                  all_combined_pics2_results$SNP3)))

## Combine subfigures A and B <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
supp_fig2_plot <- gridExtra::grid.arrange(supp_fig2A_plot,
                                          supp_fig2B_plot,
                                          nrow=2,
                                          heights = c(2.4, 2.8))
ggsave(supp_fig2_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig2.pdf"),
       height=11.8,
       width=5.9)

## Fig 3. SuSiE SNR vs Recovery plot (CS2 and CS3) =============================
## Fig 3A (Credible Set 2) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
## Gather recovery frequencies for Potential Set 2
# Top SuSiE
all_top_susie_results$PS2ContainsCausal <- mapply(
  getIntersect,
  all_top_susie_results$CausalSNPs,
  all_top_susie_results$SNP2
) 

top_susie_df_1 <- data.frame(
  SNR=sapply(c(0.05,0.1,0.2,0.4),
             function(x) {return(x/(1-x))}),
  RecProb=as.numeric(
    table(all_top_susie_results$Phi,
          all_top_susie_results$PS2ContainsCausal)[,2])/600
)
top_susie_df_1$BinomSE <- sqrt(top_susie_df_1$RecProb*
                                 (1-top_susie_df_1$RecProb)/600)
top_susie_df_1$Algorithm <- rep("Top SuSiE", nrow(top_susie_df_1))

# Stable SuSiE
all_stable_susie_results$PS2ContainsCausal <- mapply(
  getIntersect,
  all_stable_susie_results$CausalSNPs,
  all_stable_susie_results$SNP2
) 

stable_susie_df_1 <- data.frame(
  SNR=sapply(c(0.05,0.1,0.2,0.4),
             function(x) {return(x/(1-x))}),
  RecProb=as.numeric(
    table(all_stable_susie_results$Phi,
          all_stable_susie_results$PS2ContainsCausal)[,2])/600
)
stable_susie_df_1$BinomSE <- sqrt(stable_susie_df_1$RecProb*
                                    (1-stable_susie_df_1$RecProb)/600)
stable_susie_df_1$Algorithm <- rep("Stable SuSiE", 
                                   nrow(stable_susie_df_1))

## Merge tabulated arrays
merged_susie_df_ps2 <- rbind(stable_susie_df_1,
                           top_susie_df_1)

## Construct plots
supp_fig3A_plot <- ggplot(merged_susie_df_ps2,aes(x=SNR,y=RecProb)) +
  geom_point(aes(colour=Algorithm),
             alpha=0.8)+
  geom_errorbar(aes(ymin=RecProb-1.96*BinomSE,
                    ymax=RecProb+1.96*BinomSE,
                    colour=Algorithm),
                alpha=0.8,
                width=.05)+
  geom_line(lty="dashed",
            aes(colour=Algorithm),
            alpha=0.8) +
  theme_bw() +
  xlim(c(0.05,0.67)) +
  ylim(c(0,0.5))+
  scale_x_continuous(breaks=c(0.05,0.1,0.25,0.67)) +
  ylab("Recovery Frequency") +
  xlab("Signal-to-Noise Ratio (SNR)") +
  ggtitle("A. Performance of SuSiE Algorithms\n     (Credible Set 2)") +
  scale_color_manual(values=c("#2171b5","#238b45")) +
  labs(colour="Approach") + 
  theme(plot.title=element_text(face="bold",size=16),
        legend.position ='None',
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=15),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))

## Check intersecting counts
print(getSNRIntersect(all_stable_susie_results,
                      all_top_susie_results,
                      snp=2))
# 558 vs 1842

## Fig 3B (Credible Set 3) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
## Gather recovery frequencies for Potential Set 3
# Top SuSiE
all_top_susie_results$PS3ContainsCausal <- mapply(
  getIntersect,
  all_top_susie_results$CausalSNPs,
  all_top_susie_results$SNP3
) 

top_susie_df_1 <- data.frame(
  SNR=sapply(c(0.05,0.1,0.2,0.4),
             function(x) {return(x/(1-x))}),
  RecProb=as.numeric(
    table(all_top_susie_results$Phi,
          all_top_susie_results$PS3ContainsCausal)[,2])/600
)
top_susie_df_1$BinomSE <- sqrt(top_susie_df_1$RecProb*
                                 (1-top_susie_df_1$RecProb)/600)
top_susie_df_1$Algorithm <- rep("Top SuSiE", nrow(top_susie_df_1))

# Stable SuSiE
all_stable_susie_results$PS3ContainsCausal <- mapply(
  getIntersect,
  all_stable_susie_results$CausalSNPs,
  all_stable_susie_results$SNP3
) 

stable_susie_df_1 <- data.frame(
  SNR=sapply(c(0.05,0.1,0.2,0.4),
             function(x) {return(x/(1-x))}),
  RecProb=as.numeric(
    table(all_stable_susie_results$Phi,
          all_stable_susie_results$PS3ContainsCausal)[,2])/600
)
stable_susie_df_1$BinomSE <- sqrt(stable_susie_df_1$RecProb*
                                    (1-stable_susie_df_1$RecProb)/600)
stable_susie_df_1$Algorithm <- rep("Stable SuSiE", 
                                   nrow(stable_susie_df_1))

## Merge tabulated arrays
merged_susie_df_ps3 <- rbind(stable_susie_df_1,
                             top_susie_df_1)

## Construct plots
supp_fig3B_plot <- ggplot(merged_susie_df_ps3,aes(x=SNR,y=RecProb)) +
  geom_point(aes(colour=Algorithm),
             alpha=0.8)+
  geom_errorbar(aes(ymin=RecProb-1.96*BinomSE,
                    ymax=RecProb+1.96*BinomSE,
                    colour=Algorithm),
                alpha=0.8,
                width=.05)+
  geom_line(lty="dashed",
            aes(colour=Algorithm),
            alpha=0.8) +
  theme_bw() +
  xlim(c(0.05,0.67)) +
  ylim(c(0,0.5))+
  scale_x_continuous(breaks=c(0.05,0.1,0.25,0.67)) +
  ylab("Recovery Frequency") +
  xlab("Signal-to-Noise Ratio (SNR)") +
  ggtitle("B. Performance of SuSiE Algorithms\n     (Credible Set 3)") +
  scale_color_manual(values=c("#2171b5","#238b45")) +
  labs(colour="Approach") + 
  theme(plot.title=element_text(face="bold",size=16),
        legend.position ='bottom',
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=15),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))

## Check intersecting counts
print(getSNRIntersect(all_stable_susie_results,
                      all_top_susie_results,
                      snp=3))

# 319 vs 2081

## Combine subfigures A and B <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
supp_fig3_plot <- gridExtra::grid.arrange(supp_fig3A_plot,
                                          supp_fig3B_plot,
                                          nrow=2,
                                          heights = c(2.4, 2.8))
ggsave(supp_fig3_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig3.pdf"),
       height=11.8,
       width=5.9)

## Fig 4. Matching vs Non-matching causal variant recovery (Sets 2 and 3) ======
## Fig 4A (Set 2) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
## Gather PICS matching vs non-matching
pics_match_vs_nonmatch <- getMatchRecProbDF2(all_top_pics2_results,
                                             all_stable_pics2_results) %>%
  subset(PS == "Potential Set 2")
pics_match_vs_nonmatch$Algorithm <- rep("PICS", nrow(pics_match_vs_nonmatch))

# To solve jitter problem (improve visualization), manually perturb SNR by +0.02
pics_match_vs_nonmatch$SNR <- pics_match_vs_nonmatch$SNR+0.004

## Gather SuSiE matching vs non-matching
susie_match_vs_nonmatch <- getMatchRecProbDF2(all_top_susie_results,
                                              all_stable_susie_results) %>%
  subset(PS == "Potential Set 2")

susie_match_vs_nonmatch$Algorithm <- rep("SuSiE", nrow(susie_match_vs_nonmatch))

# To solve jitter problem (improve visualization), manually perturb SNR by -0.02
susie_match_vs_nonmatch$SNR <- susie_match_vs_nonmatch$SNR-0.004

## Combine PICS and SuSiE
combined_match_vs_nonmatch <- rbind(pics_match_vs_nonmatch,
                                    susie_match_vs_nonmatch)
## Construct plot
supp_fig4A_plot <- ggplot(combined_match_vs_nonmatch,
                    aes(x=SNR,y=RecProb)) +
  geom_point(aes(colour=Type,
                 shape=Algorithm),
             alpha=0.8)+
  geom_errorbar(aes(ymin=RecProb-1.96*BinomSE,
                    ymax=RecProb+1.96*BinomSE,
                    colour=Type),
                alpha=0.8,
                width=.05)+
  geom_line(aes(colour=Type,
                lty=Algorithm,
                group=interaction(Algorithm,Type)),
            alpha=0.8) +
  theme_bw() +
  xlim(c(0.05,0.67)) +
  ylim(c(0,0.75))+
  scale_x_continuous(breaks=c(0.05,0.1,0.25,0.67)) +
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.75)) +
  ylab("Recovery Frequency") +
  xlab("Signal-to-Noise Ratio (SNR)") +
  ggtitle("A. Matching vs Non-matching Variants\n     (Credible / Potential Set 2)") +
  scale_color_manual(values=c("#d7301f","#2171b5","#238b45")) +
  scale_shape_manual(values=c("PICS" = 1, "SuSiE" = 16)) +
  #guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme(plot.title=element_text(face="bold"),
        legend.position ='None',
        #legend.box = "vertical",
        #legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=15),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

## Fig 4B (Set 3) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
## Gather PICS matching vs non-matching
pics_match_vs_nonmatch <- getMatchRecProbDF2(all_top_pics2_results,
                                             all_stable_pics2_results) %>%
  subset(PS == "Potential Set 3")
pics_match_vs_nonmatch$Algorithm <- rep("PICS", nrow(pics_match_vs_nonmatch))

# To solve jitter problem (improve visualization), manually perturb SNR by +0.02
pics_match_vs_nonmatch$SNR <- pics_match_vs_nonmatch$SNR+0.004

## Gather SuSiE matching vs non-matching
susie_match_vs_nonmatch <- getMatchRecProbDF2(all_top_susie_results,
                                              all_stable_susie_results) %>%
  subset(PS == "Potential Set 3")

susie_match_vs_nonmatch$Algorithm <- rep("SuSiE", nrow(susie_match_vs_nonmatch))

# To solve jitter problem (improve visualization), manually perturb SNR by -0.02
susie_match_vs_nonmatch$SNR <- susie_match_vs_nonmatch$SNR-0.004

## Combine PICS and SuSiE
combined_match_vs_nonmatch <- rbind(pics_match_vs_nonmatch,
                                    susie_match_vs_nonmatch)
## Construct plot
supp_fig4B_plot <- ggplot(combined_match_vs_nonmatch,
                          aes(x=SNR,y=RecProb)) +
  geom_point(aes(colour=Type,
                 shape=Algorithm),
             alpha=0.8)+
  geom_errorbar(aes(ymin=RecProb-1.96*BinomSE,
                    ymax=RecProb+1.96*BinomSE,
                    colour=Type),
                alpha=0.8,
                width=.05)+
  geom_line(aes(colour=Type,
                lty=Algorithm,
                group=interaction(Algorithm,Type)),
            alpha=0.8) +
  theme_bw() +
  xlim(c(0.05,0.67)) +
  ylim(c(0,0.6))+
  scale_x_continuous(breaks=c(0.05,0.1,0.25,0.67)) +
  ylab("Recovery Frequency") +
  xlab("Signal-to-Noise Ratio (SNR)") +
  ggtitle("B. Matching vs Non-matching Variants\n     (Credible / Potential Set 3)") +
  scale_color_manual(values=c("#d7301f","#2171b5","#238b45")) +
  scale_shape_manual(values=c("PICS" = 1, "SuSiE" = 16)) +
  #guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme(plot.title=element_text(face="bold"),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=15),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

## Combine subfigures A and B <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
supp_fig4_plot <- gridExtra::grid.arrange(supp_fig4A_plot,
                                          supp_fig4B_plot,
                                          nrow=2,
                                          heights = c(2.4,2.8))
ggsave(supp_fig4_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig4.pdf"),
       height=12,
       width=5.9)


## Figs 5-7. Sensitivity of stable fine-mapping to number of potential sets ====
# Plot 3 x 2 figure, where we analyze simulations stratified by the number of 
# causal variants. For each case, inspect the causal recovery rate / distribution
# of no. causal variants recovered as a function of the number of potential sets
# included. Split by SNR.

## Stable PICS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Get causal variant recovery counts 
all_stable_pics2_results$NoCausalSNPsRecovered_1PS <- 
  mapply(countCausalSNPs,
         all_stable_pics2_results$CausalSNPs,
         all_stable_pics2_results$SNP1,
         rep(NA,nrow(all_stable_pics2_results)),
         rep(NA,nrow(all_stable_pics2_results)))
all_stable_pics2_results$NoCausalSNPsRecovered_2PS <- 
  mapply(countCausalSNPs,
         all_stable_pics2_results$CausalSNPs,
         all_stable_pics2_results$SNP1,
         all_stable_pics2_results$SNP2,
         rep(NA,nrow(all_stable_pics2_results)))
all_stable_pics2_results$NoCausalSNPsRecovered_3PS <- 
  mapply(countCausalSNPs,
         all_stable_pics2_results$CausalSNPs,
         all_stable_pics2_results$SNP1,
         all_stable_pics2_results$SNP2,
         all_stable_pics2_results$SNP3)

# Split results by number of causal variants simulated
stable_pics2_one_causal_results <- all_stable_pics2_results %>% 
  subset(S==1)
stable_pics2_two_causal_results <- all_stable_pics2_results %>% 
  subset(S==2)
stable_pics2_three_causal_results <- all_stable_pics2_results %>% 
  subset(S==3)

phi_vector <- c(0.05,0.1,0.2,0.4)

stable_pics2_trend_df <- NULL 
# for SNR ratio parameter
for (i in 1:4) {
  # for no. sets included
  for (j in 1:3) {
    # one causal
    one_causal <- table(
      stable_pics2_one_causal_results$Phi,
      stable_pics2_one_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(one_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 1",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    stable_pics2_trend_df <- rbind(stable_pics2_trend_df, ci_table)
    
    # two causal
    two_causal <- table(
      stable_pics2_two_causal_results$Phi,
      stable_pics2_two_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(two_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 2",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    stable_pics2_trend_df <- rbind(stable_pics2_trend_df, ci_table)
    
    # three causal
    three_causal <- table(
      stable_pics2_three_causal_results$Phi,
      stable_pics2_three_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(three_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 3",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    stable_pics2_trend_df <- rbind(stable_pics2_trend_df, ci_table)
  }
}

stable_pics2_trend_df$Algorithm <- rep("Stable PICS", nrow(stable_pics2_trend_df))

## Stable SuSiE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Get causal variant recovery counts 
all_stable_susie_results$NoCausalSNPsRecovered_1PS <- 
  mapply(countCausalSNPs,
         all_stable_susie_results$CausalSNPs,
         all_stable_susie_results$SNP1,
         rep(NA,nrow(all_stable_susie_results)),
         rep(NA,nrow(all_stable_susie_results)))
all_stable_susie_results$NoCausalSNPsRecovered_2PS <- 
  mapply(countCausalSNPs,
         all_stable_susie_results$CausalSNPs,
         all_stable_susie_results$SNP1,
         all_stable_susie_results$SNP2,
         rep(NA,nrow(all_stable_susie_results)))
all_stable_susie_results$NoCausalSNPsRecovered_3PS <- 
  mapply(countCausalSNPs,
         all_stable_susie_results$CausalSNPs,
         all_stable_susie_results$SNP1,
         all_stable_susie_results$SNP2,
         all_stable_susie_results$SNP3)

# Split results by number of causal variants simulated
stable_susie_one_causal_results <- all_stable_susie_results %>% 
  subset(S==1)
stable_susie_two_causal_results <- all_stable_susie_results %>% 
  subset(S==2)
stable_susie_three_causal_results <- all_stable_susie_results %>% 
  subset(S==3)

phi_vector <- c(0.05,0.1,0.2,0.4)

stable_susie_trend_df <- NULL
# for SNR ratio parameter
for (i in 1:4) {
  # for no. sets included
  for (j in 1:3) {
    # one causal
    one_causal <- table(
      stable_susie_one_causal_results$Phi,
      stable_susie_one_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(one_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 1",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    stable_susie_trend_df <- rbind(stable_susie_trend_df, ci_table)
    
    # two causal
    two_causal <- table(
      stable_susie_two_causal_results$Phi,
      stable_susie_two_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(two_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 2",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    stable_susie_trend_df <- rbind(stable_susie_trend_df, ci_table)
    
    # three causal
    three_causal <- table(
      stable_susie_three_causal_results$Phi,
      stable_susie_three_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(three_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 3",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    stable_susie_trend_df <- rbind(stable_susie_trend_df, ci_table)
  }
}

stable_susie_trend_df$Algorithm <- rep("Stable SuSiE", nrow(stable_susie_trend_df))

## Combine PICS and SuSiE trends <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
stable_trend_df <- rbind(stable_susie_trend_df,
                         stable_pics2_trend_df)
stable_trend_df$SNR <- paste0("SNR = ", round(stable_trend_df$Phi/(1-stable_trend_df$Phi),3))
stable_trend_1c_df <- stable_trend_df %>% 
  subset(NoCausalVarSim == "No. Causal Var = 1")
stable_trend_2c_df <- stable_trend_df %>% 
  subset(NoCausalVarSim == "No. Causal Var = 2")
stable_trend_3c_df <- stable_trend_df %>% 
  subset(NoCausalVarSim == "No. Causal Var = 3")

## Generate plots <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# One causal variant
supp_fig5_plot <- ggplot(stable_trend_1c_df, aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoSetsInclude~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  ggtitle("Simulations with 1 Causal Variant\n\n") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig5_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig5.pdf"),
       height=9,
       width=9)

# Two causal variants
supp_fig6_plot <- ggplot(stable_trend_2c_df, aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoSetsInclude~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  ggtitle("Simulations with 2 Causal Variants\n\n") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig6_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig6.pdf"),
       height=9,
       width=9)

# Three causal variants
supp_fig7_plot <- ggplot(stable_trend_3c_df, aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoSetsInclude~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  ggtitle("Simulations with 3 Causal Variants\n\n") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))
  
ggsave(supp_fig7_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig7.pdf"),
       height=9,
       width=9)

## Fig 8. Difference in posterior probability between non-matching variants ====
# Show the heatmaps, colored by no. causal variants simulated. 
## Potential Set 1 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_one_causal_results,
                                       other_df=top_pics2_one_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_one_causal_results,
                                    other_df=stable_pics2_one_causal_results)

post_prob_df_1c <- data.frame(Stable=stable_pics2_match_df[
  !stable_pics2_match_df$focal_SNP1_match_other_SNP1,]$PostProb1,
  Top=top_pics2_match_df[
    !stable_pics2_match_df$focal_SNP1_match_other_SNP1,]$PostProb1)
post_prob_df_1c$NoCausalVariants <- rep(1,nrow(post_prob_df_1c))

# Two Causal Variants
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_two_causal_results,
                                       other_df=top_pics2_two_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_two_causal_results,
                                    other_df=stable_pics2_two_causal_results)

post_prob_df_2c <- data.frame(Stable=stable_pics2_match_df[
  !stable_pics2_match_df$focal_SNP1_match_other_SNP1,]$PostProb1,
  Top=top_pics2_match_df[
    !stable_pics2_match_df$focal_SNP1_match_other_SNP1,]$PostProb1)
post_prob_df_2c$NoCausalVariants <- rep(2,nrow(post_prob_df_2c))

# Three Causal Variants
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_three_causal_results,
                                       other_df=top_pics2_three_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_three_causal_results,
                                    other_df=stable_pics2_three_causal_results)

post_prob_df_3c <- data.frame(Stable=stable_pics2_match_df[
  !stable_pics2_match_df$focal_SNP1_match_other_SNP1,]$PostProb1,
  Top=top_pics2_match_df[
    !stable_pics2_match_df$focal_SNP1_match_other_SNP1,]$PostProb1)
post_prob_df_3c$NoCausalVariants <- rep(3,nrow(post_prob_df_3c))

post_prob_df <- rbind(post_prob_df_1c,
                      post_prob_df_2c,
                      post_prob_df_3c)

supp_fig8A <- ggplot(post_prob_df,aes(x=Top,y=Stable)) +
  geom_point(aes(colour=factor(NoCausalVariants)),
             alpha=0.4) +
  geom_abline(slope=1, intercept=0,
              lty="dashed") +
  theme_bw() +
  ylim(c(0,1)) +
  xlim(c(0,1)) +
  xlab("Top") +
  ylab("Stable") +
  ggtitle("Potential Set 1") +
  labs(colour="Number of Causal Variants") +
  theme(plot.title=element_text(hjust=0.5),
        legend.position ='None',
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_text(size=13),
        legend.text=element_text(size=11))

## Potential Set 2 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# One Causal Variant
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_one_causal_results,
                                       other_df=top_pics2_one_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_one_causal_results,
                                    other_df=stable_pics2_one_causal_results)

post_prob_df_1c <- data.frame(Stable=stable_pics2_match_df[
  !stable_pics2_match_df$focal_SNP2_match_other_SNP2,]$PostProb2,
  Top=top_pics2_match_df[
    !stable_pics2_match_df$focal_SNP2_match_other_SNP2,]$PostProb2)
post_prob_df_1c$NoCausalVariants <- rep(1,nrow(post_prob_df_1c))

# Two Causal Variants
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_two_causal_results,
                                       other_df=top_pics2_two_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_two_causal_results,
                                    other_df=stable_pics2_two_causal_results)

post_prob_df_2c <- data.frame(Stable=stable_pics2_match_df[
  !stable_pics2_match_df$focal_SNP2_match_other_SNP2,]$PostProb2,
  Top=top_pics2_match_df[
    !stable_pics2_match_df$focal_SNP2_match_other_SNP2,]$PostProb2)
post_prob_df_2c$NoCausalVariants <- rep(2,nrow(post_prob_df_2c))

# Three Causal Variants
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_three_causal_results,
                                       other_df=top_pics2_three_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_three_causal_results,
                                    other_df=stable_pics2_three_causal_results)

post_prob_df_3c <- data.frame(Stable=stable_pics2_match_df[
  !stable_pics2_match_df$focal_SNP2_match_other_SNP2,]$PostProb2,
  Top=top_pics2_match_df[
    !stable_pics2_match_df$focal_SNP2_match_other_SNP2,]$PostProb2)
post_prob_df_3c$NoCausalVariants <- rep(3,nrow(post_prob_df_3c))

post_prob_df <- rbind(post_prob_df_1c,
                      post_prob_df_2c,
                      post_prob_df_3c)

supp_fig8B <- ggplot(post_prob_df,aes(x=Top,y=Stable)) +
  geom_point(aes(colour=factor(NoCausalVariants)),
             alpha=0.4) +
  geom_abline(slope=1, intercept=0,
              lty="dashed") +
  theme_bw() +
  ylim(c(0,1)) +
  xlim(c(0,1)) +
  xlab("Top") +
  ylab("Stable") +
  ggtitle("Potential Set 2") +
  labs(colour="Number of Causal Variants") +
  theme(plot.title=element_text(hjust=0.5),
        legend.position ='None',
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_text(size=13),
        legend.text=element_text(size=11))

## Potential Set 3 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# One Causal Variant
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_one_causal_results,
                                       other_df=top_pics2_one_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_one_causal_results,
                                    other_df=stable_pics2_one_causal_results)

post_prob_df_1c <- data.frame(Stable=stable_pics2_match_df[
  !stable_pics2_match_df$focal_SNP3_match_other_SNP3,]$PostProb3,
  Top=top_pics2_match_df[
    !stable_pics2_match_df$focal_SNP3_match_other_SNP3,]$PostProb3)
post_prob_df_1c$NoCausalVariants <- rep(1,nrow(post_prob_df_1c))

# Two Causal Variants
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_two_causal_results,
                                       other_df=top_pics2_two_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_two_causal_results,
                                    other_df=stable_pics2_two_causal_results)

post_prob_df_2c <- data.frame(Stable=stable_pics2_match_df[
  !stable_pics2_match_df$focal_SNP3_match_other_SNP3,]$PostProb3,
  Top=top_pics2_match_df[
    !stable_pics2_match_df$focal_SNP3_match_other_SNP3,]$PostProb3)
post_prob_df_2c$NoCausalVariants <- rep(2,nrow(post_prob_df_2c))

# Three Causal Variants
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_three_causal_results,
                                       other_df=top_pics2_three_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_three_causal_results,
                                    other_df=stable_pics2_three_causal_results)

post_prob_df_3c <- data.frame(Stable=stable_pics2_match_df[
  !stable_pics2_match_df$focal_SNP3_match_other_SNP3,]$PostProb3,
  Top=top_pics2_match_df[
    !stable_pics2_match_df$focal_SNP3_match_other_SNP3,]$PostProb3)
post_prob_df_3c$NoCausalVariants <- rep(3,nrow(post_prob_df_3c))

post_prob_df <- rbind(post_prob_df_1c,
                      post_prob_df_2c,
                      post_prob_df_3c)

supp_fig8C <- ggplot(post_prob_df,aes(x=Top,y=Stable)) +
  geom_point(aes(colour=factor(NoCausalVariants)),
             alpha=0.4) +
  geom_abline(slope=1, intercept=0,
              lty="dashed") +
  theme_bw() +
  ylim(c(0,1)) +
  xlim(c(0,1)) +
  xlab("Top") +
  ylab("Stable") +
  ggtitle("Potential Set 3") +
  labs(colour="No. Causal\n Variants") +
  theme(plot.title=element_text(hjust=0.5),
        legend.position ='bottom',
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_text(size=13),
        legend.text=element_text(size=11))

## Combine plots <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
library(grid)
library(ggpubr)

supp_fig8_plot <- ggarrange(supp_fig8A,supp_fig8B,supp_fig8C, 
          nrow=1,common.legend = TRUE,legend="bottom")

supp_fig8_plot<- annotate_figure(supp_fig8_plot, 
                top=text_grob(
                  "Posterior Probabilities for Non-matching Variants (Simulated Gene Expression)", 
                  size = 16))
# supp_fig8_plot <- gridExtra::grid.arrange(supp_fig8A,supp_fig8B,supp_fig8C,
#                                           nrow=1,widths = c(2.4, 2.4, 2.8),
#                                           top = textGrob("Posterior Probabilities for Non-matching Variants",
#                                                          gp=gpar(fontsize=15)))

ggsave(supp_fig8_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig8.pdf"),
       height=4.5,
       width=12)

## Fig 9. Difference in posterior probability between matching variants ========
## Potential Set 1 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_one_causal_results,
                                       other_df=top_pics2_one_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_one_causal_results,
                                    other_df=stable_pics2_one_causal_results)

post_prob_df_1c <- data.frame(Stable=stable_pics2_match_df[
  stable_pics2_match_df$focal_SNP1_match_other_SNP1,]$PostProb1,
  Top=top_pics2_match_df[
    stable_pics2_match_df$focal_SNP1_match_other_SNP1,]$PostProb1)
post_prob_df_1c$NoCausalVariants <- rep(1,nrow(post_prob_df_1c))

# Two Causal Variants
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_two_causal_results,
                                       other_df=top_pics2_two_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_two_causal_results,
                                    other_df=stable_pics2_two_causal_results)

post_prob_df_2c <- data.frame(Stable=stable_pics2_match_df[
  stable_pics2_match_df$focal_SNP1_match_other_SNP1,]$PostProb1,
  Top=top_pics2_match_df[
    stable_pics2_match_df$focal_SNP1_match_other_SNP1,]$PostProb1)
post_prob_df_2c$NoCausalVariants <- rep(2,nrow(post_prob_df_2c))

# Three Causal Variants
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_three_causal_results,
                                       other_df=top_pics2_three_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_three_causal_results,
                                    other_df=stable_pics2_three_causal_results)

post_prob_df_3c <- data.frame(Stable=stable_pics2_match_df[
  stable_pics2_match_df$focal_SNP1_match_other_SNP1,]$PostProb1,
  Top=top_pics2_match_df[
    stable_pics2_match_df$focal_SNP1_match_other_SNP1,]$PostProb1)
post_prob_df_3c$NoCausalVariants <- rep(3,nrow(post_prob_df_3c))

post_prob_df <- rbind(post_prob_df_1c,
                      post_prob_df_2c,
                      post_prob_df_3c)

supp_fig9A <- ggplot(post_prob_df,aes(x=Top,y=Stable)) +
  geom_point(aes(colour=factor(NoCausalVariants)),
             alpha=0.4) +
  geom_abline(slope=1, intercept=0,
              lty="dashed") +
  theme_bw() +
  ylim(c(0,1)) +
  xlim(c(0,1)) +
  xlab("Top") +
  ylab("Stable") +
  ggtitle("Potential Set 1") +
  labs(colour="Number of Causal Variants") +
  theme(plot.title=element_text(hjust=0.5),
        legend.position ='None',
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_text(size=13),
        legend.text=element_text(size=11))

## Potential Set 2 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# One Causal Variant
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_one_causal_results,
                                       other_df=top_pics2_one_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_one_causal_results,
                                    other_df=stable_pics2_one_causal_results)

post_prob_df_1c <- data.frame(Stable=stable_pics2_match_df[
  stable_pics2_match_df$focal_SNP2_match_other_SNP2,]$PostProb2,
  Top=top_pics2_match_df[
    stable_pics2_match_df$focal_SNP2_match_other_SNP2,]$PostProb2)
post_prob_df_1c$NoCausalVariants <- rep(1,nrow(post_prob_df_1c))

# Two Causal Variants
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_two_causal_results,
                                       other_df=top_pics2_two_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_two_causal_results,
                                    other_df=stable_pics2_two_causal_results)

post_prob_df_2c <- data.frame(Stable=stable_pics2_match_df[
  stable_pics2_match_df$focal_SNP2_match_other_SNP2,]$PostProb2,
  Top=top_pics2_match_df[
    stable_pics2_match_df$focal_SNP2_match_other_SNP2,]$PostProb2)
post_prob_df_2c$NoCausalVariants <- rep(2,nrow(post_prob_df_2c))

# Three Causal Variants
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_three_causal_results,
                                       other_df=top_pics2_three_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_three_causal_results,
                                    other_df=stable_pics2_three_causal_results)

post_prob_df_3c <- data.frame(Stable=stable_pics2_match_df[
  stable_pics2_match_df$focal_SNP2_match_other_SNP2,]$PostProb2,
  Top=top_pics2_match_df[
    stable_pics2_match_df$focal_SNP2_match_other_SNP2,]$PostProb2)
post_prob_df_3c$NoCausalVariants <- rep(3,nrow(post_prob_df_3c))

post_prob_df <- rbind(post_prob_df_1c,
                      post_prob_df_2c,
                      post_prob_df_3c)

supp_fig9B <- ggplot(post_prob_df,aes(x=Top,y=Stable)) +
  geom_point(aes(colour=factor(NoCausalVariants)),
             alpha=0.4) +
  geom_abline(slope=1, intercept=0,
              lty="dashed") +
  theme_bw() +
  ylim(c(0,1)) +
  xlim(c(0,1)) +
  xlab("Top") +
  ylab("Stable") +
  ggtitle("Potential Set 2") +
  labs(colour="Number of Causal Variants") +
  theme(plot.title=element_text(hjust=0.5),
        legend.position ='None',
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_text(size=13),
        legend.text=element_text(size=11))

## Potential Set 3 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# One Causal Variant
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_one_causal_results,
                                       other_df=top_pics2_one_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_one_causal_results,
                                    other_df=stable_pics2_one_causal_results)

post_prob_df_1c <- data.frame(Stable=stable_pics2_match_df[
  stable_pics2_match_df$focal_SNP3_match_other_SNP3,]$PostProb3,
  Top=top_pics2_match_df[
    stable_pics2_match_df$focal_SNP3_match_other_SNP3,]$PostProb3)
post_prob_df_1c$NoCausalVariants <- rep(1,nrow(post_prob_df_1c))

# Two Causal Variants
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_two_causal_results,
                                       other_df=top_pics2_two_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_two_causal_results,
                                    other_df=stable_pics2_two_causal_results)

post_prob_df_2c <- data.frame(Stable=stable_pics2_match_df[
  stable_pics2_match_df$focal_SNP3_match_other_SNP3,]$PostProb3,
  Top=top_pics2_match_df[
    stable_pics2_match_df$focal_SNP3_match_other_SNP3,]$PostProb3)
post_prob_df_2c$NoCausalVariants <- rep(2,nrow(post_prob_df_2c))

# Three Causal Variants
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_three_causal_results,
                                       other_df=top_pics2_three_causal_results)
top_pics2_match_df <- getMatchStats(focal_df=top_pics2_three_causal_results,
                                    other_df=stable_pics2_three_causal_results)

post_prob_df_3c <- data.frame(Stable=stable_pics2_match_df[
  stable_pics2_match_df$focal_SNP3_match_other_SNP3,]$PostProb3,
  Top=top_pics2_match_df[
    stable_pics2_match_df$focal_SNP3_match_other_SNP3,]$PostProb3)
post_prob_df_3c$NoCausalVariants <- rep(3,nrow(post_prob_df_3c))

post_prob_df <- rbind(post_prob_df_1c,
                      post_prob_df_2c,
                      post_prob_df_3c)

supp_fig9C <- ggplot(post_prob_df,aes(x=Top,y=Stable)) +
  geom_point(aes(colour=factor(NoCausalVariants)),
             alpha=0.4) +
  geom_abline(slope=1, intercept=0,
              lty="dashed") +
  theme_bw() +
  ylim(c(0,1)) +
  xlim(c(0,1)) +
  xlab("Top") +
  ylab("Stable") +
  ggtitle("Potential Set 3") +
  labs(colour="No. Causal\n Variants") +
  theme(plot.title=element_text(hjust=0.5),
        legend.position ='bottom',
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_text(size=13),
        legend.text=element_text(size=11))

## Combine plots <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
library(grid)
library(ggpubr)

supp_fig9_plot <- ggarrange(supp_fig9A,supp_fig9B,supp_fig9C, 
                            nrow=1,common.legend = TRUE,legend="bottom")

supp_fig9_plot<- annotate_figure(supp_fig9_plot, 
                                 top=text_grob(
                                   "Posterior Probabilities for Matching Variants (Simulated Gene Expression)", 
                                   size = 16))
ggsave(supp_fig9_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig9.pdf"),
       height=4.5,
       width=12)

## Fig 10-12. SuSiE vs PICS ====================================================
## Top SuSiE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Get causal variant recovery counts 
all_top_susie_results$NoCausalSNPsRecovered_1PS <- 
  mapply(countCausalSNPs,
         all_top_susie_results$CausalSNPs,
         all_top_susie_results$SNP1,
         rep(NA,nrow(all_top_susie_results)),
         rep(NA,nrow(all_top_susie_results)))
all_top_susie_results$NoCausalSNPsRecovered_2PS <- 
  mapply(countCausalSNPs,
         all_top_susie_results$CausalSNPs,
         all_top_susie_results$SNP1,
         all_top_susie_results$SNP2,
         rep(NA,nrow(all_top_susie_results)))
all_top_susie_results$NoCausalSNPsRecovered_3PS <- 
  mapply(countCausalSNPs,
         all_top_susie_results$CausalSNPs,
         all_top_susie_results$SNP1,
         all_top_susie_results$SNP2,
         all_top_susie_results$SNP3)

# Split results by number of causal variants simulated
top_susie_one_causal_results <- all_top_susie_results %>% 
  subset(S==1)
top_susie_two_causal_results <- all_top_susie_results %>% 
  subset(S==2)
top_susie_three_causal_results <- all_top_susie_results %>% 
  subset(S==3)

phi_vector <- c(0.05,0.1,0.2,0.4)

top_susie_trend_df <- NULL
# for SNR ratio parameter
for (i in 1:4) {
  # for no. sets included
  for (j in 1:3) {
    # one causal
    one_causal <- table(
      top_susie_one_causal_results$Phi,
      top_susie_one_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(one_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 1",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    top_susie_trend_df <- rbind(top_susie_trend_df, ci_table)
    
    # two causal
    two_causal <- table(
      top_susie_two_causal_results$Phi,
      top_susie_two_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(two_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 2",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    top_susie_trend_df <- rbind(top_susie_trend_df, ci_table)
    
    # three causal
    three_causal <- table(
      top_susie_three_causal_results$Phi,
      top_susie_three_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(three_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 3",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    top_susie_trend_df <- rbind(top_susie_trend_df, ci_table)
  }
}

top_susie_trend_df$Algorithm <- rep("Top SuSiE", nrow(top_susie_trend_df))

## Combine Top SuSiE and Stable PICS to compare <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
top_susie_vs_stable_pics2_trend_df <- rbind(top_susie_trend_df,
                                            stable_pics2_trend_df)
top_susie_vs_stable_pics2_trend_df$SNR <- paste0(
  "SNR = ", 
  round(top_susie_vs_stable_pics2_trend_df$Phi/(1-top_susie_vs_stable_pics2_trend_df$Phi),3))
top_susie_vs_stable_pics2_trend_1c_df <- top_susie_vs_stable_pics2_trend_df %>% 
  subset(NoCausalVarSim == "No. Causal Var = 1")
top_susie_vs_stable_pics2_trend_2c_df <- top_susie_vs_stable_pics2_trend_df %>% 
  subset(NoCausalVarSim == "No. Causal Var = 2")
top_susie_vs_stable_pics2_trend_3c_df <- top_susie_vs_stable_pics2_trend_df %>% 
  subset(NoCausalVarSim == "No. Causal Var = 3")

## Generate plots <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# One causal variant
supp_fig10_plot <- ggplot(top_susie_vs_stable_pics2_trend_1c_df, 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoSetsInclude~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  ggtitle("Simulations with 1 Causal Variant\n(Stable PICS vs Top SuSiE)\n\n") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig10_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig10.pdf"),
       height=9,
       width=9)

# Two causal variants
supp_fig11_plot <- ggplot(top_susie_vs_stable_pics2_trend_2c_df, 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoSetsInclude~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  ggtitle("Simulations with 2 Causal Variants\n(Stable PICS vs Top SuSiE)\n\n") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig11_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig11.pdf"),
       height=9,
       width=9)

# Three causal variants
supp_fig12_plot <- ggplot(top_susie_vs_stable_pics2_trend_3c_df, 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoSetsInclude~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  ggtitle("Simulations with 3 Causal Variants\n(Stable PICS vs Top SuSiE)\n\n") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig12_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig12.pdf"),
       height=9,
       width=9)

## Fig 13-15 Combined vs Top vs Stable PICS ====================================
# Split results by no. causal variants and SNR-stratified analysis.
# This is similar to Figs 5-7 and 10, where we split results by no. causal variants
# and then show how the performance improves as a function of the number of sets
# included. We plot this for PICS (combined vs stable vs top), then plot for 
# SuSiE (stable vs top)
## Combined PICS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Get causal variant recovery counts 
all_combined_pics2_results$NoCausalSNPsRecovered_1PS <- 
  mapply(countCausalSNPs,
         all_combined_pics2_results$CausalSNPs,
         all_combined_pics2_results$SNP1,
         rep(NA,nrow(all_combined_pics2_results)),
         rep(NA,nrow(all_combined_pics2_results)))
all_combined_pics2_results$NoCausalSNPsRecovered_2PS <- 
  mapply(countCausalSNPs,
         all_combined_pics2_results$CausalSNPs,
         all_combined_pics2_results$SNP1,
         all_combined_pics2_results$SNP2,
         rep(NA,nrow(all_combined_pics2_results)))
all_combined_pics2_results$NoCausalSNPsRecovered_3PS <- 
  mapply(countCausalSNPs,
         all_combined_pics2_results$CausalSNPs,
         all_combined_pics2_results$SNP1,
         all_combined_pics2_results$SNP2,
         all_combined_pics2_results$SNP3)

# Split results by number of causal variants simulated
combined_pics2_one_causal_results <- all_combined_pics2_results %>% 
  subset(S==1)
combined_pics2_two_causal_results <- all_combined_pics2_results %>% 
  subset(S==2)
combined_pics2_three_causal_results <- all_combined_pics2_results %>% 
  subset(S==3)

phi_vector <- c(0.05,0.1,0.2,0.4)

combined_pics2_trend_df <- NULL 
# for SNR ratio parameter
for (i in 1:4) {
  # for no. sets included
  for (j in 1:3) {
    # one causal
    one_causal <- table(
      combined_pics2_one_causal_results$Phi,
      combined_pics2_one_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(one_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 1",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    combined_pics2_trend_df <- rbind(combined_pics2_trend_df, ci_table)
    
    # two causal
    two_causal <- table(
      combined_pics2_two_causal_results$Phi,
      combined_pics2_two_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(two_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 2",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    combined_pics2_trend_df <- rbind(combined_pics2_trend_df, ci_table)
    
    # three causal
    three_causal <- table(
      combined_pics2_three_causal_results$Phi,
      combined_pics2_three_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(three_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 3",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    combined_pics2_trend_df <- rbind(combined_pics2_trend_df, ci_table)
  }
}

combined_pics2_trend_df$Algorithm <- rep("Combined PICS", nrow(combined_pics2_trend_df))

## Top PICS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Get causal variant recovery counts 
all_top_pics2_results$NoCausalSNPsRecovered_1PS <- 
  mapply(countCausalSNPs,
         all_top_pics2_results$CausalSNPs,
         all_top_pics2_results$SNP1,
         rep(NA,nrow(all_top_pics2_results)),
         rep(NA,nrow(all_top_pics2_results)))
all_top_pics2_results$NoCausalSNPsRecovered_2PS <- 
  mapply(countCausalSNPs,
         all_top_pics2_results$CausalSNPs,
         all_top_pics2_results$SNP1,
         all_top_pics2_results$SNP2,
         rep(NA,nrow(all_top_pics2_results)))
all_top_pics2_results$NoCausalSNPsRecovered_3PS <- 
  mapply(countCausalSNPs,
         all_top_pics2_results$CausalSNPs,
         all_top_pics2_results$SNP1,
         all_top_pics2_results$SNP2,
         all_top_pics2_results$SNP3)

# Split results by number of causal variants simulated
top_pics2_one_causal_results <- all_top_pics2_results %>% 
  subset(S==1)
top_pics2_two_causal_results <- all_top_pics2_results %>% 
  subset(S==2)
top_pics2_three_causal_results <- all_top_pics2_results %>% 
  subset(S==3)

phi_vector <- c(0.05,0.1,0.2,0.4)

top_pics2_trend_df <- NULL 
# for SNR ratio parameter
for (i in 1:4) {
  # for no. sets included
  for (j in 1:3) {
    # one causal
    one_causal <- table(
      top_pics2_one_causal_results$Phi,
      top_pics2_one_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(one_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 1",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    top_pics2_trend_df <- rbind(top_pics2_trend_df, ci_table)
    
    # two causal
    two_causal <- table(
      top_pics2_two_causal_results$Phi,
      top_pics2_two_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(two_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 2",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    top_pics2_trend_df <- rbind(top_pics2_trend_df, ci_table)
    
    # three causal
    three_causal <- table(
      top_pics2_three_causal_results$Phi,
      top_pics2_three_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(three_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 3",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    top_pics2_trend_df <- rbind(top_pics2_trend_df, ci_table)
  }
}

top_pics2_trend_df$Algorithm <- rep("Top PICS", nrow(top_pics2_trend_df))

## Combine all PICS results to compare <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
all_pics2_trend_df <- rbind(combined_pics2_trend_df,
                            stable_pics2_trend_df,
                            top_pics2_trend_df)
all_pics2_trend_df$SNR <- paste0(
  "SNR = ", 
  round(all_pics2_trend_df$Phi/(1-all_pics2_trend_df$Phi),3))
all_pics2_trend_1c_df <- all_pics2_trend_df %>% 
  subset(NoCausalVarSim == "No. Causal Var = 1")
all_pics2_trend_2c_df <- all_pics2_trend_df %>% 
  subset(NoCausalVarSim == "No. Causal Var = 2")
all_pics2_trend_3c_df <- all_pics2_trend_df %>% 
  subset(NoCausalVarSim == "No. Causal Var = 3")

## Generate plots <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# One causal variant
supp_fig13_plot <- ggplot(all_pics2_trend_1c_df, 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoSetsInclude~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#f38742","#4d8dc3","#4ea26a")) + 
  scale_colour_manual(values=c("#f16913","#2171b5","#238b45")) +
  ggtitle("Simulations with 1 Causal Variant (PICS Algorithms)\n\n") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig13_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig13.pdf"),
       height=9,
       width=9)

# Two causal variants
supp_fig14_plot <- ggplot(all_pics2_trend_2c_df, 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoSetsInclude~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#f38742","#4d8dc3","#4ea26a")) + 
  scale_colour_manual(values=c("#f16913","#2171b5","#238b45")) +
  ggtitle("Simulations with 2 Causal Variants (PICS Algorithms)\n\n") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig14_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig14.pdf"),
       height=9,
       width=9)

# Three causal variants
supp_fig15_plot <- ggplot(all_pics2_trend_3c_df, 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoSetsInclude~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#f38742","#4d8dc3","#4ea26a")) + 
  scale_colour_manual(values=c("#f16913","#2171b5","#238b45")) +
  ggtitle("Simulations with 3 Causal Variants (PICS Algorithms)\n\n") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig15_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig15.pdf"),
       height=9,
       width=9)

## Fig 16-18: Top vs Stable SuSiE ==============================================
## Combine Top SuSiE and Stable PICS to compare <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
top_vs_stable_susie_trend_df <- rbind(top_susie_trend_df,
                                      stable_susie_trend_df)
top_vs_stable_susie_trend_df$SNR <- paste0(
  "SNR = ", 
  round(top_vs_stable_susie_trend_df$Phi/(1-top_vs_stable_susie_trend_df$Phi),3))
top_vs_stable_susie_trend_1c_df <- top_vs_stable_susie_trend_df %>% 
  subset(NoCausalVarSim == "No. Causal Var = 1")
top_vs_stable_susie_trend_2c_df <- top_vs_stable_susie_trend_df %>% 
  subset(NoCausalVarSim == "No. Causal Var = 2")
top_vs_stable_susie_trend_3c_df <- top_vs_stable_susie_trend_df %>% 
  subset(NoCausalVarSim == "No. Causal Var = 3")

## Generate plots <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# One causal variant
supp_fig16_plot <- ggplot(top_vs_stable_susie_trend_1c_df, 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoSetsInclude~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#4d8dc3","#4ea26a")) + 
  scale_colour_manual(values=c("#2171b5","#238b45")) +
  ggtitle("Simulations with 1 Causal Variant (SuSiE Algorithms)\n\n") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig16_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig16.pdf"),
       height=9,
       width=9)

# Two causal variants
supp_fig17_plot <- ggplot(top_vs_stable_susie_trend_2c_df, 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoSetsInclude~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#4d8dc3","#4ea26a")) + 
  scale_colour_manual(values=c("#2171b5","#238b45")) +
  ggtitle("Simulations with 2 Causal Variants (SuSiE Algorithms)\n\n") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig17_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig17.pdf"),
       height=9,
       width=9)

# Three causal variants
supp_fig18_plot <- ggplot(top_vs_stable_susie_trend_3c_df, 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoSetsInclude~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#4d8dc3","#4ea26a")) + 
  scale_colour_manual(values=c("#2171b5","#238b45")) +
  ggtitle("Simulations with 3 Causal Variants (SuSiE Algorithms)\n\n") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig18_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig18.pdf"),
       height=9,
       width=9)

## Fig 19-23: Matching vs Non-matching =========================================
# Have matching, non-matching top and non-matching stable for both SuSiE and
# PICS

## One Causal Variant --- Recovery Frequency curves <<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Construct matching vs non-matching df
pics_1c_match_vs_nonmatch <- getMatchRecProbDF2(top_pics2_one_causal_results,
                                                stable_pics2_one_causal_results)
pics_1c_match_vs_nonmatch$Algorithm <- rep("PICS", 
                                           nrow(pics_1c_match_vs_nonmatch))

susie_1c_match_vs_nonmatch <- getMatchRecProbDF2(top_susie_one_causal_results,
                                                 stable_susie_one_causal_results)
susie_1c_match_vs_nonmatch$Algorithm <- rep("SuSiE", 
                                            nrow(susie_1c_match_vs_nonmatch))

# To solve jitter problem (improve visualization), manually perturb SNR by +0.004
pics_1c_match_vs_nonmatch$SNR <- pics_1c_match_vs_nonmatch$SNR+0.004

# To solve jitter problem (improve visualization), manually perturb SNR by -0.004
susie_1c_match_vs_nonmatch$SNR <- susie_1c_match_vs_nonmatch$SNR-0.004

## Combine PICS and SuSiE
combined_1c_match_vs_nonmatch <- rbind(pics_1c_match_vs_nonmatch,
                                       susie_1c_match_vs_nonmatch)

combined_1c_match_vs_nonmatch$PS <- paste0("Credible / ", combined_1c_match_vs_nonmatch$PS)
## Construct plot
supp_fig19_plot <- ggplot(combined_1c_match_vs_nonmatch,
                          aes(x=SNR,y=RecProb)) +
  geom_point(aes(colour=Type,
                 shape=Algorithm),
             alpha=0.8,
             size=2)+
  geom_errorbar(aes(ymin=RecProb-1.96*BinomSE,
                    ymax=RecProb+1.96*BinomSE,
                    colour=Type),
                alpha=0.8,
                width=.05)+
  geom_line(aes(colour=Type,
                lty=Algorithm,
                group=interaction(Algorithm,Type)),
            alpha=0.8) +
  facet_wrap(.~PS,dir="v",scales="free") +
  theme_bw() +
  xlim(c(0.05,0.67)) +
  #ylim(c(0,0.75))+
  scale_x_continuous(breaks=c(0.05,0.1,0.25,0.67)) +
  #scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.75)) +
  ylab("Recovery Frequency") +
  xlab("Signal-to-Noise Ratio (SNR)") +
  ggtitle("Simulations with 1 Causal Variant") +
  scale_color_manual(values=c("#d7301f","#2171b5","#238b45")) +
  scale_shape_manual(values=c("PICS" = 1, "SuSiE" = 16)) +
  #guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=15),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig19_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig19.pdf"),
       height=14,
       width=6)

## Two Causal Variants --- Histograms <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Construct PICS matching vs non-matching df
pics_2c_match_vs_nonmatch <- getCausalVarCountsDF(top_pics2_two_causal_results,
                                                  stable_pics2_two_causal_results)
pics_2c_match_vs_nonmatch$SNR <- 
  paste0("SNR = ", round(pics_2c_match_vs_nonmatch$Phi/
                           (1-pics_2c_match_vs_nonmatch$Phi),3))

# Construct SuSiE matching vs non-matching df 
susie_2c_match_vs_nonmatch <- getCausalVarCountsDF(top_susie_two_causal_results,
                                                   stable_susie_two_causal_results,
                                                   algo="SuSiE")
susie_2c_match_vs_nonmatch$SNR <- 
  paste0("SNR = ", round(susie_2c_match_vs_nonmatch$Phi/
                           (1-susie_2c_match_vs_nonmatch$Phi),3))

## Construct plots 
# PICS
supp_fig20_plot <- ggplot(pics_2c_match_vs_nonmatch, 
       aes(x=factor(NoCausalRec,levels=c(0,1,2)),y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Type)) +
  geom_errorbar(aes(x=factor(NoCausalRec,levels=c(0,1,2)),
                    ymin=lwr.ci,ymax=upr.ci,colour=Type),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(Potential_Set~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#df594b","#4d8dc3","#4ea26a")) + 
  scale_colour_manual(values=c("#d7301f","#2171b5","#238b45")) +
  ggtitle("Simulations with 2 Causal Variants (PICS Algorithms)") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig20_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig20.pdf"),
       height=9,
       width=9)

# SuSiE
supp_fig21_plot <- ggplot(susie_2c_match_vs_nonmatch, 
                          aes(x=factor(NoCausalRec,levels=c(0,1,2)),y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Type)) +
  geom_errorbar(aes(x=factor(NoCausalRec,levels=c(0,1,2)),
                    ymin=lwr.ci,ymax=upr.ci,colour=Type),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(Potential_Set~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#df594b","#4d8dc3","#4ea26a")) + 
  scale_colour_manual(values=c("#d7301f","#2171b5","#238b45")) +
  ggtitle("Simulations with 2 Causal Variants (SuSiE Algorithms)") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig21_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig21.pdf"),
       height=9,
       width=9)

## Three Causal Variants --- Histograms <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Construct PICS matching vs non-matching df
pics_3c_match_vs_nonmatch <- getCausalVarCountsDF3(top_pics2_three_causal_results,
                                                   stable_pics2_three_causal_results)
pics_3c_match_vs_nonmatch$SNR <- 
  paste0("SNR = ", round(pics_3c_match_vs_nonmatch$Phi/
                           (1-pics_3c_match_vs_nonmatch$Phi),3))

# Construct SuSiE matching vs non-matching df 
susie_3c_match_vs_nonmatch <- getCausalVarCountsDF3(top_susie_three_causal_results,
                                                   stable_susie_three_causal_results,
                                                   algo="SuSiE")
susie_3c_match_vs_nonmatch$SNR <- 
  paste0("SNR = ", round(susie_3c_match_vs_nonmatch$Phi/
                           (1-susie_3c_match_vs_nonmatch$Phi),3))

## Construct plots 
# PICS
supp_fig22_plot <- ggplot(pics_3c_match_vs_nonmatch, 
                          aes(x=factor(NoCausalRec,levels=c(0,1,2,3)),y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Type)) +
  geom_errorbar(aes(x=factor(NoCausalRec,levels=c(0,1,2,3)),
                    ymin=lwr.ci,ymax=upr.ci,colour=Type),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(Potential_Set~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#df594b","#4d8dc3","#4ea26a")) + 
  scale_colour_manual(values=c("#d7301f","#2171b5","#238b45")) +
  ggtitle("Simulations with 3 Causal Variants (PICS Algorithms)") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig22_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig22.pdf"),
       height=9,
       width=9)

# SuSiE
supp_fig23_plot <- ggplot(susie_3c_match_vs_nonmatch, 
                          aes(x=factor(NoCausalRec,levels=c(0,1,2,3)),y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Type)) +
  geom_errorbar(aes(x=factor(NoCausalRec,levels=c(0,1,2,3)),
                    ymin=lwr.ci,ymax=upr.ci,colour=Type),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(Potential_Set~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#df594b","#4d8dc3","#4ea26a")) + 
  scale_colour_manual(values=c("#d7301f","#2171b5","#238b45")) +
  ggtitle("Simulations with 3 Causal Variants (SuSiE Algorithms)") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig23_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig23.pdf"),
       height=9,
       width=9)

## Fig 24-26: SNR-stratified results of Plain vs Stable PICS ===================
# Plain PICS (Env Het) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
big_plain_pics_df$NoCausalSNPsRecovered_1PS <- 
  mapply(countCausalSNPs,
         big_plain_pics_df$CausalSNPs,
         big_plain_pics_df$SNP1,
         rep(NA,nrow(big_plain_pics_df)),
         rep(NA,nrow(big_plain_pics_df)))
big_plain_pics_df$NoCausalSNPsRecovered_2PS <- 
  mapply(countCausalSNPs,
         big_plain_pics_df$CausalSNPs,
         big_plain_pics_df$SNP1,
         big_plain_pics_df$SNP2,
         rep(NA,nrow(big_plain_pics_df)))
big_plain_pics_df$NoCausalSNPsRecovered_3PS <- 
  mapply(countCausalSNPs,
         big_plain_pics_df$CausalSNPs,
         big_plain_pics_df$SNP1,
         big_plain_pics_df$SNP2,
         big_plain_pics_df$SNP3)

# Split results by number of causal variants simulated
plain_pics2_env_het_one_causal_results <- big_plain_pics_df %>% 
  subset(S==1)
plain_pics2_env_het_two_causal_results <- big_plain_pics_df %>% 
  subset(S==2)
plain_pics2_env_het_three_causal_results <- big_plain_pics_df %>% 
  subset(S==3)

phi_vector <- c(0.05,0.1,0.2,0.4)

plain_pics2_env_het_trend_df <- NULL 
# for SNR ratio parameter
for (i in 1:4) {
  # for no. sets included
  for (j in 1:3) {
    # one causal
    one_causal <- table(
      plain_pics2_env_het_one_causal_results$Phi,
      plain_pics2_env_het_one_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(one_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 1",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    plain_pics2_env_het_trend_df <- rbind(plain_pics2_env_het_trend_df, ci_table)
    
    # two causal
    two_causal <- table(
      plain_pics2_env_het_two_causal_results$Phi,
      plain_pics2_env_het_two_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(two_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 2",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    plain_pics2_env_het_trend_df <- rbind(plain_pics2_env_het_trend_df, ci_table)
    
    # three causal
    three_causal <- table(
      plain_pics2_env_het_three_causal_results$Phi,
      plain_pics2_env_het_three_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(three_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 3",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    plain_pics2_env_het_trend_df <- rbind(plain_pics2_env_het_trend_df, ci_table)
  }
}

plain_pics2_env_het_trend_df$Algorithm <- rep("Plain PICS", 
                                              nrow(plain_pics2_env_het_trend_df))

# Stable PICS (Env Het) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
big_stable_pics_df$NoCausalSNPsRecovered_1PS <- 
  mapply(countCausalSNPs,
         big_stable_pics_df$CausalSNPs,
         big_stable_pics_df$SNP1,
         rep(NA,nrow(big_stable_pics_df)),
         rep(NA,nrow(big_stable_pics_df)))
big_stable_pics_df$NoCausalSNPsRecovered_2PS <- 
  mapply(countCausalSNPs,
         big_stable_pics_df$CausalSNPs,
         big_stable_pics_df$SNP1,
         big_stable_pics_df$SNP2,
         rep(NA,nrow(big_stable_pics_df)))
big_stable_pics_df$NoCausalSNPsRecovered_3PS <- 
  mapply(countCausalSNPs,
         big_stable_pics_df$CausalSNPs,
         big_stable_pics_df$SNP1,
         big_stable_pics_df$SNP2,
         big_stable_pics_df$SNP3)

# Split results by number of causal variants simulated
stable_pics2_env_het_one_causal_results <- big_stable_pics_df %>% 
  subset(S==1)
stable_pics2_env_het_two_causal_results <- big_stable_pics_df %>% 
  subset(S==2)
stable_pics2_env_het_three_causal_results <- big_stable_pics_df %>% 
  subset(S==3)

phi_vector <- c(0.05,0.1,0.2,0.4)

stable_pics2_env_het_trend_df <- NULL 
# for SNR ratio parameter
for (i in 1:4) {
  # for no. sets included
  for (j in 1:3) {
    # one causal
    one_causal <- table(
      stable_pics2_env_het_one_causal_results$Phi,
      stable_pics2_env_het_one_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(one_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 1",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    stable_pics2_env_het_trend_df <- rbind(stable_pics2_env_het_trend_df, ci_table)
    
    # two causal
    two_causal <- table(
      stable_pics2_env_het_two_causal_results$Phi,
      stable_pics2_env_het_two_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(two_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 2",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    stable_pics2_env_het_trend_df <- rbind(stable_pics2_env_het_trend_df, ci_table)
    
    # three causal
    three_causal <- table(
      stable_pics2_env_het_three_causal_results$Phi,
      stable_pics2_env_het_three_causal_results[[paste0("NoCausalSNPsRecovered_",j,"PS")]])[i,]
    ci_table <- DescTools::MultinomCI(three_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("No. Causal Var = 3",
                                   nrow(ci_table))
    ci_table$NoSetsInclude <- rep(j, nrow(ci_table))
    stable_pics2_env_het_trend_df <- rbind(stable_pics2_env_het_trend_df, ci_table)
  }
}

stable_pics2_env_het_trend_df$Algorithm <- rep("Stable PICS", 
                                              nrow(stable_pics2_env_het_trend_df))

## Combine Env Het trends to compare <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
all_pics2_env_het_trend_df <- rbind(plain_pics2_env_het_trend_df,
                                    stable_pics2_env_het_trend_df)
all_pics2_env_het_trend_df$SNR <- paste0(
  "SNR = ", 
  round(all_pics2_env_het_trend_df$Phi/(1-all_pics2_env_het_trend_df$Phi),3))
all_pics2_env_het_trend_1c_df <- all_pics2_env_het_trend_df %>% 
  subset(NoCausalVarSim == "No. Causal Var = 1")
all_pics2_env_het_trend_2c_df <- all_pics2_env_het_trend_df %>% 
  subset(NoCausalVarSim == "No. Causal Var = 2")
all_pics2_env_het_trend_3c_df <- all_pics2_env_het_trend_df %>% 
  subset(NoCausalVarSim == "No. Causal Var = 3")

## Generate plots <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# One causal variant
supp_fig24_plot <- ggplot(all_pics2_env_het_trend_1c_df, 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoSetsInclude~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#f89acb","#4d8dc3")) + 
  scale_colour_manual(values=c("#f781bf","#2171b5")) +
  ggtitle("Simulations with 1 Causal Variant (Plain vs Stable PICS)\n\n") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig24_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig24.pdf"),
       height=9,
       width=9)

# Two causal variants
supp_fig25_plot <- ggplot(all_pics2_env_het_trend_2c_df, 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoSetsInclude~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#f89acb","#4d8dc3")) + 
  scale_colour_manual(values=c("#f781bf","#2171b5")) +
  ggtitle("Simulations with 2 Causal Variants (Plain vs Stable PICS)\n\n") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig25_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig25.pdf"),
       height=9,
       width=9)

# Three causal variants
supp_fig26_plot <- ggplot(all_pics2_env_het_trend_3c_df, 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoSetsInclude~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#f89acb","#4d8dc3")) + 
  scale_colour_manual(values=c("#f781bf","#2171b5")) +
  ggtitle("Simulations with 3 Causal Variants (Plain vs Stable PICS)\n\n") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig26_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig26.pdf"),
       height=9,
       width=9)

## Fig 27-32: Scenario-stratified results of Plain vs Stable PICS ==============
## Gather histograms for each scenario <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
scenario_names <- names(plain_pics2_env_het_one_causal_results$MODEL %>% table()) 
#[1] "|i-3|" "i=3"   "t=128" "t=16"  "t=256" "t=8"  

## Plain PICS
plain_pics2_scen_env_het_trend_df <- NULL
# for SNR ratio parameter
for (i in 1:4) {
  # for each env het scenario
  for (scen in scenario_names) {
    # subset to just the relevant scenario
    scen_plain_pics2_one_causal <- plain_pics2_env_het_one_causal_results %>% 
      subset(MODEL==scen)
    scen_plain_pics2_two_causal <- plain_pics2_env_het_two_causal_results %>% 
      subset(MODEL==scen)
    scen_plain_pics2_three_causal <- plain_pics2_env_het_three_causal_results %>% 
      subset(MODEL==scen)
    
    # one causal
    one_causal <- table(
      scen_plain_pics2_one_causal$Phi,
      scen_plain_pics2_one_causal[[paste0("NoCausalSNPsRecovered_3PS")]])[i,]
    ci_table <- DescTools::MultinomCI(one_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("1 Causal Variant",
                                   nrow(ci_table))
    ci_table$Scenario <- rep(scen, nrow(ci_table))
    plain_pics2_scen_env_het_trend_df <- rbind(plain_pics2_scen_env_het_trend_df, 
                                               ci_table)
    
    # two causal
    two_causal <- table(
      scen_plain_pics2_two_causal$Phi,
      scen_plain_pics2_two_causal[[paste0("NoCausalSNPsRecovered_3PS")]])[i,]
    ci_table <- DescTools::MultinomCI(two_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("2 Causal Variants",
                                   nrow(ci_table))
    ci_table$Scenario <- rep(scen, nrow(ci_table))
    plain_pics2_scen_env_het_trend_df <- rbind(plain_pics2_scen_env_het_trend_df, 
                                               ci_table)
    
    # three causal
    three_causal <- table(
      scen_plain_pics2_three_causal$Phi,
      scen_plain_pics2_three_causal[[paste0("NoCausalSNPsRecovered_3PS")]])[i,]
    ci_table <- DescTools::MultinomCI(three_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("3 Causal Variants",
                                   nrow(ci_table))
    ci_table$Scenario <- rep(scen, nrow(ci_table))
    plain_pics2_scen_env_het_trend_df <- rbind(plain_pics2_scen_env_het_trend_df, 
                                               ci_table)
  }
}

plain_pics2_scen_env_het_trend_df$Algorithm <- 
  rep("Plain PICS", 
      nrow(plain_pics2_scen_env_het_trend_df))

## Stable PICS
stable_pics2_scen_env_het_trend_df <- NULL
# for SNR ratio parameter
for (i in 1:4) {
  # for each env het scenario
  for (scen in scenario_names) {
    # subset to just the relevant scenario
    scen_stable_pics2_one_causal <- stable_pics2_env_het_one_causal_results %>% 
      subset(MODEL==scen)
    scen_stable_pics2_two_causal <- stable_pics2_env_het_two_causal_results %>% 
      subset(MODEL==scen)
    scen_stable_pics2_three_causal <- stable_pics2_env_het_three_causal_results %>% 
      subset(MODEL==scen)
    
    # one causal
    one_causal <- table(
      scen_stable_pics2_one_causal$Phi,
      scen_stable_pics2_one_causal[[paste0("NoCausalSNPsRecovered_3PS")]])[i,]
    ci_table <- DescTools::MultinomCI(one_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("1 Causal Variant",
                                   nrow(ci_table))
    ci_table$Scenario <- rep(scen, nrow(ci_table))
    stable_pics2_scen_env_het_trend_df <- rbind(stable_pics2_scen_env_het_trend_df, 
                                               ci_table)
    
    # two causal
    two_causal <- table(
      scen_stable_pics2_two_causal$Phi,
      scen_stable_pics2_two_causal[[paste0("NoCausalSNPsRecovered_3PS")]])[i,]
    ci_table <- DescTools::MultinomCI(two_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("2 Causal Variants",
                                   nrow(ci_table))
    ci_table$Scenario <- rep(scen, nrow(ci_table))
    stable_pics2_scen_env_het_trend_df <- rbind(stable_pics2_scen_env_het_trend_df, 
                                               ci_table)
    
    # three causal
    three_causal <- table(
      scen_stable_pics2_three_causal$Phi,
      scen_stable_pics2_three_causal[[paste0("NoCausalSNPsRecovered_3PS")]])[i,]
    ci_table <- DescTools::MultinomCI(three_causal)
    ci_table <- as.data.frame(ci_table)
    ci_table$Levels <- factor(as.numeric(rownames(ci_table)))
    ci_table$Phi <- rep(phi_vector[i], nrow(ci_table))
    ci_table$NoCausalVarSim <- rep("3 Causal Variants",
                                   nrow(ci_table))
    ci_table$Scenario <- rep(scen, nrow(ci_table))
    stable_pics2_scen_env_het_trend_df <- rbind(stable_pics2_scen_env_het_trend_df, 
                                               ci_table)
  }
}

stable_pics2_scen_env_het_trend_df$Algorithm <- 
  rep("Stable PICS", 
      nrow(stable_pics2_scen_env_het_trend_df))

## Combine Env Het trends to compare <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
all_pics2_scen_env_het_trend_df <- rbind(plain_pics2_scen_env_het_trend_df,
                                         stable_pics2_scen_env_het_trend_df)
all_pics2_scen_env_het_trend_df$SNR <- paste0(
  "SNR = ", 
  round(all_pics2_scen_env_het_trend_df$Phi/(1-all_pics2_scen_env_het_trend_df$Phi),3))

## t=8 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
supp_fig27_plot <- ggplot(all_pics2_scen_env_het_trend_df %>% subset(Scenario=="t=8"), 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoCausalVarSim~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#f89acb","#4d8dc3")) + 
  scale_colour_manual(values=c("#f781bf","#2171b5")) +
  ggtitle("t=8 Simulations (Plain vs Stable PICS)") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig27_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig27.pdf"),
       height=9,
       width=9)

## t=16 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
supp_fig28_plot <- ggplot(all_pics2_scen_env_het_trend_df %>% subset(Scenario=="t=16"), 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoCausalVarSim~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#f89acb","#4d8dc3")) + 
  scale_colour_manual(values=c("#f781bf","#2171b5")) +
  ggtitle("t=16 Simulations (Plain vs Stable PICS)") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig28_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig28.pdf"),
       height=9,
       width=9)

## t=128 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
supp_fig29_plot <- ggplot(all_pics2_scen_env_het_trend_df %>% subset(Scenario=="t=128"), 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoCausalVarSim~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#f89acb","#4d8dc3")) + 
  scale_colour_manual(values=c("#f781bf","#2171b5")) +
  ggtitle("t=128 Simulations (Plain vs Stable PICS)") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig29_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig29.pdf"),
       height=9,
       width=9)

## t=256 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
supp_fig30_plot <- ggplot(all_pics2_scen_env_het_trend_df %>% subset(Scenario=="t=256"), 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoCausalVarSim~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#f89acb","#4d8dc3")) + 
  scale_colour_manual(values=c("#f781bf","#2171b5")) +
  ggtitle("t=256 Simulations (Plain vs Stable PICS)") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig30_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig30.pdf"),
       height=9,
       width=9)

## |i-3| (spiked) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
supp_fig31_plot <- ggplot(all_pics2_scen_env_het_trend_df %>% subset(Scenario=="|i-3|"), 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoCausalVarSim~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#f89acb","#4d8dc3")) + 
  scale_colour_manual(values=c("#f781bf","#2171b5")) +
  ggtitle("|i-3| Simulations (Plain vs Stable PICS)") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig31_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig31.pdf"),
       height=9,
       width=9)

## i=3 (spiked) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
supp_fig32_plot <- ggplot(all_pics2_scen_env_het_trend_df %>% subset(Scenario=="i=3"), 
                          aes(x=Levels,y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Algorithm)) +
  geom_errorbar(aes(x=Levels,ymin=lwr.ci,ymax=upr.ci,colour=Algorithm),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(NoCausalVarSim~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#f89acb","#4d8dc3")) + 
  scale_colour_manual(values=c("#f781bf","#2171b5")) +
  ggtitle("i=3 Simulations (Plain vs Stable PICS)") +
  labs(fill="Approach",colour="Approach")+
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        legend.box = "vertical",
        legend.spacing.y = unit(0.03, "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig32_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig32.pdf"),
       height=9,
       width=9)

## Table1. PP-stratified matching frequencies of Plain vs Stable PICS ===========
# Done

## Table2&3. Matching frequencies in each SNR/S/posterior probability scenario ====
# Done --- PICS and SuSiE

## Table4. Impact of including more potential sets on matching frequency ========
# Include a table reporting the number of matching variants across potential sets
# and also the fraction of these diagonally matching variants that turn out to be
# causal.
## Stable PICS as Focal Algorithm --- Done
# All simulations
test<- getMatchStats(focal_df=all_stable_pics2_results,
                     other_df=all_top_pics2_results)
# > test %>% subset(int_focal_SNP1_other_SNP1!="NA") %>% select(PS1ContainsCausal) %>% table() 
# PS1ContainsCausal
# FALSE  TRUE 
# 400  1337 
# > test %>% subset(int_focal_SNP2_other_SNP2!="NA") %>% select(PS2ContainsCausal) %>% table() 
# PS2ContainsCausal
# FALSE  TRUE 
# 780   191 
# > test %>% subset(int_focal_SNP3_other_SNP3!="NA") %>% select(PS3ContainsCausal) %>% table() 
# PS3ContainsCausal
# FALSE  TRUE 
# 554    36 


# One causal variant
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_one_causal_results,
                                       other_df=top_pics2_one_causal_results)
stable_pics2_match_df %>% subset(!focal_SNP1_match_other_SNP1) %>%
  select(c("focal_SNP1_match_other_SNP2","PS1ContainsCausal")) %>%
  table()

# Two causal variants
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_two_causal_results,
                                       other_df=top_pics2_two_causal_results)
stable_pics2_match_df %>% subset(!focal_SNP1_match_other_SNP1) %>%
  select(c("focal_SNP1_match_other_SNP2","PS1ContainsCausal")) %>%
  table()

# Three causal variants
stable_pics2_match_df <- getMatchStats(focal_df=stable_pics2_three_causal_results,
                                       other_df=top_pics2_three_causal_results)
stable_pics2_match_df %>% subset(!focal_SNP1_match_other_SNP1) %>%
  select(c("focal_SNP1_match_other_SNP2","PS1ContainsCausal")) %>%
  table()

## Top PICS as Focal Algorithm 
# One causal variant
top_pics2_match_df <- getMatchStats(other_df=stable_pics2_one_causal_results,
                                    focal_df=top_pics2_one_causal_results)
top_pics2_match_df %>% subset(!focal_SNP1_match_other_SNP1) %>%
  select(c("focal_SNP1_match_other_SNP2","PS1ContainsCausal")) %>%
  table()

# Two causal variants
top_pics2_match_df <- getMatchStats(other_df=stable_pics2_two_causal_results,
                                    focal_df=top_pics2_two_causal_results)
top_pics2_match_df %>% subset(!focal_SNP1_match_other_SNP1) %>%
  select(c("focal_SNP1_match_other_SNP2","PS1ContainsCausal")) %>%
  table()

# Three causal variants
################################ OLD STUFF #####################################
## Fig 24-26: SNR-stratified results of Plain vs Stable PICS ===================
# This version shows causal variant recovery rate of each potential set separately. 
# The version we present in the paper shows the cumulative performance as more 
# potential sets are included. This is to be consistent with how we present results
# for all other algorithms in the Supplement (with exception of matching vs non-matching
# in view that those require potential set-specific discussions).
phi_vector <- c(0.05,0.1,0.2,0.4)
## One Causal Variant <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
big_env_het_1c_df <- NULL
# For each potential set
for (snp in 1:3) {
  # Add stable and plain PICS results
  big_env_het_stable_pics_one_causal_results$NoCausalRec <- 
    mapply(countIntersect,
           big_env_het_stable_pics_one_causal_results$CausalSNPs,
           big_env_het_stable_pics_one_causal_results[[paste0("SNP",snp)]])
  
  big_env_het_plain_pics_one_causal_results$NoCausalRec <- 
    mapply(countIntersect,
           big_env_het_plain_pics_one_causal_results$CausalSNPs,
           big_env_het_plain_pics_one_causal_results[[paste0("SNP",snp)]])
  # For each value of phi / SNR
  for (i in 1:4) {
    plain_sub <- big_env_het_plain_pics_one_causal_results %>% 
      subset(Phi == phi_vector[i])
    stable_sub <- big_env_het_stable_pics_one_causal_results %>% 
      subset(Phi == phi_vector[i])
    
    plain_ci_table <- DescTools::MultinomCI(table(plain_sub$NoCausalRec))
    plain_ci_table <- as.data.frame(plain_ci_table)
    plain_ci_table$Levels <- factor(as.numeric(rownames(plain_ci_table)))
    plain_ci_table$Phi <- rep(phi_vector[i], nrow(plain_ci_table))
    plain_ci_table$PS <- rep(paste0("Potential Set ", snp), 
                             nrow(plain_ci_table))
    plain_ci_table$Approach <- rep("Plain PICS", nrow(plain_ci_table))
    big_env_het_1c_df <- rbind(big_env_het_1c_df, 
                               plain_ci_table)
    
    stable_ci_table <- DescTools::MultinomCI(table(stable_sub$NoCausalRec))
    stable_ci_table <- as.data.frame(stable_ci_table)
    stable_ci_table$Levels <- factor(as.numeric(rownames(stable_ci_table)))
    stable_ci_table$Phi <- rep(phi_vector[i], nrow(stable_ci_table))
    stable_ci_table$PS <- rep(paste0("Potential Set ", snp), 
                              nrow(stable_ci_table))
    stable_ci_table$Approach <- rep("Stable PICS", nrow(stable_ci_table))
    big_env_het_1c_df <- rbind(big_env_het_1c_df, 
                               stable_ci_table)
  }
}

big_env_het_1c_df$SNR <- 
  paste0("SNR = ", round(big_env_het_1c_df$Phi/
                           (1-big_env_het_1c_df$Phi),3))

## Construct plot
supp_fig24_plot <- ggplot(big_env_het_1c_df, 
                          aes(x=factor(Levels,levels=c(0,1)),y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Approach)) +
  geom_errorbar(aes(x=factor(Levels,levels=c(0,1)),
                    ymin=lwr.ci,ymax=upr.ci,colour=Approach),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(PS~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#f89acb","#4d8dc3")) + 
  scale_colour_manual(values=c("#f781bf","#2171b5")) +
  ggtitle("Simulations with 1 Causal Variant\n(Plain vs Stable PICS)") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig24_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig24.pdf"),
       height=9,
       width=9)

## Two Causal Variants <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
big_env_het_2c_df <- NULL
# For each potential set
for (snp in 1:3) {
  # Add stable and plain PICS results
  big_env_het_stable_pics_two_causal_results$NoCausalRec <- 
    mapply(countIntersect,
           big_env_het_stable_pics_two_causal_results$CausalSNPs,
           big_env_het_stable_pics_two_causal_results[[paste0("SNP",snp)]])
  
  big_env_het_plain_pics_two_causal_results$NoCausalRec <- 
    mapply(countIntersect,
           big_env_het_plain_pics_two_causal_results$CausalSNPs,
           big_env_het_plain_pics_two_causal_results[[paste0("SNP",snp)]])
  # For each value of phi / SNR
  for (i in 1:4) {
    plain_sub <- big_env_het_plain_pics_two_causal_results %>% 
      subset(Phi == phi_vector[i])
    stable_sub <- big_env_het_stable_pics_two_causal_results %>% 
      subset(Phi == phi_vector[i])
    
    plain_ci_table <- DescTools::MultinomCI(table(plain_sub$NoCausalRec))
    plain_ci_table <- as.data.frame(plain_ci_table)
    plain_ci_table$Levels <- factor(as.numeric(rownames(plain_ci_table)))
    plain_ci_table$Phi <- rep(phi_vector[i], nrow(plain_ci_table))
    plain_ci_table$PS <- rep(paste0("Potential Set ", snp), 
                             nrow(plain_ci_table))
    plain_ci_table$Approach <- rep("Plain PICS", nrow(plain_ci_table))
    big_env_het_2c_df <- rbind(big_env_het_2c_df, 
                               plain_ci_table)
    
    stable_ci_table <- DescTools::MultinomCI(table(stable_sub$NoCausalRec))
    stable_ci_table <- as.data.frame(stable_ci_table)
    stable_ci_table$Levels <- factor(as.numeric(rownames(stable_ci_table)))
    stable_ci_table$Phi <- rep(phi_vector[i], nrow(stable_ci_table))
    stable_ci_table$PS <- rep(paste0("Potential Set ", snp), 
                              nrow(stable_ci_table))
    stable_ci_table$Approach <- rep("Stable PICS", nrow(stable_ci_table))
    big_env_het_2c_df <- rbind(big_env_het_2c_df, 
                               stable_ci_table)
  }
}

big_env_het_2c_df$SNR <- 
  paste0("SNR = ", round(big_env_het_2c_df$Phi/
                           (1-big_env_het_2c_df$Phi),3))

## Construct plot
supp_fig25_plot <- ggplot(big_env_het_2c_df, 
                          aes(x=factor(Levels,levels=c(0,1,2)),y=est)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           alpha=0.8,
           aes(fill=Approach)) +
  geom_errorbar(aes(x=factor(Levels,levels=c(0,1,2)),
                    ymin=lwr.ci,ymax=upr.ci,colour=Approach),
                width=0.2,alpha=0.9,position=position_dodge(width=0.85)) +
  facet_grid(PS~SNR) +
  ylim(c(0,1)) +
  ylab("Probability of Recovering N Variants") +
  xlab("No. Causal Variants, N") +
  scale_fill_manual(values=c("#f89acb","#4d8dc3")) + 
  scale_colour_manual(values=c("#f781bf","#2171b5")) +
  ggtitle("Simulations with 2 Causal Variants\n(Plain vs Stable PICS)") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,size=16),
        legend.position ='bottom',
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=11))

ggsave(supp_fig25_plot,
       file=paste0(plot_dir,"Sim_Supp_Fig25.pdf"),
       height=9,
       width=9)

## Three Causal Variants <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
big_plain_pics_df 
big_stable_pics_df 

