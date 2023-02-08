# p-value Files

Under this directory, there are 9 csv files, each reporting p-values across all 378 annotations and all 3 potential sets.

- (Matching vs Non-matching Top Variants) `match_vs_top.csv`
- (Matching vs Non-matching Stable Variants) `match_vs_stable.csv`
- (Top vs Stable Variants) `top_vs_stable.csv`
- (Trend Analysis - Degree of Stability) `mod_Slice Stability.csv`
- (Trend Analysis - Population Diversity) `mod_Maximum AF Diff.csv`
- (Trend Analysis - Population Differentiation) `mod_Maximum FST.csv`
- (Trend Analysis - Inclusion of Distal Subpopulations for Top Variant) `mod_Top SNV Discovered in YRI.csv`
- (Trend Analysis - Inclusion of Distal Subpopulations for Stable Variant) `mod_Stable SNV Discovered in YRI.csv`
- (Trend Analysis - Degree of Certainty of Causality) `mod_Posterior Probability of Top SNV.csv`

## Conditional Test p-values

We also provide p-values for the conditional analyses reported in our paper. There are four such analyses in total, for which we have assigned four subdirectories. Within each subdirectory, there are 7 csv files. These are the same as the 9 enumerated above, except for the files `match_vs_stable.csv` and `match_vs_top.csv`, which we did not test. 


| Conditional Analysis Type                                                                          | Subdirectory              |
| -------------------------------------------------------------------------------------------------- | ------------------------- |
| Restriction to genes for which the top variant’s posterior probability exceeds 0.9                 | `post_prob_topsnp_0.9`  |
| Restriction to genes for which the stable variant’s posterior probability exceeds 0.9              | `post_prob_stabsnp_0.9` |
| Restriction to genes for which the top variant’s positive posterior probability support exceeds 10 | `post_prob_support_10`  |
| Restriction to genes for which the top variant’s positive posterior probability support exceeds 50 | `post_prob_support_50`  |