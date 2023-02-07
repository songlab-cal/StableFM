## Script for converting numpy files to csv

import numpy as np
import pandas as pd

# Load eight numpy files
top_ref_pred=np.load('final_ref_preds_top_snps_strand.npy')
top_alt_pred=np.load('final_alt_preds_top_snps_strand.npy')

stab_ref_pred=np.load('final_ref_preds_stab_snps_strand.npy')
stab_alt_pred=np.load('final_alt_preds_stab_snps_strand.npy')

top_ref_pred_new=np.load('final_ref_preds_top_snps_new_strand.npy')
top_alt_pred_new=np.load('final_alt_preds_top_snps_new_strand.npy')

stab_ref_pred_new=np.load('final_ref_preds_stab_snps_new_strand.npy')
stab_alt_pred_new=np.load('final_alt_preds_stab_snps_new_strand.npy')

# Compute TSS-only track predictions
top_ref_pred_tss=pd.DataFrame(top_ref_pred[:,1,:])
top_alt_pred_tss=pd.DataFrame(top_alt_pred[:,1,:])
top_diff_pred_tss=top_alt_pred_tss-top_ref_pred_tss
top_diff_pred_tss.to_csv('top_diff_pred_tss.csv')

top_ref_pred_tss_new=pd.DataFrame(top_ref_pred_new[:,1,:])
top_alt_pred_tss_new=pd.DataFrame(top_alt_pred_new[:,1,:])
top_diff_pred_tss_new=top_alt_pred_tss_new-top_ref_pred_tss_new
top_diff_pred_tss_new.to_csv('top_diff_pred_tss_new.csv')

stab_ref_pred_tss=pd.DataFrame(stab_ref_pred[:,1,:])
stab_alt_pred_tss=pd.DataFrame(stab_alt_pred[:,1,:])
stab_diff_pred_tss=stab_alt_pred_tss-stab_ref_pred_tss
stab_diff_pred_tss.to_csv('stab_diff_pred_tss.csv')

stab_ref_pred_tss_new=pd.DataFrame(stab_ref_pred_new[:,1,:])
stab_alt_pred_tss_new=pd.DataFrame(stab_alt_pred_new[:,1,:])
stab_diff_pred_tss_new=stab_alt_pred_tss_new-stab_ref_pred_tss_new
stab_diff_pred_tss_new.to_csv('stab_diff_pred_tss_new.csv')

# Compute predictions fromm average of TSS and two sides of TSS
top_ref_pred_ave=pd.DataFrame(np.mean(top_ref_pred,axis=1))
top_alt_pred_ave=pd.DataFrame(np.mean(top_alt_pred,axis=1))
top_diff_pred_ave=top_alt_pred_ave-top_ref_pred_ave
top_diff_pred_ave.to_csv('top_diff_pred_ave.csv')

top_ref_pred_ave_new=pd.DataFrame(np.mean(top_ref_pred_new,axis=1))
top_alt_pred_ave_new=pd.DataFrame(np.mean(top_alt_pred_new,axis=1))
top_diff_pred_ave_new=top_alt_pred_ave_new-top_ref_pred_ave_new
top_diff_pred_ave_new.to_csv('top_diff_pred_ave_new.csv')

stab_ref_pred_ave=pd.DataFrame(np.mean(stab_ref_pred,axis=1))
stab_alt_pred_ave=pd.DataFrame(np.mean(stab_alt_pred,axis=1))
stab_diff_pred_ave=stab_alt_pred_ave-stab_ref_pred_ave
stab_diff_pred_ave.to_csv('stab_diff_pred_ave.csv')

stab_ref_pred_ave_new=pd.DataFrame(np.mean(stab_ref_pred_new,axis=1))
stab_alt_pred_ave_new=pd.DataFrame(np.mean(stab_alt_pred_new,axis=1))
stab_diff_pred_ave_new=stab_alt_pred_ave_new-stab_ref_pred_ave_new
stab_diff_pred_ave_new.to_csv('stab_diff_pred_ave_new.csv')