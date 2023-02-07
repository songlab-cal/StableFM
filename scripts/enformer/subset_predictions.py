## Script for subsetting all track predictions to 177 LCL-relevant ones

import numpy as np
import pandas as pd

# Load all relevant CSV files
# [!] Files here not provided due to large size --- requires previous step to be run to work
relevant_metadata=pd.read_csv('enformer-data/basenji_relevant_metadata.csv')
stab_diff_pred_ave=pd.read_csv('enformer-data/stab_diff_pred_ave.csv')
stab_diff_pred_tss=pd.read_csv('enformer-data/stab_diff_pred_tss.csv')
top_diff_pred_ave=pd.read_csv('enformer-data/top_diff_pred_ave.csv')
top_diff_pred_tss=pd.read_csv('enformer-data/top_diff_pred_tss.csv')
stab_diff_pred_ave_new=pd.read_csv('enformer-data/stab_diff_pred_ave_new.csv')
stab_diff_pred_tss_new=pd.read_csv('enformer-data/stab_diff_pred_tss_new.csv')
top_diff_pred_ave_new=pd.read_csv('enformer-data/top_diff_pred_ave_new.csv')
top_diff_pred_tss_new=pd.read_csv('enformer-data/top_diff_pred_tss_new.csv')

# AAverage prediction, stab 
stab_diff_pred_ave=stab_diff_pred_ave.to_numpy()
stab_diff_pred_ave=np.transpose(stab_diff_pred_ave)
stab_diff_pred_ave_relevant=pd.DataFrame(np.transpose(stab_diff_pred_ave[(relevant_metadata['index'].to_numpy()+1),:]))
stab_diff_pred_ave_relevant.columns=relevant_metadata['identifier']
stab_diff_pred_ave_relevant.to_csv('enformer-data/stab_diff_pred_ave_relevant.csv',index=False)

stab_diff_pred_ave_new=stab_diff_pred_ave_new.to_numpy()
stab_diff_pred_ave_new=np.transpose(stab_diff_pred_ave_new)
stab_diff_pred_ave_new_relevant=pd.DataFrame(np.transpose(stab_diff_pred_ave_new[(relevant_metadata['index'].to_numpy()+1),:]))
stab_diff_pred_ave_new_relevant.columns=relevant_metadata['identifier']
stab_diff_pred_ave_new_relevant.to_csv('enformer-data/stab_diff_pred_ave_relevant_new.csv',index=False)

# Average prediction, top 
top_diff_pred_ave=top_diff_pred_ave.to_numpy()
top_diff_pred_ave=np.transpose(top_diff_pred_ave)
top_diff_pred_ave_relevant=pd.DataFrame(np.transpose(top_diff_pred_ave[(relevant_metadata['index'].to_numpy()+1),:]))
top_diff_pred_ave_relevant.columns=relevant_metadata['identifier']
top_diff_pred_ave_relevant.to_csv('enformer-data/top_diff_pred_ave_relevant.csv',index=False)

top_diff_pred_ave_new=top_diff_pred_ave_new.to_numpy()
top_diff_pred_ave_new=np.transpose(top_diff_pred_ave_new)
top_diff_pred_ave_new_relevant=pd.DataFrame(np.transpose(top_diff_pred_ave_new[(relevant_metadata['index'].to_numpy()+1),:]))
top_diff_pred_ave_new_relevant.columns=relevant_metadata['identifier']
top_diff_pred_ave_new_relevant.to_csv('enformer-data/top_diff_pred_ave_relevant_new.csv',index=False)

# TSS-only prediction, stab 
stab_diff_pred_tss=stab_diff_pred_tss.to_numpy()
stab_diff_pred_tss=np.transpose(stab_diff_pred_tss)
stab_diff_pred_tss_relevant=pd.DataFrame(np.transpose(stab_diff_pred_tss[(relevant_metadata['index'].to_numpy()+1),:]))
stab_diff_pred_tss_relevant.columns=relevant_metadata['identifier']
stab_diff_pred_tss_relevant.to_csv('enformer-data/stab_diff_pred_tss_relevant.csv',index=False)

stab_diff_pred_tss_new=stab_diff_pred_tss_new.to_numpy()
stab_diff_pred_tss_new=np.transpose(stab_diff_pred_tss_new)
stab_diff_pred_tss_new_relevant=pd.DataFrame(np.transpose(stab_diff_pred_tss_new[(relevant_metadata['index'].to_numpy()+1),:]))
stab_diff_pred_tss_new_relevant.columns=relevant_metadata['identifier']
stab_diff_pred_tss_new_relevant.to_csv('enformer-data/stab_diff_pred_tss_relevant_new.csv',index=False)

# TSS-only, top
top_diff_pred_tss=top_diff_pred_tss.to_numpy()
top_diff_pred_tss=np.transpose(top_diff_pred_tss)
top_diff_pred_tss_relevant=pd.DataFrame(np.transpose(top_diff_pred_tss[(relevant_metadata['index'].to_numpy()+1),:]))
top_diff_pred_tss_relevant.columns=relevant_metadata['identifier']
top_diff_pred_tss_relevant.to_csv('enformer-data/top_diff_pred_tss_relevant.csv',index=False)

top_diff_pred_tss_new=top_diff_pred_tss_new.to_numpy()
top_diff_pred_tss_new=np.transpose(top_diff_pred_tss_new)
top_diff_pred_tss_new_relevant=pd.DataFrame(np.transpose(top_diff_pred_tss_new[(relevant_metadata['index'].to_numpy()+1),:]))
top_diff_pred_tss_new_relevant.columns=relevant_metadata['identifier']
top_diff_pred_tss_new_relevant.to_csv('enformer-data/top_diff_pred_tss_relevant_new.csv',index=False)
