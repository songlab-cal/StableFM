# Processing Enformer Output

This subdirectory summarizes key steps to generating Enformer track annotations for our analyses. We begin with track predictions output from the Enformer. These come in _pairs_ --- one for the reference sequence and another for the perturbed sequence. As described in our manuscript, functional impact is determined by the magnitude of change in these track predictions.

Below, we break down our data processing into modular steps to facilitate exposition, at the cost of small efficiency losses. For consistency, all auxiliary or intermediate output files are located in the `enformer-data` subdirectory. (The final outputs are in the `data` subdirectory of the repo.)

## 1. Download and process track predictions

The compressed object containing track predictions is available [here](https://doi.org/10.6084/m9.figshare.22032167.v1) (note file size is slightly larger than 6GB). 

  1.1 Download the object and decompress to retrieve eight numpy files (e.g., run `unzip raw_scores.zip`) 
  1.2 Process numpy files (run `python get_csv.py`)

  _Explanation._ This step converts the numpy files to eight corresponding csv files that are easier for subsetting in the next step. The csv files contain predictions across 5,313 tracks (listed [here](https://raw.githubusercontent.com/calico/basenji/master/manuscripts/cross2020/targets_human.txt)), for a few thousand genes. 

## 2. Subset to relevant tracks

Enformer performs predictions on 5,313 tracks, but our fine-mapping experiments apply only to the lymphoblastoid cell line (GM12878), a small subset of 177 tracks. We have provided the shortlist of relevant tracks (see `enformer-data/basenji_relevant_metadata.csv`), obtained from the original [longlist](https://raw.githubusercontent.com/calico/basenji/master/manuscripts/cross2020/targets_human.txt). 

  2.1 Subset the csvs (run `python subset_predictions.py`)

## 3. Combine and compute perturbation scores

The eight numpy files are combined by grouping stable variants together and top variants together. Moreover, we convert differences in track predictions (ALT - REF) to absolute differences (|ALT - REF|), or ''perturbation scores.''

  3.1 Combine and convert differences (run `get_magnitudes.R`)

The final outputs from accomplishing the steps above are:

  - (Averaged absolute differences for stable variants) `stab_diff_ave_relevant_comb.rds` 
  - (TSS-only absolute differences for stable variants) `stab_diff_tss_relevant_comb.rds` 
  - (Averaged absolute differences for top variants) `top_diff_ave_relevant_comb.rds`
  - (TSS-only absolute differences for top variants) `top_diff_tss_relevant_comb.rds`

These objects are located under the data subdirectory of this Github repo, and can be loaded into R for statistical analyses. 