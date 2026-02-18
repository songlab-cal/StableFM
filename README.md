# Stability and Genetic Fine-Mapping

This repository contains data and scripts for the following paper:

> Aw, A.J., Jin, L.C., Ioannidis, N.M., and Song, Y.S. (2026) "[The impact of stability considerations on genetic fine-mapping](https://elifesciences.org/reviewed-preprints/88039)" *eLife* RP88039

## Fine-mapping results

Results are stored under `data/results_with_moderators`.

For convenience, csv files with both raw and Benjamini-Hochberg-adjusted p-values are also stored under `data/results_with_moderators/p_values`.

**Update (Oct 5, 2025)** Code for generating plots involving simulations have been added under `scripts/simulations` subdirectory. 
- Two files must be unzipped in `data/simulations` before running plot generating code (`gen_sim_study_figures.R`) 
- Scripts for simulating gene expression are compressed as `scripts/simulations/gene-exp-sim-scripts.zip`
- Scripts for running all fine-mapping algorithms are compressed as `scripts/simulations/fine-mapping-scripts.zip`

## Summary data files

Summary files containing the 378 functional annotations associated with each fine-mapped variant (stable or top) are stored under `data/all_func_annots`. We recommend directly working with these files to reproduce our findings and figures.

## Generating figures in paper

Code for generating figures is stored under `scripts`.

```
# Command for Executing Script
Rscript generate_paper_figs.R
```

## Interactive application

We provide a Shiny application for visualizing and interpreting our results at both single-gene and pan-gene levels. You may explore the app by clicking on the link below. 

> [https://alan-aw.shinyapps.io/stability_v1/](https://alan-aw.shinyapps.io/stability_v1/)

Open-source code for building, deploying and modifying our app for your own use is available under `ShinyApp`. 

## Running PICS2

Standalone code for running PICS2 is stored under `scripts`, along with example data. The example data uses gene expression of [PRKAA2](http://genome.cse.ucsc.edu/cgi-bin/hgGene?org=Human&hgg_chrom=none&hgg_type=knownGene&hgg_gene=uc001cyk.5) (ENSG00000162409.6), which is transcribed in Chromosome 1. 

```
# Command for Executing Script - RUN AFTER UNZIPPING toy-data/genotype_array.csv.zip
# Internal functions can also be copied and used elsewhere; see R script for details
Rscript runPICS2.R 500 3 toy-data/genotype_array.csv toy-data/phenotype_vector.txt toy-data/snp_metadata.txt toy-data/sample_metadata.txt toy-data/gene_tss_metadata.txt toy-data/covars_array.csv toy-results/
```