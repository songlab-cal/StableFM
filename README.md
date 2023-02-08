# Stability and Fine-Mapping

This repository contains data and scripts for the following paper:

> Aw, A.J., Jin, L.C., Ioannidis, N.M., and Song, Y.S. (2023+) "The impact of stability considerations on genetic fine-mapping"

## Fine-mapping results

Results are stored under `data/results_with_moderators`.

For convenience, csv files with Benjamini-Hochberg-adjusted p-values are also stored under `data/results_with_moderators/p_values`.

## Summary data files

Summary files containing the 378 functional annotations associated with each fine-mapped variant (stable or top) are stored under `data/all_func_annots`. We recommend directly working with these files to reproduce our findings and figures.

## Generating figures in paper

Code for generating figures is stored under `scripts`.

```
# Command for Executing Script
Rscript generate_paper_figs.R
```

## Shiny App

Open-source code for building and deploying our Shiny application is under `ShinyApp`. 

## Running PICS2

Standalone code for running PICS2 is stored under `scripts`, along with example data. The example data uses gene expression of [PRKAA2](http://genome.cse.ucsc.edu/cgi-bin/hgGene?org=Human&hgg_chrom=none&hgg_type=knownGene&hgg_gene=uc001cyk.5) (ENSG00000162409.6), which is transcribed in Chromosome 1. 

```
# Command for Executing Script - RUN AFTER UNZIPPING toy-data/genotype_array.csv.zip
# Internal functions can also be copied and used elsewhere; see R script for details
Rscript runPICS2.R 500 3 toy-data/genotype_array.csv toy-data/phenotype_vector.txt toy-data/snp_metadata.txt toy-data/sample_metadata.txt toy-data/gene_tss_metadata.txt toy-data/covars_array.csv toy-results/
```