# gnomAD high allele balance error mode

This repo contains code for the analysis and figures in [our recent preprint](https://www.biorxiv.org/content/10.1101/784157v1) that identifies an error mode 
in gnomAD relating to individuals with high allele balance (fraction of alternate reads) being called as heterozygotes.

### Analysis code

`find_high_ad_hets.py` requires Hail 0.2, and only the first section (`--find_high_ad_hets`) can be run without access to the restricted individual-level dataset.

The remainder of the analysis code is provided as-is for reference.

### Plotting code

`plot_high_ad_het_data.R` includes the code to parse, process, and plot the resulting data, in order to create the 3 figures in the manuscript.
