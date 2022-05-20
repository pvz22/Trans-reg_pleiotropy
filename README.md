# Trans-reg_pleiotropy

This repository includes data and all code necessary to recreate analysis and figures as published in:
"Pleiotropic effects of trans-regulatory mutations on gene expression and fitness"

Both raw and processed sequencing data are located at GEO accession #: GSE175398

A release of the repository contents (v1.0) at the time of publication is available at zenodo, under DOI 10.5281/zenodo.6567260.

Files created in "Growthrate.R" and "DESeqProcessing.R" are used in "Peio_cis_trans_analysis.R".

If you would like to begin at "Pleio_cis_trans_analysis.R", simply use the DEseq2 output files for each strain that are available at GEO.

To start at "DESeqProcessing.R", use the salmon output files available at GEO.

To re-do the salmon quantification from the raw data, use "cutmapsort.pbs" and then "salmon.pbs" along with the raw sequencing data located at GEO.
