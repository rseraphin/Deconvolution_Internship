# Deconvolution Internship

This reporory contains work done by Remi SERAPHIN during his internship from 02/2020 to 06/2020.

The objective is to infer cell composition of a bulk dataset using cell type specific expression derived using single cell dataset using deconvolution.

## Requirements

To perform data analysis R is required.

To perform deconvolution Python is required.

To perform the example Jupyter is required.

## Examples

Example of the Seurat (R) workflow to analyse data is provided in the 'Workflow_Seurat.R' file.

Example and functions for the deconvolution is provided in the 'Function_and_Examples.ipynb' notebook.

## Observations

The first simulation method "rand_bulk" and associated "rand_subset" doesn't yield good result and should be improved.

The "nnls" method seem to be the best tested deconvolution method as it is faster and more precise. It is thus the best lead to persue this work.
