# thesis

## *Comparison of statistical methods that account for unmeasured confounding in cell composition data*

This repository contains the R source code for the master's thesis: *Comparison of statistical methods that account for unmeasured confounding in cell composition data*. This thesis was written in fulfillment of the requirements for the degree of Master of Science in Statistics and Data Science at **KU Leuven**, in collaboration with **Johnson & Johnson**.

**Author:** Sarah Brosens  
**Supervisor:** Prof. A. A. Abad (KU Leuven)  
**Co-supervisor:** K. Van den Berge (Johnson & Johnson)   
**Mentor:** A. Huybrechts (Ghent University, Johnson & Johnson)  
**Academic year:** 2025-2026

## Structure

**Analysis code:**
* `01_preprocessing.R`: preprocessing of raw cell composition data into original RA data
* `02_downsampling.R`: downsampling of original RA data to create imbalanced RA data
* `03_EDA.R`: exploratory data analysis on original and imbalanced RA data
* `04_case_study.R`: cell composition analysis on original and imbalanced RA data
* `05_simulation_study.R`: simulation study on simulated RA data
* `06_robustness_test.R`: robustness test against violation of sparsity assumption on simulated RA data

**Helper functions:**
* `parametric_simulation_function.R`: parametric simulation function to simulate cell type counts
* `performance_evaluation_functions.R`: performance evaluation functions to investigate methods for bias correction in cell composition analysis

## Data availability
Cell composition data were taken from a prospective longitudinal study on Rheumatoid Arthritis (RA) [1], and are publicly available through the [Allen Institute of Immunology](https://apps.allenimmunology.org/aifi/insights/ra-progression/downloads/scrna/). You can find these data in the `data` folder.

## Prerequisites
To run this code, you need to have **R** installed, alongside several required packages. Most important is the installation of the `voomCLR` package [2], which is currently available through [GitHub](https://github.com/johnsonandjohnson/voomCLR):

```R
devtools::install_github("johnsonandjohnson/voomCLR")
``` 

## Usage
To reproduce the findings of this master's thesis, run the scripts sequentially from `01` to `06`. Ensure your working directory is set to the root of this repository first.

## References

1. He Z, Glass MC, Venkatesan P, Feser ML, Lazaro L, Okada LY, et al. Progression to rheumatoid arthritis in at-risk individuals is defined by systemic inflammation and by T and B cell dysregulation. *Sci Transl Med.* 2025 Sep 24;17(817):eadt7214.</small> 

2. Takele Assefa A, Verbist B, Van den Berge K. Assessing differential cell composition in single-cell studies using voomCLR. *Bioinformatics.* 2026 Jan 2;42(1):btaf637.

## Copyright
© Copyright by KU Leuven  

Without prior written permission from both the supervisor(s) and the author(s), copying, reproducing, using, or realizing this publication or parts thereof is prohibited. For requests or information regarding the copying and/or use and/or realization of parts of this publication, please contact KU Leuven, Faculty of Science, Celestijnenlaan 200H - box 2100, 3001 Leuven (Heverlee), Telephone +32 16 32 14 01.
