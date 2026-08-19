# BaP-Pb-additivity
Testing Additivity of Lead and Benzo[a]pyrene-induced Neurotoxicity in *Caenorhabditis elegans* Assays

This repository contains the code and tutorial for reproducing the main results from:

> [Niccolò Anceschi, Huayta, Meyer, Dunson and Herring (2026+). *Testing Additivity of Lead and Benzo[a]pyrene-induced Neurotoxicity in Caenorhabditis elegans Assays.*](https://arxiv.org/pdf/PLACEHOLDER)

## Overview

The repo implements the BaP-Pb drug-additivity analysis presented in the manuscript. The analysis uses a convex-mixture model $F(\cdot \mid d)=(1-\lambda(d))F_0+\lambda(d)F_\infty$ for ordinal damage-score distributions.
Here, $F_0$ and $F_\infty$ represent the unexposed and maximal-effect profiles, respectively, and the dose dependence is captured by the monotone weight function $\lambda(d)$.
For joint BaP-Pb exposure, the weight becomes a surface $\lambda_{12}(d_1,d_2)$ over the two-dimensional dose plane. 
Additivity is assessed by comparing the fitted joint radial-spline model against three additive null models, **Bliss**, **Loewe+** (EME), and **Hand**.
Each null model induces an additive baseline surface $\lambda_{12}^{(o)}(d_1,d_2)$ based on the corresponding single-drug fits.

The accompanying tutorial [`Analysis.hmtl`](https://niccoloanceschi.github.io/BaP-Pb-additivity/Analysis.html) walks through the complete analysis pipeline, from importing the pre-processed data and fitting the models to producing the figures reported in the manuscript and performing the bootstrap additivity tests.

## Repository structure

The main analysis is organized into the following scripts:

| File | Description |
|---|---|
| `Analysis.Rmd` | Reproducible tutorial generating the main results and figures |
| `Utils.R` | Generic utilities, including softmax, mixture probabilities, and monotone inversion |
| `Data_setup.R` | Data import and construction of spline bases, dose constants, extremal profiles, and analysis subsets |
| `Loglikelihood.R` | Single- and joint-exposure likelihoods and penalized radial-spline objectives |
| `Null_models.R` | Construction of the Bliss, Loewe+ and Hand additive null models |
| `Model_fit.R` | Fitting of the single-drug models and joint radial-spline surface |
| `Plots.R` | Functions for generating the main analysis figures |
| `Bootstrap_test.R` | Parametric conditional bootstrap procedure for testing additivity |
| `C-elgans_preprocessed_data.rds` | Pre-processed *C. elegans* assay data |

## Requirements

The analysis is implemented in **R** and uses the following packages:

```r
library(splines2)   # iSpline, bSpline
library(subplex)    # subplex() optimizer
library(ggplot2)    # plotting
library(tidyr)      # reshaping dataframes
library(tictoc)     # runtimes
library(patchwork)  # composite plots
library(here)       # locating source code
library(knitr)      # R Markdown
library(rstudioapi) # locating the active R Markdown document
```

These packages can be installed from CRAN with:

```r
install.packages(c(
  "splines2",
  "subplex",
  "ggplot2",
  "tidyr",
  "tictoc",
  "patchwork",
  "here",
  "knitr",
  "rstudioapi"
))
```

## Reproducing the analysis

The easiest way to reproduce the analysis is to open `Tutorial.Rmd` in RStudio and run in its entirety.
When the tutorial is opened in RStudio, the repository path is detected automatically from the file's location, so no manual editing of the repository path is needed:

```r
source_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
data_file <- "C-elgans_preprocessed_data.rds"
```

The tutorial sources the analysis scripts in dependency order:

1. **Utilities** - generic mathematical and computational helpers.
2. **Data setup** - imports the pre-processed data and constructs the objects required by the models.
3. **Likelihoods and radial splines** - defines the likelihood functions and penalized radial-spline objective.
4. **Additive null models** - constructs the Bliss, Loewe+ and Hand baseline surfaces.
5. **Model fitting** - fits the single-drug models and joint-exposure radial-spline surface.
6. **Visualization** - generates the figures used to examine the fitted models and compare likelihoods.

Each script contains dedicated documentation describing its contents.

## Citation

If you use this code or the associated methodology, please cite:

> Anceschi, N., Huayta, M., Meyer, D., Dunson, D. and Herring, A. (2026+).  
> *Testing Additivity of Lead and Benzo[a]pyrene-induced Neurotoxicity in Caenorhabditis elegans Assays.*

See the [arXiv manuscript](https://arxiv.org/pdf/PLACEHOLDER) for the methodological details and full description of the analysis.