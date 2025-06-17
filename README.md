# HTE‑ROSE

## Overview

This repository contains the code needed to reproduce the Heterogeneity of Treatment Effect (HTE) of Neuromuscular Blocking Agents (NMBAs) on 90-day mortality according to subphenotypes determined by clinical variables in adult patients with Acute Respiratory Distress Syndrome (ARDS) using the data from the Reevaluation of Systemic Early Neuromuscular Blockade (ROSE) clinical trial. The analysis is fully scripted in R and relies on Bayesian hierarchical logistic model. It contains two R scripts:

## Repository structure

```text
.
├── R/                      
│   ├── hterose_model_interaction.R   
│   └── hterose_model_mort.R          
├── data/
│   └── rose_subset.csv # Input dataset - not included for licensing reasons
├── renv/ # renv infrastructure for reproducible environments
├── renv.lock # Exact package versions
├── hte_rose.Rproj # RStudio project file
└── README.md
```
## Dependencies

* **R** ≥ 4.0

`renv.lock` records the complete set, but the main runtime packages are:

* **tidyverse** (data wrangling & plotting -> dplyr, ggplot2, tibble)
* **brms**, **beanz** (Bayesian HTE analysis models)
* **ggpmisc**, **gridExtra** (plot annotation)


## Quick start

```r
# 1. Clone the repo
git clone https://github.com/cincythinklab/hte_rose.git
cd hte_rose

# 2. Open in R (or RStudio) and restore the environment
install.packages("renv")      # if you do not have renv
renv::restore()               # installs the exact package versions

# 3. Place the trial extract in the data/ folder
## data/rose_subset.csv   - not included for licensing reasons

# 4. Run an analysis
source("R/hterose_model_mort.R")          
source("R/hterose_model_interaction.R")
# Or from the command line:
Rscript R/hterose_model_interaction.R
Rscript R/hterose_model_mort.R

# 5. Outputs
Figures and tables will be rendered in your R session (or saved if you modify the script to write to disk).
```






