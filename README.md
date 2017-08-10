# vanddraabe: Identification and Statistical Analysis of Conserved Waters in Proteins
[![Build Status](https://travis-ci.org/exeResearch/vanddraabe.svg?branch=master)](https://travis-ci.org/exeResearch/vanddraabe)

`vanddraabe` provides a powerful way to identify and analyze conserved waters within crystallographic protein structures and molecular dynamics simulation trajectories. Statistical parameters for each water cluster, informative graphs, and a PyMOL session file to visually explore the conserved waters and protein are returned. Hydrophilicity is the propensity of waters to congregate near specific protein atoms and is related to conserved waters. An informatics derived set of hydrophilicity values are provided based on a large, high-quality X-ray protein structure dataset.

This package includes function to:

* Perform conserved water analysis on RCSB files and molecular dynamics simulations
* Create analysis plots
* Create PyMOL session scripts to visualize the conserved waters


## Installing vanddraabe

`vanddraabe` is available here (on GitHub) and on [CRAN](https://cran.r-project.org/package=vanddraabe). To install it:

```r
# The easiest way to get vanddraabe is:
install.packages("vanddraabe")

# Or get the development version from GitHub:
install.packages("devtools")
devtools::install_github("tidyverse/vanddraabe")
```

## How to use vanddraabe

The vignette provides a detailed example of using `vanddraabe` to identify the conserved waters of Thrombin.


## Need help? Found a bug?

Please submit an [issue via GitHub](https://github.com/exeResearch/vanddraabe/issues).

