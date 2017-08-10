
<!-- README.md is generated from README.Rmd. Please edit that file -->
vanddraabe: Identification and Statistical Analysis of Conserved Waters in Proteins
===================================================================================

[![Travis-CI Build Status](https://travis-ci.org/exeResearch/vanddraabe.svg?branch=master)](https://travis-ci.org/exeResearch/vanddraabe) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/exeResearch/vanddraabe?branch=master&svg=true)](https://ci.appveyor.com/project/exeResearch/vanddraabe) [![CRAN\_Release\_Badge](http://www.r-pkg.org/badges/version-ago/vanddraabe)](https://CRAN.R-project.org/package=vanddraabe) [![CRAN\_Download\_Badge](http://cranlogs.r-pkg.org/badges/vanddraabe)](https://CRAN.R-project.org/package=vanddraabe)

`vanddraabe` provides a powerful way to identify and analyze conserved waters within crystallographic protein structures and molecular dynamics simulation trajectories. Statistical parameters for each water cluster, informative graphs, and a [PyMOL](http://www.pymol.org) session file to visually explore the conserved waters and protein are returned. Hydrophilicity is the propensity of waters to congregate near specific protein atoms and is related to conserved waters. An informatics derived set of hydrophilicity values are provided based on a large, high-quality X-ray protein structure dataset.

This package is a reimplementation and expansion of the WatCH<sup>1</sup> and PyWATER<sup>2</sup> applications and was created to provide the following abilities:

-   Perform conserved water analysis on RCSB files and molecular dynamics simulations
-   Provide access to data at each step of the analysis
-   Provide detailed statistical summaries for all waters being analyzed
-   Ability to analyze more than 65,000 waters
-   Create preformatted analysis plots
-   Create [PyMOL](http://www.pymol.org) session scripts to visualize the conserved waters
-   Write out an Excel workbook of initial, intermediate, and final results
-   Perform protein alignment using bio3d ([CRAN](https://cran.r-project.org/package=bio3d), [website](http://thegrantlab.org/bio3d/), and [BitBucket](https://bitbucket.org/Grantlab/bio3d)), an opensource applicaiton

1.  Paul C Sanschagrin and Leslie A Kuhn. Cluster analysis of consensus water sites in thrombin and trypsin shows conservation between serine proteases and contributions to ligand specificity. *Protein Science*, 1998, **7** (*10*), pp 2054-2064.
    [DOI: 10.1002/pro.5560071002](http://doi.org/10.1002/pro.5560071002)
    [PMID: 9792092](http://www.ncbi.nlm.nih.gov/pubmed/9792092)
    [WatCH webpage](http://www.kuhnlab.bmb.msu.edu/software/watch/index.html)

2.  Hitesh Patel, Bjorn A. Gruning, Stefan Gunther, and Irmgard Merfort. PyWATER: a PyMOL plug-in to find conserved water molecules in proteins by clustering. *Bioinformatics*, 2014, **30** (*20*), pp 2978-2980.
    [DOI: 10.1093/bioinformatics/btu424](http://doi.org/10.1093/bioinformatics/btu424)
    [PMID: 24990608](http://www.ncbi.nlm.nih.gov/pubmed/24990608)
    [PyWATER on GitHub](https://github.com/hiteshpatel379/PyWATER/blob/master/README.rst)

Installing vanddraabe
---------------------

`vanddraabe` is available on [GitHub](https://github.com/exeResearch/vanddraabe) and on [CRAN](https://cran.r-project.org/package=vanddraabe). To install it:

``` r
# The easiest way to get vanddraabe is:
install.packages("vanddraabe")

# Or get the development version from GitHub:
# install.packages("devtools")
devtools::install_github("exeResearch/vanddraabe")
```

How to use vanddraabe
---------------------

The vignette provided [here](http://www.exeresearch.com/vanddraabe.html) is a detailed example of using `vanddraabe` to identify the conserved waters of ten Thrombin structures.

Have a suggestion? Need help? Found a bug?
------------------------------------------

-   Contact [Emilio](https://github.com/emilioxavier) at <emilio@exeResearch.com>
-   Submit an [issue via GitHub](https://github.com/exeResearch/vanddraabe/issues)

Code of conduct
---------------

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
