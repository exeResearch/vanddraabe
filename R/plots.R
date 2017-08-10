## plots.R
##
## apr-10-2016 (exe) created
## feb-14-2017 (exe) added b-value, occupany, normalized b-value, and mobility
##                   summary histograms
## mar-16-2017 (exe) updated the documentation (rMarkdown)
## apr-12-2017 (exe) updated plotting functions to take results from ConservedWaters()
## apr-12-2017 (exe) added multi-page facet plot examples
## apr-13-2017 (exe) updated plotting functions to take results from CleanProteinStructures()
## apr-16-2017 (exe) updated documentation of ClusterSummaryPlots() sub-plots
## jul-07-2017 (exe) updated functions to accept ConservedWaters() output
## jul-25-2017 (exe) updated documentation
## jul-29-2017 (exe) MobNormBvalEvalPlots() uses multiple facets to create plot
## jul-31-2017 (exe) updated *.summ() plot examples for ggforce
## jul-31-2017 (exe) updated *.summ() plot examples with \dontrun{}
## jul-31-2017 (exe) updated CleanProteinStructure() -> CleanProteinStructures() in documentation
## jul-31-2017 (exe) updated ClusterSummaryPlots() documentation
## aug-04-2017 (exe) added protein mobility to BoundWaterEnvSummaryPlot()
## aug-04-2017 (exe) removed protein hydrophilicity from BoundWaterEnvSummaryPlot()
##
## Please direct all questions to Emilio Xavier Esposito, PhD
## exeResearch LLC, East Lansing, Michigan 48823 USA
## http://www.exeResearch.com
## emilio AT exeResearch DOT com
## emilio DOT esposito AT gmail DOT com
##
## Copyright (c) 2017, Emilio Xavier Esposito
##
## Permission is hereby granted, free of charge, to any person obtaining
## a copy of this software and associated documentation files (the
## "Software"), to deal in the Software without restriction, including
## without limitation the rights to use, copy, modify, merge, publish,
## distribute, sublicense, and/or sell copies of the Software, and to
## permit persons to whom the Software is furnished to do so, subject to
## the following conditions:
##
## The above copyright notice and this permission notice shall be
## included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
## EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
## MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
## LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
## OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
## WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##
## Based on http://opensource.org/licenses/MIT
##


## SUMMARY BARPLOTS (log10) ----------------------------------------------------

## Occupancy Summary Barplots (log10) docs -------------------------------------
#' @title Occupancy Summary Barplots
#' @description Occupancy summary barplots for the PDB structures. The plots are
#'   faceted and displays the binned occupancy values for all the structures.
#'   The counts are presented on a `log10` scale.
#'
#' @param data The results from the [CleanProteinStructures()] function. Will
#'   use the binned occupancy data.
#'
#' @export
#'
#' @import ggplot2
#' @import reshape2
#'
#' @examples
#'   \dontrun{
#'   OccupancyBarplot.summ(data)
#'
#'   ##----- multiple pages
#'   library(ggforce)
#'   occ.barplots.summary <- OccupancyBarplot.summ(data)
#'   num.pages <- ceiling(nrow(data$occupancy.counts) / 10)
#'
#'   pdf(file="multiple_pages.pdf", height=11, width=8.5)
#'   for (page in seq_len(num.pages)) {
#'     print(occ.barplots.summary +
#'           ggforce::facet_wrap_paginate(~PDBid,
#'                                        ncol = 2, nrow = 5, page = page) )
#'   }
#'   dev.off()
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family plots
#'
OccupancyBarplot.summ <- function(data) {

  ##_ extract occupancy.counts from cleaning results -----
  occupancy.counts <- data$occupancy.counts

  ##_ melt the data -----
  occ.summ.melt <- reshape2::melt(occupancy.counts, id.var = "PDBid")

  ##_ construct the plot -----
  occ.summ.plots <- ggplot2::ggplot(data = occ.summ.melt,
                                    aes_string(x = "variable", y = "value")) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_wrap(~PDBid, ncol = 2) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(title = "Occupancy Summary (Initial Structures)",
                  x = "Occupancy Values",
                  y = "Counts (log10)") +
    ggplot2::scale_x_discrete(labels = seq(from = 0, to = 1, by = 0.1)) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = element_text(angle = 90,
                                              size = 10,
                                              hjust = 1,
                                              vjust = 0.5))

  ##_ return plot -----
  occ.summ.plots

}


## B-value Summary Barplots (log10) docs -------------------------------------
#' @title B-value Summary Barplots
#' @description B-value summary barplots for the PDB structures. The plots are
#'   faceted and displays the binned B-value values for all the structures.
#'   The counts are presented on a `log10` scale.
#'
#' @param data The results from the [CleanProteinStructures()] function. Will use
#'   the binned B-value data.
#'
#' @export
#'
#' @import ggplot2
#' @import reshape2
#'
#' @examples
#'   \dontrun{
#'   BvalueBarplot.summ(data)
#'
#'   ##----- multiple pages
#'   library(ggforce)
#'   Bvalue.barplots.summary <- BvalueBarplot.summ(data)
#'   num.pages <- ceiling(nrow(data$Bvalue.counts) / 10)
#'
#'   pdf(file="multiple_pages.pdf", height=11, width=8.5)
#'   for (page in seq_len(num.pages)) {
#'     print(Bvalue.barplots.summary +
#'           ggforce::facet_wrap_paginate(~PDBid,
#'                                        ncol = 2, nrow = 5, page = page) )
#'   }
#'   dev.off()
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family plots
#'
BvalueBarplot.summ <- function(data) {

  ##_ extract occupancy.counts from cleaning results -----
  Bvalue.counts <- data$Bvalue.counts

  ##_ melt the data -----
  bvalue.summ.melt <- reshape2::melt(Bvalue.counts, id.var = "PDBid")

  ##_ construct the plot -----
  bvalue.summ.plots <- ggplot2::ggplot(data=bvalue.summ.melt,
                                       aes_string(x = "variable", y = "value")) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_wrap(~PDBid, ncol = 2) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(title = "B-value Summary (Initial Structures)",
                  x = "B-values",
                  y = "Counts (log10)") +
    ggplot2::scale_x_discrete(labels = seq(from = 0, to = 100, by = 5)) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = element_text(angle = 90, size = 10,
                                              hjust = 1, vjust = 0.5))

  ##_ return plot -----
  bvalue.summ.plots

}


## Normalized B-value Summary Barplots (log10) docs ---------------------------
#' @title Normalized B-value Summary Barplots
#' @description B-value summary barplots for the PDB structures. The plots are
#'   faceted and displays the binned B-value values for all the structures.
#'   The counts are presented on a `log10` scale.
#'
#' @param data The results from the [CleanProteinStructures()] function. Will
#'   use the binned normalized B-value data.
#'
#' @export
#'
#' @import ggplot2
#' @import reshape2
#'
#' @examples
#'   \dontrun{
#'   normBvalueBarplot.summ(data)
#'
#'   ##----- multiple pages
#'   library(ggforce)
#'   nBvalue.barplots.summary <- normBvalueBarplot.summ(data)
#'   num.pages <- ceiling(nrow(data$normBvalue.counts) / 10)
#'
#'   pdf(file="multiple_pages.pdf", height=11, width=8.5)
#'   for (page in seq_len(num.pages)) {
#'     print(nBvalue.barplots.summary +
#'           ggforce::facet_wrap_paginate(~PDBid,
#'                                        ncol = 2, nrow = 5, page = page) )
#'   }
#'   dev.off()
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family plots
#'
normBvalueBarplot.summ <- function(data) {

  ##----- extract occupancy.counts from cleaning results
  normBvalue.counts <- data$normBvalue.counts

  ##----- sum the less populated z-scores
  ##--- -7 to -2.5 z-scores
  zScores.lower.sums <- rowSums(normBvalue.counts[, 2:47])
  ##--- 3 to 7 z-scores
  zScores.upper.sums <- rowSums(normBvalue.counts[, 102:142])
  ##--- construct a new data.frame
  nBvalue.counts.short <- data.frame(PDBid=normBvalue.counts[, 1],
                                     zScores.lower.sums,
                                     normBvalue.counts[, 48:101],
                                     zScores.upper.sums,
                                     stringsAsFactors = FALSE)
  ##--- construct x-axis tick labels
  x.axis.tick.labels <- round(seq(from = -2.5, to = 3, by = 0.1), digits = 1)
  x.axis.tick.minor.tf <- rep(c(FALSE, TRUE, TRUE, TRUE, TRUE), 11)
  x.axis.tick.labels[x.axis.tick.minor.tf] <- ""
  x.axis.tick.labels[1] <- expression(paste('\u2264', "-2.5", sep = ""))
  x.axis.tick.labels[length(x.axis.tick.labels)] <-
    expression(paste('\u2265', "3", sep = ""))

  ##----- melt the data
  nBvalue.counts.melt <- reshape2::melt(nBvalue.counts.short, id.var = "PDBid")

  ##----- construct the plot
  nBvalue.counts.plots <- ggplot2::ggplot(data = nBvalue.counts.melt,
                                          aes_string(x = "variable", y = "value")) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_wrap(~PDBid, ncol = 2) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(title = "Normalized B-value Summary (Initial Structures)",
                  x = "Z-Score",
                  y = "Counts (log10)") +
    ggplot2::scale_x_discrete(labels = x.axis.tick.labels) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = element_text(size = 10))

  ##----- return plot
  nBvalue.counts.plots

}


## Mobility Summary Barplots (log10) docs --------------------------------------
#' @title Mobility Summary Barplots
#' @description Mobility summary barplots for PDB structures. The plots are
#'   faceted and displays the binned B-value values for all the structures. The
#'   counts are presented on a `log10` scale. The function will automatically
#'   plot ten plots per page.
#'
#' @param data The results from the [CleanProteinStructures()] function. Will
#'   use the binned mobility data.
#'
#' @export
#'
#' @import ggplot2
#' @import reshape2
#'
#' @examples
#'   \dontrun{
#'   MobilityBarplot.summ(data)
#'
#'   ##----- multiple pages
#'   library(ggforce)
#'   mob.barplots.summary <- MobilityBarplot.summ(data)
#'   num.pages <- ceiling(nrow(data$mobility.counts) / 10)
#'
#'   pdf(file="multiple_pages.pdf", height=11, width=8.5)
#'   for (page in seq_len(num.pages)) {
#'     print(mob.barplots.summary +
#'           ggforce::facet_wrap_paginate(~PDBid,
#'                                        ncol = 2, nrow = 5, page = page) )
#'   }
#'   dev.off()
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family plots
#'
MobilityBarplot.summ <- function(data) {

  ##----- extract occupancy.counts from cleaning results
  mobility.counts <- data$mobility.counts

  ##----- sum the less populated z-scores
  ##--- 3 to 6 z-scores
  zScores.upper.sums <- rowSums(mobility.counts[, 32:62])
  ##--- construct a new data.frame
  mobility.counts.short <- data.frame(PDBid=mobility.counts[, 1],
                                      mobility.counts[, 2:31],
                                      zScores.upper.sums,
                                      stringsAsFactors = FALSE)
  ##--- construct x-axis tick labels
  x.axis.tick.labels <- round(seq(from = 0, to = 3, by = 0.1), digits = 1)
  x.axis.tick.minor.tf <- rep(c(FALSE, TRUE, TRUE, TRUE, TRUE), 6)
  x.axis.tick.labels[x.axis.tick.minor.tf] <- ""
  x.axis.tick.labels[length(x.axis.tick.labels)] <-
    expression(paste('\u2265', "3", sep = ""))

  ##----- melt the data
  mobility.counts.melt <- reshape2::melt(mobility.counts.short, id.var = "PDBid")

  ##----- construct the plot
  mobility.counts.plots <- ggplot2::ggplot(data = mobility.counts.melt,
                                           aes_string(x = "variable",
                                                      y = "value")) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_wrap(~PDBid, ncol = 2) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(title = "Mobility Summary (Initial Structures)",
                  x = "Z-Score",
                  y = "Counts (log10)") +
    ggplot2::scale_x_discrete(labels = x.axis.tick.labels) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = element_text(size = 10
                                              # , angle = 90,
                                              # hjust = 1,
                                              # vjust = 0.5
                   ) )

  ##----- return plot
  mobility.counts.plots

}


## RESULTS FOCUSED PLOTS -------------------------------------------------------

## Conservation Plot docs ------------------------------------------------------
#' @title Conservation Plot (Number of Waters Per Cluster Histogram)
#' @description Histogram and density plots for number of cluster with number of
#'   atoms
#' @details Constructs a histogram for the number of waters per cluster.
#'   Clusters with less than 50\% conservation are light grey, clusters with 50
#'   to 69\% water conseration are dark red, clusters with 70 to 79\%
#'   conservation are red, 80 to 89\% conservation are light blue, 90 to 99\%
#'   conservation are blue, and 100\% conservation (waters from all structures)
#'   are dark blue.
#'
#' @param data The `h2o.clusters.summary` data.frame from the
#'   [ClusterWaters()] function containing the `num.waters` information.
#'   The `num.waters` values are integers.
#' @param passed.waters Logical indicator to plot results for waters **passing**
#'   [Mobility()] and [NormalizedBvalue()] _**OR**_ using **all** waters within
#'   the `PDB` files.
#'
#' @export
#'
#' @import ggplot2
#'
#' @examples
#'   \dontrun{
#'   Conservation.plot <- ConservationPlot(data=thrombin10.conservedWaters,
#'                                         passed.waters=TRUE)
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Results Visualization" "Results Plots"
#'
#' @references
#'   Paul C Sanschagrin and Leslie A Kuhn. Cluster analysis of
#'   consensus water sites in thrombin and trypsin shows conservation between
#'   serine proteases and contributions to ligand specificity. _Protein
#'   Science_, 1998, **7** (_10_), pp 2054-2064.
#'   [DOI: 10.1002/pro.5560071002](http://doi.org/10.1002/pro.5560071002)
#'   [PMID: 9792092](http://www.ncbi.nlm.nih.gov/pubmed/9792092)
#'   [WatCH webpage](http://www.kuhnlab.bmb.msu.edu/software/watch/index.html)
#'
ConservationPlot <- function(data, passed.waters=TRUE) {

  ##_ get the PASSED or ALL waters plot data -----
  if (passed.waters == TRUE) {
    plot.data <- data$h2o.cluster.passed$h2o.clusters.summary
  } else {
    plot.data <- data$h2o.cluster.all$h2o.clusters.summary
  }

  ##_ add the conservation sets if not present -----
  if ( !any(colnames(plot.data) == "conserve.set") ) {
    plot.data$conserve.set <- ConservationSet(plot.data$pct.conserved)
  }

  ##_ determine binwidth based on the largest number of waters in a cluster -----
  num.waters.max <- max(plot.data$num.waters)
  if ( num.waters.max >= 100 ) bin.width <- 5
  if ( num.waters.max <  100 ) bin.width <- 2
  if ( num.waters.max <=  25 ) bin.width <- 1

  ##_ create the plot -----
  hist.plot <- ggplot2::ggplot(data = plot.data,
                               aes_string(x = "num.waters",
                                          fill = "conserve.set")) +
    ggplot2::geom_histogram(stat = "bin", binwidth = bin.width) +
    ##__ set the fill and line colors -----
    ggplot2::scale_fill_manual(values = cons.color6,
                               name = cons.color6.legend,
                               breaks = cons.color6.breaks,
                               labels = cons.color6.labels) +
    ##__ format axis labels -----
    # hist.plot <- hist.plot + scale_y_continuous(breaks=NULL)
    ggplot2::labs(title = "Number of Waters per Cluster",
                  x="Number of Protein Structures",
                  y="Number of Water Clusters") +
    ggplot2::theme(legend.position = "bottom")

  ##_ return plot -----
  hist.plot

}


## Occupancy Barplot docs ------------------------------------------------------
#' @title Occupancy Barplots
#' @description Occupancy Barplots for Cluster with at least 50\% Conservation
#' @details Constructs a barplot with corresponding density plot for the mean
#'   occupancy value for all water within each cluster with at least 50\% water
#'   conservation. Clusters with 50 to 69\% water conseration are dark red,
#'   clusters with 70 to 79\% conservation are red, 80 to 89\% conservation are
#'   light blue, 90 to 99\% conservation are blue, and 100\% conservation (waters
#'   from all structures) are dark blue.
#'
#'   This plot was inspired by Figure 1 of Sanschagrin and Kuhn (1998).
#'
#' @param data The `h2o.clusters.summary` data.frame from the
#'   [ClusterWaters()] function containing the `o.mu` information.
#' @param passed.waters Logical indicator to plot results for waters **passing**
#'   [Mobility()] and [NormalizedBvalue()] _**OR**_ using **all** waters within
#'   the `PDB` files.
#'
#' @export
#'
#' @import ggplot2
#'
#' @examples
#'   \dontrun{
#'   occupancy.plot <- OccupancyBarplot(data=thrombin10.conservedWaters,
#'                                      passed.waters=TRUE)
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family plots
#'
#' @references
#'   Paul C Sanschagrin and Leslie A Kuhn. Cluster analysis of
#'   consensus water sites in thrombin and trypsin shows conservation between
#'   serine proteases and contributions to ligand specificity. _Protein
#'   Science_, 1998, **7** (_10_), pp 2054-2064.
#'   [DOI: 10.1002/pro.5560071002](http://doi.org/10.1002/pro.5560071002)
#'   [PMID: 9792092](http://www.ncbi.nlm.nih.gov/pubmed/9792092)
#'   [WatCH webpage](http://www.kuhnlab.bmb.msu.edu/software/watch/index.html)
#'
OccupancyBarplot <- function(data, passed.waters=TRUE) {

  ##_ get the PASSED or ALL waters plot data -----
  if (passed.waters == TRUE) {
    plot.data <- data$h2o.cluster.passed$h2o.clusters.summary
  } else {
    plot.data <- data$h2o.cluster.all$h2o.clusters.summary
  }

  ##_ only the waters with 50% or greater conservation ------
  plot.data <- plot.data[plot.data$pct.conserved >= 50, ]

  ##_ determine the minimum and maximum x-axis values -----
  x.axis.min <- floor(min(plot.data$o.mu))
  x.axis.max <- ceiling(max(plot.data$o.mu))

  ##_ add the conservation sets if not present -----
  if ( !any(colnames(plot.data) == "conserve.set") ) {
    plot.data$conserve.set <- ConservationSet(plot.data$pct.conserved)
  }

  ##_ construct the plot -----
  occ.plot <- ggplot2::ggplot(data = plot.data,
                              aes_string(x = "o.mu",
                                         fill = "conserve.set",
                                         colour = "conserve.set")) +
    ggplot2::geom_density(alpha = 0.2, size = 1) +
    ggplot2::geom_histogram(aes_string(y = "..density.."),
                            position = "dodge", stat = "bin",
                            binwidth = 0.05, colour = "white") +
    ##__ set the fill and line colors -----
    ggplot2::scale_fill_manual(values = cons.color5,
                               name = cons.color5.legend,
                               breaks = cons.color5.breaks,
                               labels = cons.color5.labels) +
    ggplot2::scale_colour_manual(values=cons.color5,
                                 name = cons.color5.legend,
                                 breaks = cons.color5.breaks,
                                 labels = cons.color5.labels) +
    ##__ format axis labels -----
    # occ.plot <- occ.plot + scale_y_continuous(breaks=NULL)
    ggplot2::labs(title = "Occupancy",
                  x = "Mean Occupancy Value", y = "") +
    ggplot2::theme(legend.position = "bottom")

  ##_ return plot -----
  occ.plot

}


## Mobility Barplot docs --------------------------------------------------
#' @title Mobility Barplots
#' @description Mobility Barplots for Cluster with at least 50\% Conservation
#' @details Constructs a barplot with corresponding density plot for the mean
#'   mobility value for all water within each cluster with at least 50\% water
#'   conservation. Clusters with 50 to 69\% water conseration are dark red,
#'   clusters with 70 to 79\% conservation are red, 80 to 89\% conservation are
#'   light blue, 90 to 99\% conservation are blue, and 100\% conservation (waters
#'   from all structures) are dark blue.
#'
#'   The mobility values are calculated by the [Mobility()] function.
#'
#'   This plot was inspired by Figure 1 of Sanschagrin and Kuhn (1998).
#'
#' @param data The `h2o.clusters.summary` data.frame from the
#'   [ClusterWaters()] function containing the `mobility.mu`
#'   information.
#' @param passed.waters Logical indicator to plot results for waters **passing**
#'   [Mobility()] and [NormalizedBvalue()] _**OR**_ using **all** waters within
#'   the `PDB` files.
#'
#' @export
#'
#' @import ggplot2
#'
#' @examples
#'   \dontrun{
#'   mobility.plot <- MobilityBarplot(data=thrombin10.conservedWaters,
#'                                    passed.waters=TRUE)
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family plots
#'
#' @references
#'   Paul C Sanschagrin and Leslie A Kuhn. Cluster analysis of
#'   consensus water sites in thrombin and trypsin shows conservation between
#'   serine proteases and contributions to ligand specificity. _Protein
#'   Science_, 1998, **7** (_10_), pp 2054-2064.
#'   [DOI: 10.1002/pro.5560071002](http://doi.org/10.1002/pro.5560071002)
#'   [PMID: 9792092](http://www.ncbi.nlm.nih.gov/pubmed/9792092)
#'   [WatCH webpage](http://www.kuhnlab.bmb.msu.edu/software/watch/index.html)
#'
MobilityBarplot <- function(data, passed.waters=TRUE) {

  ##_ get the PASSED or ALL waters plot data -----
  if (passed.waters == TRUE) {
    plot.data <- data$h2o.cluster.passed$h2o.clusters.summary
  } else {
    plot.data <- data$h2o.cluster.all$h2o.clusters.summary
  }

  ##_ only the waters with 50% or greater conservation ------
  plot.data <- plot.data[plot.data$pct.conserved >= 50, ]

  ##_ determine the minimum and maximum x-axis values -----
  x.axis.min <- floor(min(plot.data$mobility.mu))
  x.axis.max <- ceiling(max(plot.data$mobility.mu))

  ##_ add the conservation sets if not present -----
  if ( !any(colnames(plot.data) == "conserve.set") ) {
    plot.data$conserve.set <- ConservationSet(plot.data$pct.conserved)
  }

  ##_ construct the plot -----
  mobility.plot <- ggplot2::ggplot(data = plot.data,
                                   aes_string(x = "mobility.mu",
                                              fill = "conserve.set",
                                              colour = "conserve.set")) +
    ggplot2::geom_density(alpha = 0.4, size = 1) +
    ggplot2::geom_histogram(aes_string(y = "..density.."),
                            position = "dodge", stat = "bin",
                            binwidth = 0.10, colour = "white") +
    ##__ set the fill and line colors -----
    ggplot2::scale_fill_manual(values = cons.color5,
                               name = cons.color5.legend,
                               breaks = cons.color5.breaks,
                               labels = cons.color5.labels) +
    ggplot2::scale_colour_manual(values = cons.color5,
                                 name = cons.color5.legend,
                                 breaks = cons.color5.breaks,
                                 labels = cons.color5.labels) +
    ##__ format axis labels ----- -----
    # mobility.plot <- mobility.plot + scale_y_continuous(breaks=NULL)
    ggplot2::scale_x_continuous(breaks = seq(x.axis.min, x.axis.max, 0.5)) +
    ggplot2::labs(title = "Mobility",
                  x = "Mean Mobility Value", y = "") +
    ggplot2::theme(legend.position = "bottom")

  ##_ return the plot -----
  mobility.plot

}


## B-value Barplot docs -------------------------------------------------------
#' @title B-value Barplots
#' @description B-value Barplots for Cluster with at least 50\% Conservation
#' @details Constructs a barplot with corresponding density plot for the mean
#'   B-value value for all water within each cluster with at least 50\% water
#'   conservation. Clusters with 50 to 69\% water conseration are dark red,
#'   clusters with 70 to 79\% conservation are red, 80 to 89\% conservation are
#'   light blue, 90 to 99\% conservation are blue, and 100\% conservation (waters
#'   from all structures) are dark blue.
#'
#'   This plot was inspired by Figure 1 of Sanschagrin and Kuhn (1998).
#'
#' @param data The `h2o.clusters.summary` data.frame from the
#'   [ClusterWaters()] function containing the `b.exp.mu`
#'   information.
#' @param passed.waters Logical indicator to plot results for waters **passing**
#'   [Mobility()] and [NormalizedBvalue()] _**OR**_ using **all** waters within
#'   the `PDB` files.
#' @param calc.values Plot the calculated B-values the mean experimental B-values; default: TRUE
#'
#' @export
#'
#' @import ggplot2
#'
#' @examples
#'   \dontrun{
#'   Bvalue.plot <- BvalueBarplot(data=thrombin10.conservedWaters,
#'                                passed.waters=TRUE)
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family plots
#'
#' @references
#'   Paul C Sanschagrin and Leslie A Kuhn. Cluster analysis of
#'   consensus water sites in thrombin and trypsin shows conservation between
#'   serine proteases and contributions to ligand specificity. _Protein
#'   Science_, 1998, **7** (_10_), pp 2054-2064.
#'   [DOI: 10.1002/pro.5560071002](http://doi.org/10.1002/pro.5560071002)
#'   [PMID: 9792092](http://www.ncbi.nlm.nih.gov/pubmed/9792092)
#'   [WatCH webpage](http://www.kuhnlab.bmb.msu.edu/software/watch/index.html)
#'
BvalueBarplot <- function(data, passed.waters=TRUE, calc.values=TRUE) {

  ##_ get the PASSED or ALL waters plot data -----
  if (passed.waters == TRUE) {
    plot.data <- data$h2o.cluster.passed$h2o.clusters.summary
  } else {
    plot.data <- data$h2o.cluster.all$h2o.clusters.summary
  }

  ##_ experimental mean or calculated B-values -----
  if (calc.values == TRUE) {
    x.column.name <- "b.calc"
    x.axis.name <- "Calculated B-value"
    title <- "B-value (Calculated from cluster's RMSF)"
  } else {
    x.column.name <- "b.exp.mu"
    x.axis.name <- "Mean B-value"
    title <- "B-value (Mean calculated from cluster's experimental B-values)"
  }

  ##_ only the waters with 50% or greater conservation ------
  plot.data <- plot.data[plot.data$pct.conserved >= 50, ]

  ##----- determine the minimum and maximum x-axis values
  # x.axis.min <- floor(min(plot.data[[x.column.name]]))
  # x.axis.max <- ceiling(max(plot.data[[x.column.name]]))

  ##_ add the conservation sets if not present -----
  if ( !any(colnames(plot.data) == "conserve.set") ) {
    plot.data$conserve.set <- ConservationSet(plot.data$pct.conserved)
  }

  ##----- construct the plot
  Bvalue.plot <- ggplot2::ggplot(data = plot.data,
                                 aes_string(x = x.column.name,
                                            fill = "conserve.set",
                                            colour = "conserve.set")) +
    ggplot2::geom_density(alpha = 0.2, size = 1) +
    ggplot2::geom_histogram(aes_string(y = "..density.."),
                            position = "dodge", stat = "bin",
                            binwidth = 5, colour = "white") +
    ##--- set the fill and line colors
    ggplot2::scale_fill_manual(values = cons.color5,
                               name = cons.color5.legend,
                               breaks = cons.color5.breaks,
                               labels = cons.color5.labels) +
    ggplot2::scale_colour_manual(values = cons.color5,
                                 name = cons.color5.legend,
                                 breaks = cons.color5.breaks,
                                 labels = cons.color5.labels) +
    ##--- format axis labels
    # scale_y_continuous(breaks = NULL) +
    # scale_x_continuous(breaks = seq(x.axis.min, x.axis.max, 0.250)) +
    ggplot2::labs(title = title,
                  x = x.axis.name, y = "") +
    ggplot2::theme(legend.position = "bottom")

  ##----- return the plot
  Bvalue.plot

}


## Normalized B-value Barplot docs --------------------------------------------
#' @title Normalized B-value Barplots
#' @description Normalized B-value Barplots for Cluster with at least 50\%
#'   Conservation
#' @details Constructs a barplot with corresponding density plot for the mean
#'   normalized B-value value for all water within each cluster with at least
#'   50\% water conservation. Clusters with 50 to 69\% water conseration are dark
#'   red, clusters with 70 to 79\% conservation are red, 80 to 89\% conservation
#'   are light blue, 90 to 99\% conservation are blue, and 100\% conservation
#'   (waters from all structures) are dark blue.
#'
#'   The normalized B-value values are calculated by the
#'   [NormalizedBvalue()] function.
#'
#'   This plot was inspired by Figure 1 of Sanschagrin and Kuhn (1998).
#'
#' @param data The `h2o.clusters.summary` data.frame from the
#'   [ClusterWaters()] function containing the `nBvalue.mu`
#'   information.
#' @param passed.waters Logical indicator to plot results for waters **passing**
#'   [Mobility()] and [NormalizedBvalue()] _**OR**_ using **all** waters within
#'   the `PDB` files.
#'
#' @export
#'
#' @import ggplot2
#'
#' @examples
#'   \dontrun{
#'   nBvalue.plot <- nBvalueBarplot(data=thrombin10.conservedWaters,
#'                                  passed.waters=TRUE)
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family plots
#'
#' @references
#'   Paul C Sanschagrin and Leslie A Kuhn. Cluster analysis of
#'   consensus water sites in thrombin and trypsin shows conservation between
#'   serine proteases and contributions to ligand specificity. _Protein
#'   Science_, 1998, **7** (_10_), pp 2054-2064.
#'   [DOI: 10.1002/pro.5560071002](http://doi.org/10.1002/pro.5560071002)
#'   [PMID: 9792092](http://www.ncbi.nlm.nih.gov/pubmed/9792092)
#'   [WatCH webpage](http://www.kuhnlab.bmb.msu.edu/software/watch/index.html)
#'
nBvalueBarplot <- function(data, passed.waters=TRUE) {

  ##_ get the PASSED or ALL waters plot data -----
  if (passed.waters == TRUE) {
    plot.data <- data$h2o.cluster.passed$h2o.clusters.summary
  } else {
    plot.data <- data$h2o.cluster.all$h2o.clusters.summary
  }

  ##_ only the waters with 50% or greater conservation ------
  plot.data <- plot.data[plot.data$pct.conserved >= 50, ]

  ##_ determine the minimum and maximum x-axis values -----
  x.axis.min <- floor(min(plot.data$nBvalue.mu))
  x.axis.max <- ceiling(max(plot.data$nBvalue.mu))

  ##_ add the conservation sets if not present -----
  if ( !any(colnames(plot.data) == "conserve.set") ) {
    plot.data$conserve.set <- ConservationSet(plot.data$pct.conserved)
  }

  ##_ construct the plot -----
  nBvalue.plot <- ggplot2::ggplot(data = plot.data,
                                  aes_string(x = "nBvalue.mu",
                                             fill = "conserve.set",
                                             colour = "conserve.set")) +
    ggplot2::geom_density(alpha = 0.2, size =1) +
    ggplot2::geom_histogram(aes_string(y = "..density.."),
                            position = "dodge", stat = "bin",
                            binwidth = 0.25, colour = "white") +
    ##--- set the fill and line colors
    ggplot2::scale_fill_manual(values = cons.color5,
                               name = cons.color5.legend,
                               breaks = cons.color5.breaks,
                               labels = cons.color5.labels) +
    ggplot2::scale_colour_manual(values = cons.color5,
                                 name = cons.color5.legend,
                                 breaks = cons.color5.breaks,
                                 labels = cons.color5.labels) +
    ##__ format axis labels -----
    # nBvalue.plot <- nBvalue.plot + scale_y_continuous(breaks=NULL)
    ggplot2::scale_x_continuous(breaks = seq(x.axis.min, x.axis.max, 0.5)) +
    ggplot2::labs(title = "Normalized B-value",
                  x = "Mean Normalized B-value", y = "") +
    ggplot2::theme(legend.position = "bottom")

  ##_ return the plot -----
  nBvalue.plot

}


## Cluster Summary Plots docs --------------------------------------------------
#' @title Cluster Summary Plots
#' @description Collection of cluster summary plots.
#' @details The Number of Water Cluster (see [ConservationPlot()]), Occupancy
#'   (see [OccupancyBarplot()]), Mobility (see [MobilityBarplot()]), B-value
#'   (see [BvalueBarplot()]), and Normalized B-value (see [nBvalueBarplot()])
#'   plots are combined into a single plot image. The ability to label each plot
#'   with capital letters (upper-case) or lower-case letters is available.
#'
#' @param data The results from the [ConservedWaters()] function.
#' @param passed.waters Logical indicator to plot results for waters **passing**
#'   [Mobility()] and [NormalizedBvalue()] _**OR**_ using **all** waters within
#'   the `PDB` files.
#' @param plot.labels Using the same options as [cowplot::plot_grid()] plus
#'   `NULL`. The option `"AUTO"` labels each plot with upper-case letters
#'   (_e.g._, A, B, C, D, E), `"auto"` labels each plot with lower-case letters
#'   (_e.g._, a, b, c, d, e), and `NULL` returns plots without labels. Default
#'   is `NULL`.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom cowplot plot_grid
#'
#' @examples
#'   \dontrun{
#'   cluster.summary.plot <- ClusterSummaryPlots(data=thrombin10.conservedWaters,
#'                                              passed.waters=TRUE,
#'                                              labels=NULL)
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family plots
#'
#' @references
#'   Paul C Sanschagrin and Leslie A Kuhn. Cluster analysis of
#'   consensus water sites in thrombin and trypsin shows conservation between
#'   serine proteases and contributions to ligand specificity. _Protein
#'   Science_, 1998, **7** (_10_), pp 2054-2064.
#'   [DOI: 10.1002/pro.5560071002](http://doi.org/10.1002/pro.5560071002)
#'   [PMID: 9792092](http://www.ncbi.nlm.nih.gov/pubmed/9792092)
#'   [WatCH webpage](http://www.kuhnlab.bmb.msu.edu/software/watch/index.html)
#'
ClusterSummaryPlots <- function(data, passed.waters=TRUE, plot.labels=NULL) {

  ##_ check plot.labels parameter -----
  plot.labels.available <- c("AUTO", "auto")
  if ( !any(plot.labels.available %in% plot.labels) & !is.null(plot.labels) ) {
    mess <- paste("The provided plot.labels value of", plot.labels,
                  "is not valid and has been replace with NULL. Please",
                  "consider using \"AUTO\", \"auto\", or NULL.", sep = " ")
    message(mess)
    plot.labels <- NULL
  }

  ##_ construct the individual plots -----
  h2o.per.cluster.hist <- ConservationPlot(data=data, passed.waters=passed.waters)
  occupancy.plot <- OccupancyBarplot(data=data, passed.waters=passed.waters)
  mobility.plot <- MobilityBarplot(data=data, passed.waters=passed.waters)
  Bvalue.plot <- BvalueBarplot(data=data, passed.waters=passed.waters)
  nBvalue.plot <- nBvalueBarplot(data=data, passed.waters=passed.waters)

  ##_ combine the plots -----
  cluster.summary.plots <- cowplot::plot_grid(
    h2o.per.cluster.hist + ggplot2::theme(legend.position="none"),
    occupancy.plot + ggplot2::theme(legend.position="none"),
    mobility.plot + ggplot2::theme(legend.position="none"),
    Bvalue.plot + ggplot2::theme(legend.position="none"),
    nBvalue.plot + ggplot2::theme(legend.position="none"),
    align = 'v',
    labels = plot.labels,
    hjust = -1,
    nrow = 5)

  ##_ legend at the bottom of all plots -----
  ##  extract the legend from the water per cluster histogram
  grobs <- ggplot2::ggplotGrob(h2o.per.cluster.hist +
                                 # guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
                                 ggplot2::theme(legend.position = "bottom",
                                                legend.box = "horizontal"))$grobs
  legend.bottom <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
  ##__ add the legend under the last plot in the column. the legend is given a -----
  ##   of 10% of the height of a plot row using the rel_heights option
  cluster.summary.plots <- cowplot::plot_grid(cluster.summary.plots,
                                              legend.bottom,
                                              ncol=1,
                                              rel_heights=c(1, .1))

  ##_ return the plot -----
  cluster.summary.plots

}


## Mobility and Normalized B-values evaluation plots -----------------------------
#' @title Mobility and Normalized B-values Evaluation Plots
#' @description Mean bound water environment summary per percent conservation
#' @details Constructs a series of scatterplots illustrating the relationship
#'   between mobility and normalized B-values and (i) percent water
#'   conservation, (ii) mean distance between waters in a cluster, and (iii)
#'   mean distance between waters in a cluster and the cluster's centroid. The
#'   dots are colored based on water cluster "percent conservation":
#'   - **light grey dots**: less than 50\% conservation
#'   - **dark red dots**: 50\% to 69\% conservations
#'   - **red dots**: 70\% to 79\% conservation
#'   - **light blue dots**: 80\% to 89\% conservation
#'   - **blue dots**: 90\% to 99\% conservation
#'   - **dark blue dots**: 100\% conservation (all structures contribute to the
#'     water cluster).
#'
#'   The mean distance plots will have a column of dots at a distance of 0.0 if
#'   there are clusters composed of a single water molecule. Thus, these
#'   clusters have a zero distance between and to other waters in their cluster
#'   because there are _**no other waters**_ in their cluster.
#'
#'   This plot was inspired by Figure 2 of Ogata and Wodak (2002).
#'
#' @param data The `h2o.clusters.summary` data.frame from the
#'   [ConservedWaters()] function containing the `nBvalue.mu`
#'   information. This data.frame is found within the `h2o.cluster.passed`
#'   and `h2o.cluster.all`
#' @param passed.waters Logical indicator to plot results for waters **passing**
#'   [Mobility()] and [NormalizedBvalue()] _**OR**_ using **all** waters within
#'   the `PDB` files.
#' @param title The title for the plot
#'
#' @export
#'
#' @import ggplot2
#' @import reshape2
#' @import scales
#' @importFrom cowplot ggdraw draw_label plot_grid
#'
#' @examples
#'   \dontrun{
#'   bwe.summary.plot <- MobNormBvalEvalPlots(data=thrombin10.conservedWaters,
#'                                            passed.waters=TRUE,
#'                                            title="Mobility and Normalized B-value Evaluation")
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family plots
#'
#' @references
#'   Koji Ogata and Shoshana J Wodak. Conserved water molecules in MHC
#'   class-I molecules and their putative structural and functional roles.
#'   _Protein Engineering_, 2002, **15** (_8_), pp 697-705.
#'   [DOI: 10.1093/protein/15.8.697](http://doi.org/10.1093/protein/15.8.697)
#'   [PMID: 12364585](http://www.ncbi.nlm.nih.gov/pubmed/12364585)
#'
MobNormBvalEvalPlots <- function(data,
                                 passed.waters = TRUE,
                                 title = "Mobility and Normalized B-value Evaluation") {

  ##_ get the PASSED or ALL waters plot data -----
  if (passed.waters == TRUE) {
    plot.data <- data$h2o.cluster.passed$h2o.clusters.summary
  } else {
    plot.data <- data$h2o.cluster.all$h2o.clusters.summary
  }

  ##_ add the conservation sets if not present -----
  if ( !any(colnames(plot.data) == "conserve.set") ) {
    plot.data$conserve.set <- ConservationSet(plot.data$pct.conserved)
  }

  ##_ convert percent to fraction -----
  # plot.data$pct.conserved <- (plot.data$pct.conserved / 100)

  ##_ setup the names for the facet labeling -----
  conserve.set.names <- c("set0", "set1", "set2", "set3", "set4", "set5")
  facet.labels <- c("Mobility", "Normalized B-values")

  ##_ create the needed data.frames -----
  ##__ percent conservation -----
  df.conserv <- plot.data[, c("pct.conserved",
                              "mobility.mu",
                              "nBvalue.mu",
                              "conserve.set")]
  df.conserv$group <- "PercentConserved"
  colnames(df.conserv)[1] <- "X.values"
  ##__ mean distance between waters -----
  df.dist <- plot.data[, c("dist.mu",
                           "mobility.mu",
                           "nBvalue.mu",
                           "conserve.set")]
  df.dist$group <- "DistWaters"
  colnames(df.dist)[1] <- "X.values"
  ##__ mean distance between the water and the centroid -----
  df.cent <- plot.data[, c("dist.centroid.mu",
                           "mobility.mu",
                           "nBvalue.mu",
                           "conserve.set")]
  df.cent$group <- "DistCentroid"
  colnames(df.cent)[1] <- "X.values"

  ##__ combine all the data -----
  df.plot.data <- rbind(df.conserv,
                        df.dist,
                        df.cent)

  ##_ melt the plot.data -----
  df.plot.data.melt <- reshape2::melt(df.plot.data[c(nrow(df.plot.data):1), ],
                                      id.vars = c("X.values", "conserve.set", "group") )
  df.plot.data.melt$conserve.set <- factor(df.plot.data.melt$conserve.set,
                                           levels = conserve.set.names)
  df.plot.data.melt$variable <- factor(df.plot.data.melt$variable,
                                       labels = c(expression(paste("Mobility", sep="")),
                                                  expression(paste("Normalized B-values", sep=""))
                                       )
  )
  df.plot.data.melt$group <- factor(df.plot.data.melt$group,
                                    levels = c("PercentConserved",
                                               "DistWaters",
                                               "DistCentroid")
                                    )
  df.plot.data.melt$group <- factor(df.plot.data.melt$group,
                                       labels = c(expression(paste("Percent Conservation", sep="")),
                                                  expression(paste("Mean Distance Between Waters (", ring(A), ")", sep = "")),
                                                  expression(paste("Mean Distance to Centroid (", ring(A), ")", sep = ""))
                                                  )
                                    )

  ##_ make the plot -----
  summary.plots <- ggplot2::ggplot(df.plot.data.melt,
                                   aes_string(x = "X.values", y = "value")) +
    ggplot2::geom_point(aes_string(colour = "conserve.set"), size = 3) +
    ggplot2::facet_grid(variable ~ group, scales = "free", labeller = label_parsed) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::theme(legend.position = "bottom", legend.box = "horizontal") +
    ggplot2::scale_colour_manual(values = cons.color6,
                                 name = cons.color6.legend,
                                 breaks = cons.color6.breaks,
                                 labels = cons.color6.labels)

  ##_ add the title -----
  title.main <- cowplot::ggdraw() + cowplot::draw_label(title, fontface = 'bold')
  ## rel_heights values control title margins
  summary.plots <- cowplot::plot_grid(title.main, summary.plots,
                                      ncol = 1, rel_heights = c(0.04, 1))


  ##_ return the plot -----
  summary.plots

}


## bound water environment plots -----------------------------------------------
#' @title Bound Water Environment Barplots
#' @description Normalized B-value Barplots for Cluster with at least 50\%
#'   Conservation
#' @details Constructs a barplot with corresponding density plot for the mean
#'   normalized B-value value for all water within each cluster with at least
#'   50\% water conservation. Clusters with 50 to 69\% water conseration are dark
#'   red, clusters with 70 to 79\% conservation are red, 80 to 89\% conservation
#'   are light blue, 90 to 99\% conservation are blue, and 100\% conservation
#'   (waters from all structures) are dark blue.
#'
#'   The normalized B-value values are calculated by the
#'   \code{NormalizedBvalue} function.
#'
#'   This plot was inspired by Figure 1 of Sanschagrin and Kuhn (1998).
#'
#' @param data The \code{h2o.clusters.summary} data.frame from the
#'   \code{ClusterWaters} function containing the \code{nBvalue.mu}
#'   information.
#' @param passed.waters Logical indicator to plot results for waters **passing**
#'   [Mobility()] and [NormalizedBvalue()] _**OR**_ using **all** waters within
#'   the `PDB` files.
#' @param pct.conserved.gte minimum percent conservation within a water cluster;
#'   default: `50.0`If the number of clusters is less than the number of
#'   clusters defined by `num.clusters`, then the number of clusters defined by
#'   `pct.conserved.gte` is displayed.
#' @param num.clusters number (integer) of clusters to display. If the number of
#'   clusters is less than the number of clusters defined by
#'   `pct.conserved.gte`, then the number of clusters defined by `num.clusters`
#'   is displayed. A value of `NULL` results in the provided value for
#'   `pct.conserved.gte` being used.
#'
#' @export
#'
#' @import ggplot2
#' @import reshape2
#' @import scales
#'
#' @examples
#'   \dontrun{
#'   bwe.plots <- BoundWaterEnvPlots(data=thrombin10.conservedWaters,
#'                                   passed.waters=TRUE,
#'                                   pct.conserved.gte = 50.0,
#'                                   num.clusters = 50)
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family plots
#'
#' @references
#'   Paul C Sanschagrin and Leslie A Kuhn. Cluster analysis of
#'   consensus water sites in thrombin and trypsin shows conservation between
#'   serine proteases and contributions to ligand specificity. _Protein
#'   Science_, 1998, **7** (_10_), pp 2054-2064.
#'   [DOI: 10.1002/pro.5560071002](http://doi.org/10.1002/pro.5560071002)
#'   [PMID: 9792092](http://www.ncbi.nlm.nih.gov/pubmed/9792092)
#'   [WatCH webpage](http://www.kuhnlab.bmb.msu.edu/software/watch/index.html)
#'
BoundWaterEnvPlots <- function(data,
                               passed.waters = TRUE,
                               pct.conserved.gte = 50.0,
                               num.clusters = 50) {

  ##_ get the PASSED or ALL waters plot data -----
  if (passed.waters == TRUE) {
    plot.data <- data$h2o.cluster.passed$h2o.clusters.summary
  } else {
    plot.data <- data$h2o.cluster.all$h2o.clusters.summary
  }

  ##_ check user provided parameters -----
  ##__ percentage of conserved waters -----
  if ( pct.conserved.gte <= 1.0) {
    mess <- paste("Looks like you provided a fraction and not a percentage. ",
                  "The provided value (", pct.conserved.gte,
                  ") has been converted to ", pct.conserved.gte * 100, "%.",
                  sep = "")
    message(mess)
    pct.conserved.gte <- pct.conserved.gte * 100
  }

  ##__ number of clusters -----
  num.clusters <- as.integer(num.clusters)
  if ( num.clusters <= 1 ) {
    mess <- paste("The provided num.cluster value of", num.clusters,
                  "is less than or equal to 1 and has been changed to NULL.",
                  "The percent conversed value will be used instead.",
                  sep = " ")
    message(mess)
    num.clusters <- NULL
  }

  ##_ determine the clusters passing the pct.conserved.gte value -----
  pct.conserved.gte.tf <- plot.data$pct.conserved >= pct.conserved.gte
  num.pct.conserved.gte <- sum(pct.conserved.gte.tf)

  if ( !is.null(num.clusters) && (num.clusters <= num.pct.conserved.gte) ) {
    plot.data <- plot.data[1:num.clusters, ]
    mess <- paste("The number of requested conserved water clusters (",
                  num.clusters,
                  ") is less than the number (", num.pct.conserved.gte,
                  ") of conserved water clusters with ",
                  pct.conserved.gte, "% conservation. Thus, only ",
                  "the results for the top ",
                  num.clusters, " conserved water clusters will be plotted.",
                  sep = "")
    message(mess)
  }

  if ( is.null(num.clusters) ) {
    plot.data <- plot.data[pct.conserved.gte.tf, ]
    mess <- paste("The number of conserved water clusters with ",
                  pct.conserved.gte, "% conservation is ",
                  num.pct.conserved.gte, ".", sep = "")
    message(mess)
  }

  ##_ bound water environment plots -----
  cluster.name <- paste("cluster_", plot.data$cluster, sep="")
  num.clusters <- length(cluster.name)
  cluster.breaks <- seq(from = 10, to = num.clusters, by = 10)
  cluster.breaks <- c(1, cluster.breaks)
  cluster.numbering <- seq(from = 1, to = num.clusters)
  cluster.numbering <- rep("", num.clusters)
  cluster.numbering[cluster.breaks] <- cluster.breaks

  ##_ create the plot's data.frame -----
  df.prot <- data.frame(name = cluster.name,
                        adn = plot.data$prot.adn / plot.data$num.waters,
                        ahp = plot.data$prot.ahp.mu,
                        hbonds = plot.data$prot.hbonds / plot.data$num.waters,
                        mobiity = plot.data$prot.mobility.mu,
                        nBvalue = plot.data$prot.nBvalue.mu,
                        type = "Protein",
                        stringsAsFactors = FALSE)

  df.h2o <- data.frame(name = cluster.name,
                       adn = plot.data$h2o.adn / plot.data$num.waters,
                       ahp = NA,
                       hbonds = plot.data$h2o.hbonds / plot.data$num.waters,
                       mobiity = plot.data$h2o.mobility.mu,
                       nBvalue = plot.data$h2o.nBvalue.mu,
                       type = "Water",
                       stringsAsFactors = FALSE)

  df.het <- data.frame(name = cluster.name,
                       adn = plot.data$het.adn / plot.data$num.waters,
                       ahp = NA,
                       hbonds = plot.data$het.hbonds / plot.data$num.waters,
                       mobiity = plot.data$het.mobility.mu,
                       nBvalue = plot.data$het.nBvalue.mu,
                       type = "Other",
                       stringsAsFactors = FALSE)

  ##__ combine the rows -----
  df.bwe <- rbind(df.prot, df.h2o, df.het)

  ##_ modify the data to make the plots look nice -----
  mobility <- plot.data$call$mobility
  nBvalue <- plot.data$call$nBvalue
  df.bwe$mobiity[df.bwe$mobiity >= mobility] <- mobility + 0.1
  df.bwe$nBvalue[df.bwe$nBvalue >= nBvalue] <- nBvalue + 0.1

  ##_ rename the columns -----
  colnames(df.bwe) <- c("name", "Atomic Density", "Hydrophilicity",
                        "Hydrogen Bonds", "Mobility", "Norm B-value",
                        "type")

  ##_ melt the data.frame and order the values using factor() -----
  df.bwe.melt <- reshape2::melt(df.bwe, id.vars = c("name", "type") )
  df.bwe.melt$name <- factor(df.bwe.melt$name, levels=cluster.name)
  df.bwe.melt$type <- factor(df.bwe.melt$type, levels=c("Protein",
                                                        "Water",
                                                        "Other"))

  ##_ color -----
  bwe.colors <- c("#33a02c", ## medium green
                  "#1f78b4", ## medium blue
                  "#b2df8a"  ## light green
  )

  ##_ make the facet plot -----
  bwe.plots <- ggplot2::ggplot(data=df.bwe.melt,
                               aes_string(x = "name",
                                          y = "value",
                                          fill = "type")) +
    ggplot2::geom_bar(position = "dodge", stat = "identity") +
    ggplot2::scale_fill_manual(values = bwe.colors, name = "") +
    ggplot2::theme(axis.text.x = element_text(angle = 90,
                                              hjust = 1, vjust = 0.5)) +
    ggplot2::facet_grid(variable ~ ., scales = "free_y") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(x = "Water Cluster", y = "Mean Value") +
    ggplot2::scale_x_discrete(labels = cluster.numbering)

  ##_ return the plot -----
  bwe.plots

}


## bound water environment summary plot ----------------------------------------
#' @title Bound Water Environment Summary Plot
#' @description Mean bound water environment summary per percent conservation
#' @details Constructs a line plot with the bound water environment measures for
#'   the nearby protein and water atoms. The protein atomic density (ADN),
#'   hydrophilicity, mobility, normalized B-values, and potential hydrogen bonds
#'   are summarized for protein heavy atoms with 3.6 Angstroms along with the
#'   mobility, normalized B-values, and hydrogen bonds are summarized for the
#'   waters within 3.6 Angstroms of the protein and water atoms of interest,
#'   respectively. The raw values are scaled to values between 0 and 1 and
#'   plotted for each of the percent conservation available. Thus if there are
#'   ten structures being analyzed the percent conservation can range from 10 to
#'   100\% in 10\% increments. The protein related values are shown as solid
#'   lines and the water related values are shown as dotted lines.
#'
#'   _Interpreting the plot_
#'   - **dark green**: protein atom density
#'   - **medium green**: protein atom hydrophilicity
#'   - **green**: protein mobility
#'   - **pale green**: protein nBvalue
#'   - **light green**: protein hydrogen bonds
#'   - **dark blue**: water mobility
#'   - **medium blue**: water nBvalue
#'   - **blue**: water hydrogen bonds
#'
#'   This plot is based on Figure 3 of Sanschagrin and Kuhn (1998). Please note
#'   the B-value have been replaced with normalized B-values and hydrophilicity
#'   has been removed. Hydrophilicity was removed because the range between
#'   average hydrophilicity values for the percent conservations would likely be
#'   narrow. Due to the way scaling works, the lowest value is scaled to zero
#'   and the greatest value is scaled to one. Scaling the mean hydrophilicity
#'   values works against our goal of showing an overall tread and instead
#'   creates confusion about the values.
#'
#' @param data The `h2o.clusters.summary` data.frame from the `ClusterWaters`
#'   function containing the `nBvalue.mu` information. This data.frame is found
#'   within the `h2o.cluster.passed` and `h2o.cluster.all`
#' @param passed.waters Logical indicator to plot results for waters **passing**
#'   [Mobility()] and [NormalizedBvalue()] _**OR**_ using **all** waters within
#'   the `PDB` files.
#' @param title The title for the plot
#'
#' @export
#'
#' @import ggplot2
#' @import reshape2
#' @import scales
#' @importFrom stats aggregate
#'
#' @examples
#'   \dontrun{
#'   bwe.summary.plot <- BoundWaterEnvSummaryPlot(data=thrombin10.conservedWaters,
#'                                                passed.waters=TRUE,
#'                                                title="Bound Water Environment per Conservation")
#'  }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family plots
#'
#' @references
#'   Paul C Sanschagrin and Leslie A Kuhn. Cluster analysis of
#'   consensus water sites in thrombin and trypsin shows conservation between
#'   serine proteases and contributions to ligand specificity. _Protein
#'   Science_, 1998, **7** (_10_), pp 2054-2064.
#'   [DOI: 10.1002/pro.5560071002](http://doi.org/10.1002/pro.5560071002)
#'   [PMID: 9792092](http://www.ncbi.nlm.nih.gov/pubmed/9792092)
#'   [WatCH webpage](http://www.kuhnlab.bmb.msu.edu/software/watch/index.html)
#'
BoundWaterEnvSummaryPlot <- function(data,
                                     passed.waters = TRUE,
                                     title = "Bound Water Environment per Conservation") {

  ##_ get the PASSED or ALL waters plot data -----
  if (passed.waters == TRUE) {
    plot.data <- data$h2o.cluster.passed$h2o.clusters.summary
  } else {
    plot.data <- data$h2o.cluster.all$h2o.clusters.summary
  }

  ##_ construct the data.frame for the plot -----
  ##__ calculate the mean values per the percent conservation -----
  list.pct.conserved <- list(plot.data$pct.conserved)
  mean.pct.conserved <- aggregate(plot.data, list.pct.conserved, mean, na.rm = TRUE)
  ##__ get the needed values
  df.plot <- mean.pct.conserved[, c("pct.conserved",
                                    "prot.adn",
                                    "prot.mobility.mu", "prot.nBvalue.mu",
                                    "prot.hbonds",
                                    "mobility.mu", "nBvalue.mu", "h2o.hbonds")]
  ##__ normalize the values -----
  df.plot.norm <- apply(df.plot[, -1], 2, RescaleValues)
  ##__ construct the normalized data.frame
  df.plot.norm <- data.frame(pct.conserved = (df.plot[,1] / 100),
                             df.plot.norm,
                             stringsAsFactors = FALSE)

  ##_ information for the plot -----
  ##__ make easy to read legend labels -----
  legend.labels <- c("Protein Atom Density",
                     "Protein Mobility", "Protein nBvalues", "Protein H-bonds",
                     "Water Mobility", "Water nBvalues",
                     "Water H-bonds"
  )
  ##__ color values | colorbrewer2.org; single hue; 9 data classes; only 7 used -----
  prot.h2o.color <- c("#006d2c", ## dark green (protein atom density)
                      "#238b45", ## medium green (protein mobility)
                      "#41ab5d", ## green (protein nBvalue mu)
                      "#74c476", ## pale green (protein hydrogen bonds)
                      "#08519c", ## dark blue (water mobility mu)
                      "#2171b5", ## medium blue (water nBvalue mu)
                      "#4292c6"  ## blue (water hydrogen bonds)
  )
  ##__ define the line types -----
  prot.h2o.line.types <- c(rep("solid", 5), rep("dotted", 3))

  ##_ melt the data for ggplot2 -----
  ##__ protein values -----
  df.plot.prot <- df.plot.norm[, c("pct.conserved",
                                   "prot.adn",
                                   "prot.mobility.mu", "prot.nBvalue.mu",
                                   "prot.hbonds")]
  ##__ add the type column -----
  df.plot.prot <- data.frame(df.plot.prot, type = "protein",
                             stringsAsFactors = FALSE)
  ##__ melt the data.frame -----
  df.plot.prot.melt <- reshape2::melt(df.plot.prot,
                                      id.vars = c("pct.conserved", "type") )
  ##__ easy to read legend labels -----
  df.plot.prot.melt$variable <- factor(df.plot.prot.melt$variable,
                                       labels = legend.labels[1:4])
  ##__ water values -----
  df.plot.h2o <- df.plot.norm[, c("pct.conserved",
                                  "mobility.mu", "nBvalue.mu",
                                  "h2o.hbonds")]
  ##__ add the type column -----
  df.plot.h2o <- data.frame(df.plot.h2o, type = "water",
                            stringsAsFactors = FALSE)
  ##__ melt the data.frame -----
  df.plot.h2o.melt <- reshape2::melt(df.plot.h2o,
                                     id.vars = c("pct.conserved", "type") )
  ##__ easy to read legend labels -----
  df.plot.h2o.melt$variable <- factor(df.plot.h2o.melt$variable,
                                      labels = legend.labels[5:7])
  ##__ combine the protein and water data -----
  df.plot.norm.melt <- rbind(df.plot.prot.melt, df.plot.h2o.melt)

  ##_ make the plot -----
  bwe.plot <- ggplot2::ggplot(data = df.plot.norm.melt,
                              aes_string(x = "pct.conserved",
                                         y = "value",
                                         colour = "variable",
                                         linetype = "variable")) +
    ggplot2::geom_line(size = 2, lineend = "round") +
    ggplot2::scale_x_continuous(labels = scales::percent) +
    ggplot2::scale_colour_manual(values = prot.h2o.color,
                                 name = "",
                                 labels = legend.labels) +
    ggplot2::scale_linetype_manual(values = prot.h2o.line.types,
                                   name = "",
                                   labels = legend.labels) +
    ggplot2::labs(x = "Percent Conservation", y = "Scaled Values",
                  title = title) +
    ggplot2::theme(legend.position = "bottom")

  ##_ return the plot -----
  bwe.plot

}
