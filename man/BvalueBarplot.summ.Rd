% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{BvalueBarplot.summ}
\alias{BvalueBarplot.summ}
\title{B-value Summary Barplots}
\usage{
BvalueBarplot.summ(data)
}
\arguments{
\item{data}{The results from the \code{\link[=CleanProteinStructures]{CleanProteinStructures()}} function. Will use
the binned B-value data.}
}
\description{
B-value summary barplots for the PDB structures. The plots are
faceted and displays the binned B-value values for all the structures.
The counts are presented on a \code{log10} scale.
}
\examples{
  \dontrun{
  BvalueBarplot.summ(data)

  ##----- multiple pages
  library(ggforce)
  Bvalue.barplots.summary <- BvalueBarplot.summ(data)
  num.pages <- ceiling(nrow(data$Bvalue.counts) / 10)

  pdf(file="multiple_pages.pdf", height=11, width=8.5)
  for (page in seq_len(num.pages)) {
    print(Bvalue.barplots.summary +
          ggforce::facet_wrap_paginate(~PDBid,
                                       ncol = 2, nrow = 5, page = page) )
  }
  dev.off()
  }

}
\seealso{
Other plots: \code{\link{BoundWaterEnvPlots}},
  \code{\link{BoundWaterEnvSummaryPlot}},
  \code{\link{BvalueBarplot}},
  \code{\link{ClusterSummaryPlots}},
  \code{\link{MobNormBvalEvalPlots}},
  \code{\link{MobilityBarplot.summ}},
  \code{\link{MobilityBarplot}},
  \code{\link{OccupancyBarplot.summ}},
  \code{\link{OccupancyBarplot}},
  \code{\link{nBvalueBarplot}},
  \code{\link{normBvalueBarplot.summ}}
}
\author{
Emilio Xavier Esposito \email{emilio@exeResearch.com}
}
\concept{plots}
