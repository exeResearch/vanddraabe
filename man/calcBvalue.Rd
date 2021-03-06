% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BoundWaterEnvironments.R
\name{calcBvalue}
\alias{calcBvalue}
\title{Calculate B-value}
\usage{
calcBvalue(rmsfValue)
}
\arguments{
\item{rmsfValue}{rmsf value calculated by \code{\link[bio3d:rmsf]{bio3d::rmsf()}}}
}
\value{
B-value (aka B-factor) in Angstroms^2^
}
\description{
Calculate the B-value for an atom.
}
\details{
The B-value (aka B-factor) is calcualted from the rmsf from a
collection of atoms. The rmsf is calculated using \code{\link[bio3d:rmsf]{bio3d::rmsf()}}.

\deqn{B-value = rmsf^{2} * 8 * {pi}^{2}}{B-value = rmsf^2 * 8 * pi^2}

The calculated B-values are returned within the \code{\link[=BoundWaterEnvironment]{BoundWaterEnvironment()}}
results and used to define the size of conserved waters for the depiction
of MDS conserved waters.
}
\examples{
  calcBvalue(rmsfValue=0.25)
  # [1] 4.935
  calcBvalue(rmsfValue=0.50)
  # [1] 19.74
  calcBvalue(rmsfValue=0.75)
  # [1] 44.41
  calcBvalue(rmsfValue=1.0)
  # [1] 78.96
  calcBvalue(rmsfValue=1.25)
  # [1] 123.4

}
\references{
Eaton E Lattman & Patrick J Loll. \emph{Protein Crystallography: A Concise
Guide}. Baltimore, Maryland, USA: The Johns Hopkins University Press, 2008.
QP551.L345 2008. ISBN: 978-0-8018-8808-3
\href{https://jhupbooks.press.jhu.edu/content/protein-crystallography}{website}
}
\seealso{
Other "Bound Water Environment": \code{\link{BoundWaterEnvironment.interact}},
  \code{\link{BoundWaterEnvironment.quality}},
  \code{\link{BoundWaterEnvironment}},
  \code{\link{Mobility}}, \code{\link{NormalizedBvalue}},
  \code{\link{calcNearbyHydrationFraction}},
  \code{\link{calcNumHydrogenBonds}}
}
\author{
Emilio Xavier Esposito \email{emilio@exeResearch.com}
}
\concept{"Bound Water Environment"}
