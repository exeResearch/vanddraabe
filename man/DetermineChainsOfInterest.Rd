% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChainsOfInterest.R
\name{DetermineChainsOfInterest}
\alias{DetermineChainsOfInterest}
\title{Determine Chains Of Interest}
\usage{
DetermineChainsOfInterest(chains.to.explore)
}
\arguments{
\item{chains.to.explore}{\emph{\strong{NOTE}}: \code{"first"} is alphabetically first. Thus
if the order within the original PDB file is \code{L} and then \code{H}, this
function will return \code{H} because it is alphebetically first.}
}
\value{
string indicating which chain designation (\emph{e.g.}, \code{"first"} chain,
\code{"all"} chains, or \code{"user"} defined) to include in the conserved water
analysis
}
\description{
Determine the chains identification
}
\details{
Standardizes user provided chain(s) of interest. This function
simply standardizes the user provided chains of interest. Acceptable values
are: - \strong{first}: alphabetically the first chain - \strong{all}: all chains
within a structure file - \strong{user defined}: a single letter or a set of
letters; \emph{e.g.}; \code{"A"} or \code{c("H", "L")}

\emph{\strong{NOTE}}: This is a \emph{\strong{non-public}} function and is \emph{\strong{NOT}} available
for general use. Please contact the author if you believe this function
should be available for general use.
}
\examples{
  \dontrun{
  DetermineChainsOfInterest("first")
  # [1] "first"
  DetermineChainsOfInterest("ALL")
  # [1] "all"
  DetermineChainsOfInterest("D")
  # [1] "user"
  DetermineChainsOfInterest(c("H", "L"))
  # [1] "user"
  DetermineChainsOfInterest("vanddraabe")
  # The provided chain ID VANDDRAABE is not valid and the first chain will
  # be used; likely chain A.
  # [1] "first"
  }

}
\seealso{
Other utilities: \code{\link{ConservationSet}},
  \code{\link{ExtractFileTimeStamp}},
  \code{\link{ExtractPDBids}}, \code{\link{FileTimeStamp}},
  \code{\link{HasXWaters}}, \code{\link{Nearby}},
  \code{\link{ProtHetWatIndices}},
  \code{\link{RescaleValues}},
  \code{\link{RetainChainsOfInterest}},
  \code{\link{ReturnPDBfullPath}},
  \code{\link{StandardizeAsparticAcidNames}},
  \code{\link{StandardizeCysteineNames}},
  \code{\link{StandardizeGlutamicAcidNames}},
  \code{\link{StandardizeHistidineNames}},
  \code{\link{StandardizeLysineNames}},
  \code{\link{TimeSpan}}, \code{\link{UniqueAtomHashes}},
  \code{\link{aaStandardizeNames}},
  \code{\link{getAtomTypeCounts}},
  \code{\link{getResTypeCounts}}, \code{\link{res2xyz}},
  \code{\link{resAtomType2AtomClass}},
  \code{\link{write.basic.pdb}},
  \code{\link{write.conservedWaters.pdb}}
}
\author{
Emilio Xavier Esposito \email{emilio@exeResearch.com}
}
\concept{utilities}
