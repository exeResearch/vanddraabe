% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{TimeSpan}
\alias{TimeSpan}
\title{Time Span}
\usage{
TimeSpan(time.start)
}
\arguments{
\item{time.start}{The start time determined using the \code{\link[base:Sys.time]{base::Sys.time()}}}
}
\value{
character string of the calculation duration
}
\description{
Calculate the duration of a set of calculations.
}
\details{
Using the time a set of calculations started, the duration of the
calculations is returned.

\emph{\strong{NOTE}}: This is a non-public function.
}
\examples{
  \dontrun{
  time.start <- Sys.time() - 25   ## subtract 25 seconds from time.start
  TimeSpan(time.start)
  # [1] "00:00:25"
  }

}
\seealso{
Other utilities: \code{\link{ConservationSet}},
  \code{\link{DetermineChainsOfInterest}},
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
  \code{\link{UniqueAtomHashes}},
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
