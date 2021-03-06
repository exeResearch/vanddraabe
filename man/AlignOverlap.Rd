% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AlignOverlap.R
\name{AlignOverlap}
\alias{AlignOverlap}
\title{Alignment Overlap Check}
\usage{
AlignOverlap(aligned.dir = "blast_fitlsq_1.0ang", out.dir = "blast",
  ref.PDBid = "1hai", overlap = 0.7, removal = 0.1, CA.dist = 1.25,
  filename = "ProteinSystem")
}
\arguments{
\item{aligned.dir}{Directory with aligned structures}

\item{out.dir}{Directory prefix for the correctly and incorrectly aligned
structures. The out.dir variable is also used to construct the Excel
workbook with the summary of the alignment evaluation}

\item{ref.PDBid}{Reference structure PDB ID, four character ID, used to
compare all other aligned structures}

\item{overlap}{The ratio of overlapping C-alpha atoms; default 0.70}

\item{removal}{The ratio of overlapping C-alpha atoms to remove a chain from
a collection of chains passing the \code{overlap} requirement; default:
0.10}

\item{CA.dist}{The minimum distance between C-alpha atoms for the two C-alpha
atoms to be considered aligned; default: 1.25}

\item{filename}{The filename of the Excel workbook containing all the
results from the analysis.}
}
\value{
This function returns:
\itemize{
\item \strong{Overlapping sturctures}: PDB structures satisfying the overlap
requirements are written to the \code{out.dir_alignedGood} directory
\item \strong{Non-Overlapping sturctures}: PDB structures \emph{not} satisfying
the overlap requirement are written to \code{out.dir_alignedPoor}
\item \strong{AlignOverlap.summary}: \code{data.frame} of the information written
to the Excel workbook
\item \strong{call}: The user provided parameters for the function
}
}
\description{
Determine if two protein structures are aligned using
C-alpha atoms.
}
\details{
Using the C-alpha atoms of two aligned proteins, the amount of
atomic overlap is determined and the overlapped chains are written to
individual PDB files in the \code{NAME_alignedGood} directory. The PDB files have
the \code{PDBID_aligned_pruned.pdb} naming convention where the PDBID is the RCSB
four-character identification code. Structures not meeting the user defined
overlap ratio are written to the \code{NAME_alignedPoor} directory. The structures are written using the
\code{\link[bio3d:write.pdb]{bio3d::write.pdb()}} function of the \link{bio3d} package.
}
\examples{
 \dontrun{
  ## example from the Thrombin vignette
  AlignOverlap(aligned.dir = "",
               out.dir = "OVERLAP",
               ref.PDBid = "1hai",
               overlap = 0.70, removal = 0.10,
               CA.dist = 1.25,
               filename = "Thrombin")
 }

}
\seealso{
Other "Alignment Overlap": \code{\link{CalcAlignOverlap}}
}
\author{
Emilio Xavier Esposito \email{emilio@exeResearch.com}
}
\concept{"Alignment Overlap"}
