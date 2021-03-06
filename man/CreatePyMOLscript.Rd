% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreatePyMOLscript.R
\name{CreatePyMOLscript}
\alias{CreatePyMOLscript}
\title{Create PyMOL Script File}
\usage{
CreatePyMOLscript(conservedWaters.data, passed.waters = TRUE,
  PDBid.ref = "1hai", LigResname.ref = NULL, hbond = 3.75,
  lig.carbon.color = "cyan", filename = "thrombin10")
}
\arguments{
\item{conservedWaters.data}{The \code{h2o.clusters.summary} data.frame from the
\code{\link[=ConservedWaters]{ConservedWaters()}} function containing the \code{nBvalue.mu}
information. This data.frame is found within the \code{h2o.cluster.passed}
and \code{h2o.cluster.all}}

\item{passed.waters}{Logical indicator to plot results for waters \strong{passing}
\code{\link[=Mobility]{Mobility()}} and \code{\link[=NormalizedBvalue]{NormalizedBvalue()}} \emph{\strong{OR}} using \strong{all} waters within
the \code{PDB} files.}

\item{PDBid.ref}{name for reference structure in PyMOL; \emph{e.g.},
\code{"1hai"}}

\item{LigResname.ref}{PDB residue code for reference ligand; \emph{e.g.}, \code{"0g6"}}

\item{hbond}{The minimum distance between hydrogen bond acceptor and donor;
default: \code{3.75}}

\item{lig.carbon.color}{One of the ten pre-defined carbon-color options using
PyMOL's \code{util.cbaX} command. The \code{X} represents the user defined
color of carbon atoms. \code{X} can be \code{g}: green; \code{c}: cyan;
\code{m}: magenta; \code{y}: yellow; \code{s}: salmon; \code{w}: grey;
\code{b}: slate; \code{o}: orange; \code{p}: purple; and \code{k}: pink;
default: \code{cyan}}

\item{filename}{Prefix for the PyMOL script files. It is probably best to use
the initial portion of the conserved waters PDB filename; \emph{e.g.},
\code{"thrombin10"}}
}
\description{
Create PyMOL script file to visualize conserved waters
}
\details{
The ability to visualize the conserved waters is important and their
surroundings is when exploring conserved water results.

Conserved waters within 6 Angstroms of the PyMOL identified ligands are
displayed. The conserved waters are colored based on their percent
conservation range using the same color scheme as the Percent Conservation
plot. Waters conserved less than 50\% are colored light grey, 50-69\% are
red, 70-79\% are dark red, 80-89\% are light blue, 90-99\% are medium blue,
and 100\% are dark blue. The conserved waters are labeled using their
ranking based on percent conseration.

This function creates \emph{\strong{two}} PyMOL script files; one with a
black background and another with a white background. The color of the
pocket residues is changed based on the background. The pocket residues are
colored light-grey for the black background and dark-grey for the white
background. The ligand is assigned the user-defined color for both
representations. Pocket residues -- and associated molecular surface -- are
defined as those within 5 Angstroms of the conserved waters. The depicted
cartoon representation is for residues within 15 Angstroms of the
ligand(s).

The potential hydrogen bonds are depicted between:
\itemize{
\item conserved waters and ligand: orange dashed line
\item conserved waters and protein: green dashed line
\item conserved waters: blue dashed line
}
}
\examples{
  \dontrun{
  current.time <- Sys.time()
  CreatePyMOLfile(PDBid.ref = "Thrombin_initial10_alignedGood/1hai_aligned_pruned.pdb",
                  PDBid.ref = "1hai",
                  LigResname.ref = "0g6",
                  conserved.waters = "Thrombin_initial10_ConservedWaters_PASSED_mar292017_1535.pdb",
                  hbond = 3.75,
                  lig.carbon.color = "cyan",
                  filename = "thrombin10_ConservedWaters_PASSED")
  }

}
\author{
Emilio Xavier Esposito \email{emilio@exeResearch.com}
}
\concept{Visualization}
