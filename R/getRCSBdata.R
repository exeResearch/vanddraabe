## getRCSBdata.R
##
## oct-05-2016 (exe) created
## oct-06-2016 (exe) working
## feb-09-2017 (exe) updated the documentation (rMarkdown)
## mar-28-2017 (exe) fixed file copy error (charmatch to grep)
## jul-25-2017 (exe) updated documentation
## jul-31-2017 (exe) updated getRCSBdata() documentation
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


## getRCSBdata docs -----------------------------------------------
#' @title Clean RCSB Dataset
#' @description Clean the protein dataset based on quality values.
#' @details The provided protein models determined by X-ray crystallography
#'   and downloaded from the RCSB include structure quality measures.
#'   The resolution, rObservation, and rFree are the three commonly used
#'   and referenced evaluation measures.
#'
#'   The B-value normalization exclusion value is user defined within the main
#'   [ConservedWaters()] function but has a default value of 1.0.
#'
#' @param prefix Directory of aligned structures; string.
#' @param resolution Structures with a `resolution` value greater than this
#'   value are removed from analysis; default: `3.0`
#' @param rFree Structures with a `rFree` values greater than this value are
#'   removed from analysis; default: `0.26`
#' @param rObserved Structures with a `rObserved` values greater than this value
#'   are removed from analysis; default: 0.20
#' @param filename The filename prefix for the returned results. Default is
#'   "ProteinSystem"
#'
#' @examples
#'   \dontrun{
#'   proteins.info <- getRCSBdata(prefix="./thrombin_fitlsq_1.0ang/",
#'                                resolution=3.0, rFree=NULL, rObserved= 0.20,
#'                                filename="ProteinSystem")
#'   }
#'
#' @return
#'   This function returns:
#'   * **PDB.info**: RCSB provided information for all protein structures
#'   * **PDB.info.passed**: RCSB provided information for all protein
#'     structures _**passing**_ the user defined parameters
#'   * **PDB.info.rejected**: RCSB provided information for all protein
#'     structures _**failing**_ the user defined parameters
#'   * **call**: parameters provided by the user
#'   * **Excel workbook**: containing the PDB.info, PDB.info.passed, and
#'     PDB.info.rejected data as individual tabs
#'
#' @export
#'
#' @import bio3d
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family Clean RCSB dataset
#'
getRCSBdata <- function (prefix = "./alignTesting",
                              resolution = 3.0, rFree = 0.26, rObserved = 0.20,
                              filename = "ProteinSystem") {

  ##----- the provided call
  the.call <- match.call()

  ##----- create the date and time portion of the filenames
  current.time <- Sys.time()
  now.date.time <- FileTimeStamp(current.time)

  ##----- construct Excel sheetnames
  PDBs.info.all      <- paste("PDBsInfo_all", now.date.time, sep="_")
  PDBs.info.passed   <- paste("PDBsInfo_pass", now.date.time, sep="_")
  PDBs.info.rejected <- paste("PDBsInfo_reject", now.date.time, sep="_")

  ##----- rename user defined variables
  cutoff.resolution <- resolution
  cutoff.rObserved <- rObserved
  cutoff.rFree <- rFree

  ##----- check the user provided cutoff values
  if ( is.logical(cutoff.resolution) ) {
    cutoff.resolution <- NULL
  }
  if ( is.logical(cutoff.rObserved) ) {
    cutoff.rObserved <- NULL
  }
  if ( is.logical(cutoff.rFree) ) {
    cutoff.rFree <- NULL
  }

  ##----- was a resolution cutoff provided?
  if ( is.null(cutoff.resolution) == TRUE ) {
    message.res.cutoff <- paste("A maximum resolution value was NOT provided!",
                                "A default value of 3.0 will be used.",
                                "This value reflects the size of a water molecule.",
                                sep=" ")
    message(message.res.cutoff)
    cutoff.resolution <- 3.0
  }
  ##--- is the provided resolution cutoff greater than the default value of 3.0?
  if ( cutoff.resolution > 3.0 ) {
    message.res.cutoff <- paste("The typical resolution cutoff value is 3.0 is used",
                                "to resolve the backbone. The provided value of",
                                cutoff.resolution, "Angstroms will be used but",
                                "it is STRONGLY recommented to use a value of",
                                "3.0 Angstroms.", sep=" ")
    message(message.res.cutoff)
  }


  ##----- determine the number of cutoffs to include
  num.RCSB.cutoffs <- ( (!is.null(cutoff.resolution)) +
                         (!is.null(cutoff.rObserved)) +
                         (!is.null(cutoff.rFree)) )


  ##----- get list of PDB files within prefix
  pdb.location <- ReturnPDBfullPath(prefix)
  ##--- determine PDB IDs
  pdb.ids <- ExtractPDBids(pdb.location)
  ##--- number of PDB IDs
  num.pdb.structures.initial <- length(pdb.ids)

  ##----- remove structures with a resolution greater than the cutoff
  ##--- get the resolution, rObserved, and rFree values
  message("Please be patient... Getting PDB information from www.rcsb.org")
  pdbs.info.orig <- bio3d::pdb.annotate(pdb.ids, unique = TRUE)
  ##--- re-arrange the pdbs.information data.frame
  pdbs.info.col.order <- c("structureId", "resolution", "rObserved", "rFree",
                           "chainId", "experimentalTechnique", "ligandId",
                           "ligandName", "source", "scopDomain",
                           "classification", "compound", "title",
                           "citationAuthor", "journalName", "publicationYear",
                           "citation", "structureTitle", "depositionDate",
                           "structureMolecularWeight", "macromoleculeType",
                           "entityId", "sequence", "chainLength", "db_id",
                           "db_name")

  pdbs.info <- pdbs.info.orig[, pdbs.info.col.order]


  ##----- determine if the RCSB values meet the cutoff requirements
  pdbs.resolution <- as.numeric(pdbs.info$resolution)
  pdbs.rObserved <- as.numeric(pdbs.info$rObserved)
  pdbs.rFree <- as.numeric(pdbs.info$rFree)
  ##--- cutoff T/F vectors
  ##- resolution T/F
  if ( !is.null(cutoff.resolution) ){
    resolution.keep.tf <- pdbs.resolution <= cutoff.resolution
  } else {
    resolution.keep.tf <- rep(FALSE, num.pdb.structures.initial)
  }
  ##- rObserved T/F
  if ( !is.null(cutoff.rObserved) ){
    rObserved.keep.tf <- pdbs.rObserved <= cutoff.rObserved
  } else {
    rObserved.keep.tf <- rep(FALSE, num.pdb.structures.initial)
  }
  ##- rFree T/F
  if ( !is.null(cutoff.rFree) ){
    rFree.keep.tf <- pdbs.rFree <= cutoff.rFree
  } else {
    rFree.keep.tf <- rep(FALSE, num.pdb.structures.initial)
  }

  ##--- convert NA values in T/F vectors to FALSE
  resolution.keep.tf[is.na(resolution.keep.tf)] <- FALSE
  rObserved.keep.tf[is.na(rObserved.keep.tf)] <- FALSE
  rFree.keep.tf[is.na(rFree.keep.tf)] <- FALSE
  ##--- T/F vector to indicate the structures to keep
  reso.rObs.rFree.keep.tf <- (resolution.keep.tf +
                                rObserved.keep.tf +
                                rFree.keep.tf) >= num.RCSB.cutoffs
  reso.rObs.rFree.keep.num <- sum(reso.rObs.rFree.keep.tf)
  reso.rObs.rFree.remove.num <- sum(!reso.rObs.rFree.keep.tf)
  ##--- passed and removed PDB IDs
  pdb.ids.passed <- pdb.ids[reso.rObs.rFree.keep.tf]
  pdb.ids.removed <- pdb.ids[!reso.rObs.rFree.keep.tf]
  if ( length(pdb.ids.removed) == 0 ) {
    pdb.ids.removed <- NULL
  }


  ##----- create data.frames of PDB information based on PASSED and REJECTED
  df.pdbs.info.PASSED <- pdbs.info[reso.rObs.rFree.keep.tf, ]
  if ( reso.rObs.rFree.remove.num > 0 ){
    df.pdbs.info.REJECTED <- pdbs.info[!reso.rObs.rFree.keep.tf, ]
  } else {
    df.pdbs.info.REJECTED <- NULL
  }


  ## INFORM THE USER ABOUT REMOVED STRUCTURES ----------------------------------
  ##----- inform the user about removed structures based on resolution
  if ( sum(!resolution.keep.tf) > 0 ) {
    ##--- the message about removed structures
    message.reso.removed <- paste("The following", sum(!resolution.keep.tf),
                                  "structure(s) was NOT retained because its",
                                  "resolution is greater than the cutoff value of",
                                  sprintf("%.1f", cutoff.resolution),
                                  "Angstroms.", sep=" ")
    message(message.reso.removed)
    pdbs.removed.reso <- pdb.ids[!resolution.keep.tf]
    reso.removed.pdbs <- pdbs.resolution[!resolution.keep.tf]
    message.pdbs.removed <- paste(pdbs.removed.reso, "==>",
                                  sprintf("%.1f", reso.removed.pdbs),
                                  sep=" ", collapse="\n")
    message("PDBid    Resolution")
    message(message.pdbs.removed)
    message("\n")
  }


  ##----- inform the user about removed structures based on rObserved
  if ( (sum(!rObserved.keep.tf) > 0) & (!is.null(cutoff.rObserved)) ) {
    ##--- the message about removed structures
    message.rObserved.removed <- paste("The following", sum(!rObserved.keep.tf),
                                       "structure(s) was NOT retained because",
                                       "its R-observed is greater than the",
                                       "cutoff value of",
                                       sprintf("%.3f", cutoff.rObserved),
                                       sep=" ")
    message(message.rObserved.removed)
    pdbs.removed.rObs <- pdb.ids[!rObserved.keep.tf]
    rObs.removed.pdbs <- pdbs.rObserved[!rObserved.keep.tf]
    message.pdbs.removed <- paste(pdbs.removed.rObs, "==>",
                                  sprintf("%.3f", rObs.removed.pdbs),
                                  sep=" ", collapse="\n")
    message("PDBid    rObserved")
    message(message.pdbs.removed)
    message("\n")
  }
  if ( is.null(cutoff.rObserved) ) {
    message.removed <- paste("The R-observed cutoff is set to \"NULL\" and",
                             "is not a factor in evaluating structures for",
                             "removal.\n", sep=" ")
    message(message.removed)
  }


  ##----- inform the user about removed structures based on rFree
  if ( (sum(!rFree.keep.tf) > 0) & (!is.null(cutoff.rFree)) ) {
    ##--- the message about removed structures
    message.rFree.removed <- paste("The following", sum(!rFree.keep.tf),
                                   "structure(s) was NOT retained because its",
                                   "R-free is greater than the cutoff value of",
                                   sprintf("%.3f", cutoff.rFree), sep=" ")
    message(message.rFree.removed)
    pdbs.removed.rFree <- pdb.ids[!rFree.keep.tf]
    rFree.removed.pdbs <- pdbs.rFree[!rFree.keep.tf]
    message.pdbs.removed <- paste(pdbs.removed.rFree, "==>",
                                  sprintf("%.3f", rFree.removed.pdbs),
                                  sep=" ", collapse="\n")
    message("PDBid    rFree")
    message(message.pdbs.removed)
    message("\n")
  } else {
    pdbs.removed.rFree <- NULL
  }
  if ( is.null(cutoff.rFree) ) {
    message.removed <- paste("The R-free cutoff is set to \"NULL\" and is not",
                             "a factor in evaluating structures for removal.\n",
                             sep=" ")
    message(message.removed)
  }


  ## COPYING RETAINED (PASSED) AND REMOVED (REJECTED) STRUCTURES ---------------
  ##----- copy the structures that do NOT pass the resolution, rObserved,
  ##      and/or rFree cutoffs to the rejected folder
  if ( reso.rObs.rFree.remove.num > 0 ) {
    dir.pdbs.rejected <- paste0(filename, "_RCSB_rejected")
    dir.create(dir.pdbs.rejected, showWarnings = FALSE)
    pdb.ids.removed.idc <- unlist(lapply(X=pdb.ids.removed,
                                         FUN=grep,
                                         pdb.location))
    blah <- file.copy(from = pdb.location[pdb.ids.removed.idc],
                      to = dir.pdbs.rejected,
                      overwrite = FALSE, copy.date = TRUE, copy.mode = TRUE)
  }


  ##----- copy retained/used PDBs into the passed folder
  if ( reso.rObs.rFree.keep.num > 0 ) {
    dir.pdbs.passed <- paste0(filename, "_RCSB_passed")
    dir.create(dir.pdbs.passed, showWarnings = FALSE)
    pdb.ids.passed.idc <- unlist(lapply(X=pdb.ids.passed,
                                        FUN=grep,
                                        pdb.location))
    blah <- file.copy(from = pdb.location[pdb.ids.passed.idc],
                      to = dir.pdbs.passed,
                      overwrite = FALSE, copy.date = TRUE, copy.mode = TRUE)
  }


  ## WRITE RCSB DATA TO AN EXCEL WORKBOOK --------------------------------------
  ##----- construct workbook name
  filename.xlsx <- paste(filename, "_DATA_RESULTS.xlsx", sep="")
  ##--- open existing excel workbook if available
  if ( file.exists(filename.xlsx) == TRUE ) {
    results.wb <- openxlsx::loadWorkbook(file=filename.xlsx)
  } else {  ##--- construct the workbook
    results.wb <- openxlsx::createWorkbook()
  }

  ##--- RCSB/PDB structure information
  ##- ALL provided structures
  results.wb <- oxRCSBinfoSheet(wb.name=results.wb,
                                sheet.name=PDBs.info.all,
                                df=pdbs.info)
  ##- PASSED structures
  results.wb <- oxRCSBinfoSheet(wb.name=results.wb,
                                sheet.name=PDBs.info.passed,
                                df=df.pdbs.info.PASSED)
  ##- REJECTED structures
  if ( reso.rObs.rFree.remove.num > 0 ) {
    results.wb <- oxRCSBinfoSheet(wb.name=results.wb,
                                  sheet.name=PDBs.info.rejected,
                                  df=df.pdbs.info.REJECTED)
  }

  ##--- write the workbook
  openxlsx::saveWorkbook(results.wb, filename.xlsx, overwrite=TRUE)


  ## SUMMARY TO USER -----------------------------------------------------------
  message("\n----- getRCSBdata SUMMARY _____\n")
  message("getRCSBdata is DONE! \n")
  message(paste("RCSB information for each PDB structure was written to the ",
                "Excel workbook: ",
                filename.xlsx, "\n", sep="") )

  if ( reso.rObs.rFree.remove.num == 0 ) {
    summary.mess <- paste("All structures (", num.pdb.structures.initial,
                          ") PASSED the structure evaluation requirements ",
                          "and were copied to the \"",
                          dir.pdbs.passed, "\" folder.", sep="")
  } else {
    summary.mess <- paste(reso.rObs.rFree.keep.num,
                          " of the initial ", num.pdb.structures.initial,
                          " PDB structures from the provided alignment folder \"",
                          prefix, "\" PASSED the structure evaluation ",
                          "requirements and were copied to the \"",
                          dir.pdbs.passed, "\" folder. The ",
                          reso.rObs.rFree.remove.num,
                          " rejected structures were copied to the \"",
                          dir.pdbs.rejected, "\" folder.\n", sep="")
  }
  message(summary.mess)


  ## RETURN THE RESULTS --------------------------------------------------------
  list(PDB.info=pdbs.info,
       PDB.info.passed=df.pdbs.info.PASSED,
       PDB.info.rejected=df.pdbs.info.REJECTED,
       call=the.call
  )

}
