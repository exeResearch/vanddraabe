## AlignOverlap.R
##
## apr-22-2016 (exe) created
## apr-23-2016 (exe) working
## apr-25-2016 (exe) added output to Excel workbook
## may-04-2016 (exe) added removal of chains with minimum overlap from
##                   collection of acceptable chains
## oct-06-2016 (exe) ability to add results to existing workbook
## jan-20-2017 (exe) corrected formatting based on lintr
## jan-23-2017 (exe) update user provided parameter checks
## jan-23-2017 (exe) update documentation
## feb-09-2017 (exe) updated the documentation (rMarkdown)
## feb-27-2017 (exe) created CalcAlignOverlap function to aid testing
## feb-27-2017 (exe) updated the documentation
## mar-28-2017 (exe) removed option to clean structures
## jul-25-2017 (exe) updated documentation
## aug-08-2017 (exe) updated PASSED and FAILED messages in AlignOverlap()
## aug-10-2017 (exe) added @importFrom stats ...
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


## alignment overlap docs ---------------------------------------------
#' @title Alignment Overlap Check
#' @description Determine if two protein structures are aligned using
#'   C-alpha atoms.
#' @details Using the C-alpha atoms of two aligned proteins, the amount of
#'   atomic overlap is determined and the overlapped chains are written to
#'   individual PDB files in the `NAME_alignedGood` directory. The PDB files have
#'   the `PDBID_aligned_pruned.pdb` naming convention where the PDBID is the RCSB
#'   four-character identification code. Structures not meeting the user defined
#'   overlap ratio are written to the `NAME_alignedPoor` directory. The structures are written using the
#'   [bio3d::write.pdb()] function of the [bio3d] package.
#'
#' @param aligned.dir Directory with aligned structures
#' @param out.dir Directory prefix for the correctly and incorrectly aligned
#'   structures. The out.dir variable is also used to construct the Excel
#'   workbook with the summary of the alignment evaluation
#' @param ref.PDBid Reference structure PDB ID, four character ID, used to
#'   compare all other aligned structures
#' @param overlap The ratio of overlapping C-alpha atoms; default 0.70
#' @param removal The ratio of overlapping C-alpha atoms to remove a chain from
#'   a collection of chains passing the `overlap` requirement; default:
#'   0.10
#' @param CA.dist The minimum distance between C-alpha atoms for the two C-alpha
#'   atoms to be considered aligned; default: 1.25
#' @param filename The filename of the Excel workbook containing all the
#'   results from the analysis.
#'
#' @return
#'   This function returns:
#'   * **Overlapping sturctures**: PDB structures satisfying the overlap
#'     requirements are written to the `out.dir_alignedGood` directory
#'   * **Non-Overlapping sturctures**: PDB structures _not_ satisfying
#'     the overlap requirement are written to `out.dir_alignedPoor`
#'   * **AlignOverlap.summary**: `data.frame` of the information written
#'     to the Excel workbook
#'   * **call**: The user provided parameters for the function
#'
#' @export
#'
#' @import bio3d
#'
#' @examples
#'  \dontrun{
#'   ## example from the Thrombin vignette
#'   AlignOverlap(aligned.dir = "",
#'                out.dir = "OVERLAP",
#'                ref.PDBid = "1hai",
#'                overlap = 0.70, removal = 0.10,
#'                CA.dist = 1.25,
#'                filename = "Thrombin")
#'  }
#'
#' @family "Alignment Overlap"
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
AlignOverlap <- function(aligned.dir = "blast_fitlsq_1.0ang",
                         out.dir = "blast",
                         ref.PDBid = "1hai",
                         overlap = 0.70, removal = 0.10,
                         CA.dist = 1.25,
                         filename = "ProteinSystem") {

  ##----- check the input
  if ( !dir.exists(aligned.dir) ) {
    stop("Directory of aligned structures does not exist.")
  }
  if ( (overlap <= 0.0) || (overlap > 1.0) ) {
    stop("Overlap must be positive and greater than 0.0 and less than or equal to 1.0.")
  }
  if ( (removal <= 0.0) || (removal > 1.0) ) {
    stop("Removal must be positive and greater than 0.0 and less than or equal to 1.0.")
  }
  if (CA.dist <= 0.0) {
    stop("CA.dist must be positive and greater than 0.0.")
  }
  if (removal >= overlap){
    stop("Removal must be less than overlap.")
  }

  ##----- the provided call
  the.call <- match.call()


  ##----- create the date and time portion of the filenames
  current.time <- Sys.time()
  now.date.time <- FileTimeStamp(current.time)


  ##----- examine the the aligned.dir for structures
  PDB.path.files <- list.files(path = aligned.dir,
                               pattern = ".pdb",
                               full.names = TRUE)
  PDB.files <- basename(PDB.path.files)
  PDB.path <- dirname(PDB.path.files)
  PDB.path.file.clean <- paste(PDB.path, PDB.files, sep = "/")
  num.PDB.files <- length(PDB.files)
  PDBids <- as.vector(sapply(PDB.files,
                             function(x) {
                               paste(unlist(strsplit(x, ""))[1:4],
                                     sep = "",
                                     collapse = "")
                               }
                             )
                      )
  # PDBids <- PDBids[seq(1, (num.PDB.files * 3), by = 3)]


  ##----- get the reference structure
  ref.PDBid <- tolower(ref.PDBid)
  ref.PDB.file.idx <- grep(pattern=ref.PDBid, tolower(PDB.files))
  if (length(ref.PDB.file.idx) == 0) {
    stop.message <- paste0("The reference structure ", ref.PDBid,
                           " is not in the ", aligned.dir, " directory. ",
                           "Please make sure the reference structure is present.")
    stop(stop.message)
  }
  ref.PDB.file <- PDB.path.files[ref.PDB.file.idx]

  ref.PDB <- bio3d::read.pdb2(file = ref.PDB.file)
  ref.atom <- ref.PDB$atom[ref.PDB$atom$type == "ATOM", ]
  ref.ca <- ref.atom[ref.atom$elety == "CA", ]
  ref.num.atoms <- dim(ref.ca)
  ref.idc <- 1:ref.num.atoms[1]


  ##----- output directories
  ##--- create the names
  out.dir.good <- paste0(out.dir, "_alignedGood")
  out.dir.poor <- paste0(out.dir, "_alignedPoor")
  ##--- check if they exists
  out.dir.good.tf <- dir.exists(out.dir.good)
  out.dir.poor.tf <- dir.exists(out.dir.poor)

  ##--- create the good alignments directory
  if ( out.dir.good.tf ) {
    mess.txt <- paste0(out.dir.good,
                       " already exists. Contents will be overwritten or added to.")
    message(mess.txt)
  } else {
    dir.create(out.dir.good)
    mess.txt <- paste0("Created ", out.dir.good)
    message(mess.txt)
  }
  ##--- create the poor alignments directory
  if ( out.dir.poor.tf ) {
    mess.txt <- paste0(out.dir.poor,
                       " already exists. Contents will be overwritten or added to.")
    message(mess.txt)
  } else {
    dir.create(out.dir.poor)
    mess.txt <- paste0("Created ", out.dir.poor)
    message(mess.txt)
  }


  ##----- make data.frame for results
  df.results <- data.frame(PDBids = PDBids,
                           reference = FALSE,
                           chains.initial = NA,
                           chains.retained = NA,
                           passed = NA,
                           pct.overlap = NA,
                           pct.overlap.min = NA,
                           pct.overlap.max = NA,
                           stringsAsFactors = FALSE)


  ## COMPARE THE ALIGNMENT OF ALL TO A REFERENCE STRUCTURE ---------------------
  for (soi in 1:num.PDB.files) {

    ##----- is this structure the reference?
    if ( ref.PDBid == PDBids[soi] ) {
      ref.PDBid.tf <- TRUE
    } else {
      ref.PDBid.tf <- FALSE
    }

    ##----- get the structure of interest (soi)
    soi.PDB <- bio3d::read.pdb2(file = PDB.path.file.clean[soi])

    ##----- calculate the amount of alignment overlap
    AO.values <- CalcAlignOverlap(ref.num.atoms,
                                  ref.ca,
                                  ref.idc,
                                  soi.PDB,
                                  CA.dist)

    ##----- update the results data.frame
    ##--- overlapping C-alpha atoms percentage
    ratio.intersection <- AO.values$ratio.intersection
    passed.overlap.tf <- ratio.intersection >= overlap
    soi.chain.overlap.unique <- unique(AO.values$soi.chain.overlap)
    pct.intersec.clean.lt.removal.tf <- ratio.intersection > removal
    ##--- determine the chains to keep
    intersec.gt.removal.tf <- ratio.intersection > removal
    # pct.intersec.clean.lt.removal.tf <- ratio.intersection > removal
    ##--- make into percent with one and three decimal place
    pct.intersec.1 <- round(ratio.intersection * 100, digits = 1)
    pct.intersec.3 <- round(ratio.intersection * 100, digits = 3)
    if ( length(pct.intersec.3) > 1 ) {
      pct.intersec.3 <- paste(pct.intersec.3,
                              sep = "", collapse = ", ")
    }

    ##--- initial and retained chains
    soi.chain.uniq <- unique(AO.values$soi.chain)
    soi.chain.overlap.uniq <- unique(AO.values$soi.chain.overlap)

    chains.initial <- paste(soi.chain.uniq, collapse = ", ")
    chains.retained <- paste(soi.chain.overlap.uniq[intersec.gt.removal.tf], collapse = ", ")
    ##----- update the results data.frame
    df.results[soi, ] <- c(PDBids[soi], ## PDBids
                           ref.PDBid.tf, ## reference structure
                           chains.initial, ## chains.initial
                           chains.retained, ## chains.retains
                           any(passed.overlap.tf), ## passed (TRUE/FALSE)
                           pct.intersec.3, ## pct.overlap
                           min(pct.intersec.1), ## minimum pct.overlap
                           max(pct.intersec.1)) ## maximum pct.overlap


    ##----- write out the PDB of the aligned chains
    ##--- new name
    if ( any(passed.overlap.tf) ) {
      ##--- message
      mess.txt <- paste("  PASSED -->> ", PDBids[soi],
                        " Overlap: ", pct.intersec.1, "%", sep = "")
      message(mess.txt)

      ##--- create the filename with directory
      outfile.pass <- paste0(out.dir.good, "/",
                             PDBids[soi], "_aligned_pruned.pdb")

      ##--- write the overlapping portions to a PDB file
      soi.atom.all <- soi.PDB$atom
      soi.keep.tf <- soi.atom.all$chain %in% soi.chain.overlap.unique[pct.intersec.clean.lt.removal.tf]
      soi.overlapped <- soi.atom.all[soi.keep.tf, ]
      bio3d::write.pdb(pdb = NULL,
                       file = outfile.pass,
                       type = soi.overlapped[, "type"],
                       eleno = soi.overlapped[, "eleno"],
                       elety = soi.overlapped[, "elety"],
                       resid = soi.overlapped[, "resid"],
                       chain = soi.overlapped[, "chain"],
                       resno = soi.overlapped[, "resno"],
                       xyz = as.vector(unlist(t(soi.overlapped[, c("x", "y", "z")]))),
                       o = soi.overlapped[, "o"],
                       b = soi.overlapped[, "b"],
                       elesy = soi.overlapped[, "elesy"]
      )
    } else {
      ##--- message
      mess.txt <- paste("  <<<|||>>> FAILED -->> ", PDBids[soi],
                        " Overlap: ", pct.intersec.1, "%", sep = "")
      message(mess.txt)

      ##--- create the filename with directory
      outfile.fail <- paste0(out.dir.poor, "/", PDBids[soi], "_poorAlignment.pdb")

      ##--- write the poorly aligned structure to a PDB file
      bio3d::write.pdb(pdb = soi.PDB,
                       file = outfile.fail)
    }

  }


  ## WRITE RCSB DATA TO AN EXCEL WORKBOOK --------------------------------------
  ##----- construct workbook name
  filename.xlsx <- paste(filename, "_DATA_RESULTS.xlsx", sep = "")
  ##--- open existing excel workbook if available
  if ( file.exists(filename.xlsx) == TRUE ) {
    results.wb <- openxlsx::loadWorkbook(file = filename.xlsx)
  } else {  ##--- construct the workbook
    results.wb <- openxlsx::createWorkbook()
  }

  ##--- create the date and time portion of the filenames
  sheet.name <- paste(ref.PDBid, "AliOver", now.date.time, sep = "_")

  ##--- add the AlignOverlap summary
  results.wb <- oxAlignOverlapSheet(wb.name = results.wb,
                                    sheet.name = sheet.name,
                                    df = df.results)

  ##--- write the workbook
  openxlsx::saveWorkbook(results.wb, filename.xlsx, overwrite = TRUE)

  message("----- Results written to Excel workbook _____\n")

  ##----- done!
  message("done!!")

  ##----- return the results
  list(AlignOverlap.summary = df.results,
       call = the.call)

}


## calculate align overlap docs ------------------------------------------------
#' @title Calculate Alignment Overlap
#' @description Calculate the amount of alignment overlap between two protein
#'   structures using C-alpha atoms.
#' @details Using the C-alpha atoms of two aligned proteins, the amount of
#'   atomic overlap is determined. This function is within the [AlignOverlap]
#'   function.
#'
#'   This is a _**non-public**_ function and is _**NOT**_ available for general
#'   use. Please contact the author if you believe this function should be
#'   available for general use.
#'
#' @param ref.num.atoms Number of atoms in the reference structure
#' @param ref.ca PDB formatted `data.frame` containing only C-alpha atoms
#' @param ref.idc The indices of the reference structure atoms; from `1` to the
#'   number of atoms in the reference structure
#' @param soi.PDB The structure of interest (SoI) being compared to the
#'   reference structure. This is the full PDB structure read into `R` using the
#'   [bio3d::read.pdb2()] function
#' @param CA.dist The minimum distance between C-alpha atoms for the two C-alpha
#'   atoms to be considered aligned; default: 1.25
#'
#' @importFrom stats dist
#'
#' @return
#'   This function returns:
#'   * **ratio.intersection**: fraction of SOI overlapping with the reference structure
#'   * **soi.chain**: Chain letter designations for the aligned SOI
#'   * **soi.chain.overlap**: Unique chain letter designations for the aligned SOI
#'   These values are then used within the [AlignOverlap()] function to determine
#'   if the structures are adequately aligned.
#'
#' @family "Alignment Overlap"
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
CalcAlignOverlap <- function(ref.num.atoms,
                             ref.ca,
                             ref.idc,
                             soi.PDB,
                             CA.dist) {

  ##----- extract structure-of-interest data
  soi.atom.all <- soi.PDB$atom
  soi.atom <- soi.atom.all[soi.atom.all$type == "ATOM", ]
  soi.ca <- soi.atom[soi.atom$elety == "CA", ]
  soi.chain <- soi.ca$chain
  soi.size <- dim(soi.ca)

  ##----- construct the soi indices
  soi.idc <- (ref.num.atoms[1] + 1):(soi.size[1] + ref.num.atoms[1])


  ##----- create the coordinate data.frame
  df.refVSsoi <- rbind(ref.ca[, c("x","y","z")],
                       soi.ca[, c("x","y","z")],
                       stringsAsFactors = FALSE)

  ##----- calculate the distances
  refVSsoi.dist <- as.matrix(dist(df.refVSsoi,
                                  method = "euclidean",
                                  diag = TRUE, upper = TRUE))

  ##----- distances between reference and structure of interest
  refVSsoi.intersect <- refVSsoi.dist[ref.idc, soi.idc]

  ##----- C-alpha distances passing distance test
  refVSsoi.intersect.tf <- refVSsoi.intersect <= CA.dist
  ca.overlap.tf <- as.logical(colSums(refVSsoi.intersect.tf))

  ##----- overlapping chains
  soi.chain.overlap <- soi.chain[ca.overlap.tf]
  soi.chain.overlap.num <- as.integer(table(soi.chain.overlap))
  soi.chain.overlap.unique <- unique(soi.chain.overlap)
  soi.chain.overlap.full <- soi.chain[soi.chain %in% soi.chain.overlap.unique]
  soi.chain.overlap.full.num <- as.integer(table(soi.chain.overlap.full))

  ##----- overlapping C-alpha atoms percentage
  ratio.intersection <- soi.chain.overlap.num / soi.chain.overlap.full.num

  ##----- convert ratio values to percentages
#  pct.intersection.clean <- formatC(ratio.intersection * 100, digits = 3)
#  if ( length(pct.intersection.clean) > 1 ) {
#    pct.intersection.clean <- paste(pct.intersection.clean,
#                                    sep = "", collapse = ", ")
#  }

  ##----- return results
  list(ratio.intersection = ratio.intersection,
#       pct.intersection.clean = pct.intersection.clean,
       soi.chain = soi.chain,
       soi.chain.overlap = soi.chain.overlap
#       soi.chain.overlap.unique = soi.chain.overlap.unique,
#       pct.intersection.clean = pct.intersection.clean)
      )

}
