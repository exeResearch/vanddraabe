## utilities.R
##
## apr-10-2016 (exe) created
## may-04-2016 (exe) added FileTimeStamp() function
## oct-24-2016 (exe) added DOD (deuterated water) resid to ProtHetWatIndices
## nov-23-2016 (exe) added getAtomTypeCounts and getResTypeCounts functions
## nov-26-2016 (exe) updated documentation
## dec-14-2016 (exe) added examples
## feb-01-2017 (exe) fixed res2xyz bug
## feb-02-2017 (exe) added examples
## feb-03-2017 (exe) updated documentation
## feb-06-2017 (exe) alphabetized functions
## feb-07-2017 (exe) updated documentation
## mar-23-2017 (exe) added AtomType2AtomClass()
## apr-17-2017 (exe) added ConservationSet()
## apr-18-2017 (exe) updated documentation
## apr-27-2017 (exe) write.conservedWaters.pdb() writes calculated B-values
##                   instead of mean B-values
## jul-25-2017 (exe) updated documentation
## jul-28-2017 (exe) updated write.conservedWaters.pdb documentation
## jul-31-2017 (exe) updated FileTimeStamp() and ExtractFileTimeStamp() documentation
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


## aaStandardizeNames docs -----------------------------------------------------
#' @title Standardize Amino Acid Names
#' @description Standardize the various three-letter amino acid residue names.
#' @details The various three-letter amino acid residue names used to indicate
#'   protonation state or uncommon sidechain bonding (ligatation) are converted
#'   to the standard amino acid residue name.
#'
#' @param residue.names A vector of strings containing the three-letter residue
#'   names (strings)
#'
#' @return vector of _standardized_ amino acid residue names
#'
#' @examples
#'   residue.names <- c("HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP",
#'                      "CYS", "CYM", "CYX", "ASP", "ASH", "GLU", "GLH",
#'                      "LYS", "LYN")
#'   aaStandardizeNames(residue.names)
#'   # [1] "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "CYS" "CYS" "CYS"
#'   #     "ASP" "ASP" "GLU" "GLU" "LYS" "LYS"
#'
#' @export
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family utilities
#'
aaStandardizeNames <- function(residue.names) {

  ##----- find and convert various charged form names to standard names
  residue.names <- StandardizeHistidineNames(residue.names)
  residue.names <- StandardizeCysteineNames(residue.names)
  residue.names <- StandardizeAsparticAcidNames(residue.names)
  residue.names <- StandardizeGlutamicAcidNames(residue.names)
  residue.names <- StandardizeLysineNames(residue.names)

  ##----- return updated residue names
  return(residue.names)

}


## StandardizeAsparticAcidNames docs -------------------------------------------
#' @title Standardize Aspartic Acid Names
#' @description Standardize the protonated aspartic acid three-letter residue
#'   name to ASP.
#' @details The the protonated aspartic acid three-letter residue name (ASH) is
#'   converted to the standard "ASP" residue name. This function is part of the
#'   [aaStandardizeNames()].
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param residue.names A vector of strings containing the three-letter residue
#'   names (strings)
#'
#' @return vector of three-letter residue names with _standardized_
#'   aspartic acid residue names
#'
#' @examples
#'   \dontrun{
#'   residue.names <- c("HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP",
#'                      "CYS", "CYM", "CYX", "ASP", "ASH", "GLU", "GLH",
#'                      "LYS", "LYN")
#'   StandardizeAsparticAcidNames(residue.names)
#'   # [1] "HIS" "HID" "HIE" "HIP" "HSD" "HSE" "HSP" "CYS" "CYM" "CYX"
#'   #     "ASP" "ASP" "GLU" "GLH" "LYS" "LYN"
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family utilities
#'
StandardizeAsparticAcidNames <- function(residue.names) {

  ##----- find and replace "ASH" with "ASP"
  residue.names[which(residue.names %in% "ASH")] <- "ASP"

  ##----- return updated residue names
  return(residue.names)

}


## StandardizeCysteineNames docs -----------------------------------------------
#' @title Standardize Cysteine Names
#' @description Standardize the two three-letter cysteine residue names to CYS.
#' @details The two three-letter cysteine residue names used to indicate the
#'   different cystine states (CYM: deprotonated cysteine and CYX: no proton,
#'   neutral charge,  part of a disulfide bridge) are converted to the standard
#'   "CYS" (protonated) residue name. This function is part of the
#'   [aaStandardizeNames()].
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param residue.names A vector of strings containing the three-letter residue
#'   names (strings)
#'
#' @return vector of three-letter residue names with _standardized_
#'   cysteine residue names
#'
#' @examples
#'   \dontrun{
#'   residue.names <- c("HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP",
#'                      "CYS", "CYM", "CYX", "ASP", "ASH", "GLU", "GLH",
#'                      "LYS", "LYN")
#'   StandardizeCysteineNames(residue.names)
#'   # [1] "HIS" "HID" "HIE" "HIP" "HSD" "HSE" "HSP" "CYS" "CYS" "CYS"
#'   #     "ASP" "ASH" "GLU" "GLH" "LYS" "LYN"
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family utilities
#'
StandardizeCysteineNames <- function(residue.names) {

  ##----- other cysteine residue three-letter abbreviations
  CYS.other <- c("CYM", "CYX")

  ##----- find and replace with "CYS"
  residue.names[which(residue.names %in% CYS.other)] <- "CYS"

  ##----- return updated residue names
  return(residue.names)

}


## StandardizeGlutamicAcidNames docs -------------------------------------------
#' @title Standardize Glutamic Acid Names
#' @description Standardize the protonated glutamic acid three-letter residue
#'   name to GLU.
#' @details The the protonated glutamic acid three-letter residue name (GLH) is
#'   converted to the standard "GLU" residue name. This function is part of the
#'   [aaStandardizeNames()].
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param residue.names A vector of strings containing the three-letter residue
#'   names (strings)
#'
#' @return vector of three-letter residue names with _standardized_
#'   glutamic acid residue names
#'
#' @examples
#'   \dontrun{
#'   residue.names <- c("HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP",
#'                      "CYS", "CYM", "CYX", "ASP", "ASH", "GLU", "GLH",
#'                      "LYS", "LYN")
#'   StandardizeGlutamicAcidNames(residue.names)
#'   # [1] "HIS" "HID" "HIE" "HIP" "HSD" "HSE" "HSP" "CYS" "CYM" "CYX"
#'   #     "ASP" "ASH" "GLU" "GLU" "LYS" "LYN"
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family utilities
#'
StandardizeGlutamicAcidNames <- function(residue.names) {

  ##----- find and replace "GLH" with "GLU"
  residue.names[which(residue.names %in% "GLH")] <- "GLU"

  ##----- return updated residue names
  return(residue.names)

}


## StandardizeHistidineNames docs ----------------------------------------------
#' @title Standardize Histidine Names
#' @description Standardize the various three-letter histine residue names to
#'   HIS.
#' @details The various three-letter histidine residue names ("HID", "HIE",
#'   "HIP", "HSD", "HSE", "HSP") used to indicate the different protonation
#'   states are converted to the standard "HIS" residue name. This function is
#'   part of the [aaStandardizeNames()].
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param residue.names A vector of strings containing the three-letter residue
#'   names (strings)
#'
#' @return vector of three-letter residue names with _standardized_
#'   histidine residue names
#'
#' @examples
#'   \dontrun{
#'   residue.names <- c("HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP",
#'                      "CYS", "CYM", "CYX", "ASP", "ASH", "GLU", "GLH",
#'                      "LYS", "LYN")
#'   StandardizeHistidineNames(residue.names)
#'   # [1] "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "CYS" "CYM" "CYX"
#'   #     "ASP" "ASH" "GLU" "GLH" "LYS" "LYN"
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family utilities
#'
StandardizeHistidineNames <- function(residue.names) {

  ##----- other histidine residue three-letter abbreviations
  HIS.other <- c("HID", "HIE", "HIP", "HSD", "HSE", "HSP")

  ##----- find and replace with "HIS"
  residue.names[which(residue.names %in% HIS.other)] <- "HIS"

  ##----- return updated residue names
  return(residue.names)

}


## StandardizeLysineAcidNames docs ---------------------------------------------
#' @title Standardize Lysine Names
#' @description Standardize the de-protonated lysine three-letter residue name
#'   to LYS.
#' @details The the de-protonated lysine three-letter residue name (LYN) is
#'   converted to the standard "LYS" residue name. This function is part of the
#'   [aaStandardizeNames()].
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param residue.names A vector of strings containing the three-letter residue
#'   names (strings)
#'
#' @return vector of three-letter residue names with _standardized_ lysine
#'   residue names
#'
#' @examples
#'   \dontrun{
#'   residue.names <- c("HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP",
#'                      "CYS", "CYM", "CYX", "ASP", "ASH", "GLU", "GLH",
#'                      "LYS", "LYN")
#'   StandardizeLysineNames(residue.names)
#'   # [1] "HIS" "HID" "HIE" "HIP" "HSD" "HSE" "HSP" "CYS" "CYM" "CYX"
#'   #     "ASP" "ASH" "GLU" "GLH" "LYS" "LYS"
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family utilities
#'
StandardizeLysineNames <- function(residue.names) {

  ##----- find and replace "LYN" with "LYS"
  residue.names[which(residue.names %in% "LYN")] <- "LYS"

  ##----- return updated residue names
  return(residue.names)

}


## ConservationSet docs --------------------------------------------------------
#' @title Conservation Set
#' @description Assign the percent conservation to a "set#" for plotting.
#' @details Several of the plots color-code conserved water clusters based on
#'   percent conservation (see [ClusterSummaryPlots()] for color-coding) and is
#'   controlled by a `conserve.set` column. This function assigns less than 50\%
#'   conservation to `set0`, 50 to 69\% `set1`, 70 to 79\% `set2`, 80 to 89\%
#'   `set3`, 90 to 99\% `set4`, and eqaul to 100\% `set5`,
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param pct.conserved A `vector` from containing the [ConservedWaters()]
#'   function containing the percent conservation (`pct.conserved`)
#'
#' @return vector indicating the conservation set
#'
#' @examples
#'   \dontrun{
#'   pct.conserved <- c(100, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50,
#'                       45, 40, 35, 30, 25, 20, 15, 10, 10)
#'   ConservationSet(pct.conserved)
#'   # [1] "set5" "set4" "set4" "set3" "set3" "set2" "set2" "set1" "set1" "set1"
#'   # "set1" "set0" "set0" "set0" "set0" "set0" "set0" "set0" "set0" "set0"
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family utilities
#'
ConservationSet <- function(pct.conserved) {

  ##_ get percent conserved -----
  conserve.set <- rep("set0", length(pct.conserved))

  ##_ assign the corresponding set -----
  conserve.set[pct.conserved >=  50.0] <- "set1"
  conserve.set[pct.conserved >=  70.0] <- "set2"
  conserve.set[pct.conserved >=  80.0] <- "set3"
  conserve.set[pct.conserved >=  90.0] <- "set4"
  conserve.set[pct.conserved >= 100.0] <- "set5"

  ##_ return the data.frame -----
  return(conserve.set)

}


## ExtractPDBids docs ----------------------------------------------------------
#' @title Extract PDB IDs
#' @description Extract the four (4) character PDB identifier from the file name
#' @details The first four (4) characters of the file name -- typically the PDB
#'   ID is placed at the beginning of the file name -- are extracted and assumed
#'   to be the unique PDB ID.
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param pdb.location A collection of string values with the complete
#'   (normalized) path for each PDB file within the provided directory/folder
#'   obtained with the [ReturnPDBfullPath()].
#'
#' @return a vector of strings containing the PDB identifiers for the protein
#'   structures
#'
#' @examples
#'   \dontrun{
#'   ExtractPDBids("1hai.pdb")
#'   # [1] "1hai"
#'   ExtractPDBids("/home/someuser/pdbs/1hai.pdb")
#'   # [1] "1hai"
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family utilities
#'
ExtractPDBids <- function(pdb.location) {

  ##----- only the file name
  pdb.files <- basename(pdb.location)

  ##----- determine PDB IDs
  pdb.ids <- substring(pdb.files, 1, 4)

  ##----- return the PDB IDs
  return(pdb.ids)

}



## extract filename time stamp -------------------------------------------------
#' @title Extract Filename Time Stamp
#' @description Extract date & time stamp from a file
#' @details Create a date-time string to append to filenames to try and make
#'   them unique. The date-time string has the format
#'   month-day-year_hour-minute; for example, May 4, 2016 at 12:34pm is
#'   represented as `may042016_1234`.
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param filename String of the file name to extract the _**FileTimeStamp**_
#'   information
#'
#' @return A string with the date and time.
#'
#' @examples
#'   \dontrun{
#'   filename <- "ConservedWaters_PASSED_may042016_1234.pdb"
#'   ExtractFileTimeStamp(filename)
#'   # [1] may042016_1234
#'   }
#'
#' @family utilities
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
ExtractFileTimeStamp <- function(filename) {

  ##----- remove file extension
  main.name <- unlist(strsplit(x = filename, split = ".", fixed = TRUE))[1]

  ##----- split filename based on underscores "_"
  main.name.parts <- unlist(strsplit(x = main.name, split = "_"))

  ##----- number of filename sections/parts
  main.name.num <- length(main.name.parts)

  ##----- FileTimeStamp should be the last two parts (date & time)
  date.time <- main.name.parts[c(main.name.num - 1, main.name.num)]

  ##----- paste the date & time together
  now.date.time <- paste(date.time, collapse = "_")

  ##----- return original FileTimeStamp
  return(now.date.time)

}


## filename time stamp ---------------------------------------------------------
#' @title Filename Time Stamp
#' @description Date-time string to make file names unique
#' @details Create a date-time string to append to filenames to try and make
#'   them unique. The date-time string has the format
#'   month-day-year_hour-minute; for example, May 4, 2016 at 12:34pm is
#'   represented as `may042016_1234`.
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param current.time The current time determined with [base::as.POSIXct()]
#'
#' @return A string with the date and time.
#'
#' @examples
#'   \dontrun{
#'   current.time <- as.POSIXct("2016-05-04 12:34:56.78", tz = "UTC")
#'   FileTimeStamp(current.time)
#'   # [1] may042016_1234
#'   }
#'
#' @family utilities
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
FileTimeStamp <- function(current.time) {

  ##----- create the time stamp
  file.time.stamp <- paste(tolower(format(current.time, "%b%d%Y")),
                           format(current.time, "%H%M"),
                           sep = "_")

  ##----- return the time stamp
  file.time.stamp

}


## getAtomTypeCounts docs ------------------------------------------------------
#' @title Get AtomType Counts
#' @description Counts the number of AtomTypes within the provided string.
#' @details This is a wrapper using the [base::table()] function. The
#'   vector of AtomTypes (`strings`) are passed to the function,
#'   non-standard AtomTypes are removed, the AtomTypes are counted, and the
#'   counts are ordered based on the `names.AtomTypes` constant.
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param atom.types A vector of strings containing a combination of the 167
#'   AtomTypes.
#'
#' @return a vector of numbers indicating the counts of each AtomType. The
#'   vector is ordered based on the `names.AtomTypes` with AtomTypes not
#'   included assigned a value of zero (0).
#'
#' @examples
#'   \dontrun{
#'   set.seed(13)
#'   num.AtomTypes <- sample(1:10, 30, replace = TRUE)
#'   atom.types <- rep(sample(names.res.AtomTypes, 30), num.AtomTypes)
#'   getAtomTypeCounts(atom.types)
#'   # [1]  0 0 0 1 0 0 0 0 8 0 7 0 0 1 0 0 0 1 0 0 1 0 4 0 0 0 0 0 0 0 0 0 4
#'   #      0 0 9 0 0 0 0 0 0 0 0 0 0 0 5 0 6 0 0 0 0 0 0 0 0 6 0 0 0 0 0 7 0
#'   #      0 3 0 0 8 0 0 6 0 0 0 0 0 0 0 0 2 0 0 0 0 0 10 0 0 0 0 6 7 0 0 6
#'   #      0 0 7 0 0 0 0 0 0 0 0 9 0 0 0 0 0 0 0 0 0 0 9 0 0 0 0 0 0 0 0 0 0
#'   #      0 1 0 0 0 0 0 0 4 0 0 0 6 0 0 0 0 0 9 0 4 0 0 0 0 0 7 0 0 0 0 0 0
#'   #      0 0 0
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family utilities
#'
getAtomTypeCounts <- function(atom.types) {

  ##----- create empty vector for counts
  counts.vector <- rep(0, length(names.res.AtomTypes))

  ##----- retain only valid AtomTypes
  atom.types.cleaned <- atom.types[which(atom.types %in% names.res.AtomTypes)]

  ##----- count the number of each AtomType
  atom.types.counts <- t(table(atom.types.cleaned))

  ##----- determine the column order for the evaluated AtomTypes
  atom.types.counts.idc <- which(names.res.AtomTypes %in%
                                   colnames(atom.types.counts))

  ##----- add AtomTypes counts to the empty count vector
  counts.vector[atom.types.counts.idc] <- atom.types.counts

  ##----- return values
  return(counts.vector)

}


## getResTypeCounts docs -------------------------------------------------------
#' @title Get ResType Counts
#' @description Counts the number of ResType within the provided string.
#' @details This is a wrapper for the [base::table()] function. The
#'   vector of ResType are passed to the function, non-standard ResType are
#'   removed, the ResType are counted, and the counts are ordered based on the
#'   `names.ResTypes` constant.
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param res.types A vector of strings containing a combination of the 20
#'   ResTypes.
#'
#' @return a vector of numbers indicating the counts of each ResType. The
#'   vector is ordered based on the `names.ResTypes` with ResTypes not
#'   included assigned a value of zero (0).
#'
#' @examples
#'   \dontrun{
#'   set.seed(13)
#'   num.ResTypes <- sample(1:10, 20, replace = TRUE)
#'   res.types <- rep(names.residues, num.ResTypes)
#'   getResTypeCounts(res.types)
#'   # [1] 8  3  4  1 10  1  6  8  9  1  7  9  9  6  6  4  4  6  9  7
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family utilities
#'
getResTypeCounts <- function(res.types) {

  ##----- create empty vector for counts
  counts.vector <- rep(0, length(names.residues))

  ##----- retain only valid AtomTypes
  res.types.cleaned <- res.types[which(res.types %in% names.residues)]

  ##----- count the number of each AtomType
  res.types.counts <- t(table(res.types.cleaned))

  ##----- determine the column order for the evaluated AtomTypes
  res.types.counts.idc <- which(names.residues %in% colnames(res.types.counts))

  ##----- add AtomTypes counts to the empty count vector
  counts.vector[res.types.counts.idc] <- res.types.counts

  ##----- return values
  return(counts.vector)

}


## HasXWaters docs -------------------------------------------------------------
#' @title Has "X" Waters
#' @description Determines if PDB structure has water molecules.
#' @details Determine if the PDB structure has at least the user defined number
#'   of water oxygen atoms. The number of water oxygen atoms is returned along
#'   with a logical value indicating if the structure satisfies the user
#'   defined minimum.
#'
#'   Waters are identified using the three water three-letter residue
#'   names: HOH, WAT, and DOD.
#'
#' @param atoms.oi.resid vector of character strings containing the standardized
#'   three-letter amino acid residue names
#' @param min.num.h2o numeric value indicating the minimum number of water
#'   molecules required to return a `TRUE` logical value
#'
#' @return logical indicating if the PDB structure has the minimum user defined
#'   number of waters
#' @return numeric value indicating the number of water oxygen atoms within the
#'   PDB structure
#'
#' @export
#'
#' @examples
#'   resids <- c("ALA", "HOH", "WAT", "ALA", "HOH", "DOD", "ALA", "HOH")
#'   HasXWaters(resids, min.num.h2o = 4)
#'   # $has.h2o.tf
#'   # [1] TRUE
#'   #
#'   # $num.water
#'   # [1] 5
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family utilities
#'
HasXWaters <- function(atoms.oi.resid, min.num.h2o = 20) {

  ##----- possible water residue id values
  water.resids <- c("HOH", "WAT", "DOD")

  ##----- determine if waters are present
  water.tf <- atoms.oi.resid %in% water.resids
  num.water <- sum(water.tf, na.rm = TRUE)
  has.h2o.tf <- num.water >= min.num.h2o

  ##----- return TRUE or FALSE
  list(has.h2o.tf = has.h2o.tf,
       num.water = num.water)

}


## nearby docs -----------------------------------------------------------------
#' @title Nearby
#' @description Determine the entities near the entity of interest using a
#'   distance matrix.
#' @details Identify the entity, or entities, near an entity or collection of
#'   entites of interest. The previously calculated distance matrix, set of
#'   indicies, and a user defined radius are required.
#'
#'   _**NOTE**_: This function is designed to work with
#'   [BoundWaterEnvironment()] and the [base::apply()] function processing rows
#'   (the `MARGIN = 1` option). For this reason it is **NOT** a public function.
#'
#' @param distances Vector of distance values; see above note.
#' @param set.idc Vector of indices (as integers) indicating the entities of
#'   interest. This set of entities corresponds to the columns of the distance
#'   matrix because the provided distance matrix should be square. No check is
#'   performed on the squareness of the distance matrix because it is calculated
#'   _**within**_ the ConservedWaters function.
#' @param radius Numerical value indicating the distance to look for neighboring
#'   entities; default: 3.6
#'
#' @examples
#'   \dontrun{
#'   ##----- determine atom indices
#'   ProtHetWat.idc <- ProtHetWatIndices(thrombin.1hai$atom)
#'   prot.idc <- ProtHetWat.idc$prot.idc
#'   het.idc <- ProtHetWat.idc$het.idc
#'   h2o.idc <- ProtHetWat.idc$h2o.idc
#'
#'   ##----- calculate the distances
#'   atoms.dist <- as.matrix(dist(thrombin.1hai$atom[, c("x","y","z")],
#'                                method = "euclidean",
#'                                diag = TRUE, upper = TRUE))
#'   diag(atoms.dist) <- NA
#'   atom.idc <- sort(c(prot.idc, het.idc, h2o.idc))
#'   atoms.dist <- atoms.dist[atom.idc, atom.idc]
#'
#'   ##----- determine nearby atoms
#'   nearby.prot.idc <- Nearby(distances = atoms.dist[h2o.idc[1], ],
#'                             set.idc = prot.idc,
#'                             radius = 3.6)
#'   nearby.prot.idc
#'   # [1] 571
#'   atoms.dist[h2o.idc[1], nearby.prot]
#'   # [1] 3.571
#'   }
#'
#' @return Vector of indicies.
#'
#' @family utilities
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
Nearby <- function(distances, set.idc, radius = 3.6) {

  ##----- determine nearby indices
  ##      designed to work with apply(FUNCTION, MARGIN=1)
  nearby.idc <- set.idc[which(distances[set.idc] <= radius)]

  ##----- return nearby indicies
  return(nearby.idc)

}


## protein, het, and water atom indices ----------------------------------------
#' @title Protein, HET, and Water Atom Indices
#' @description Indices for the protein, HET-atom, and water atoms
#' @details Returns individual numerical vectors for the protein, HET-atom, and
#'   water atoms from the atom [base::data.frame()] of a PDB.
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param data The atom data.frame of the PDB read into the `R` session
#'   using the function [bio3d::read.pdb()].
#'
#' @return Individual vectors for the indices of the protein, HET-atom, and
#'   water atoms for a PDB file.
#'
#' @examples
#'   \dontrun{
#'   ProtHetWatIndices(thrombin.1hai$atom[c(1:10, 2341:2350, 2385:2394), ])
#'   # $prot.idc
#'   # [1]  1  2  3  4  5  6  7  8  9 10
#'   #
#'   # $het.idc
#'   # [1] 11 12 13 14 15 16 17 18 19 20
#'   #
#'   # $h2o.idc
#'   # [1] 21 22 23 24 25 26 27 28 29 30
#'   }
#'
#' @family utilities
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
ProtHetWatIndices <- function(data) {

  ##----- type classification
  type.ATOM.tf <- data$type == "ATOM"
  type.HETAM.tf <- data$type == "HETATM"

  ##----- residue name
  resid.HOH.tf <- data$resid == "HOH"
  resid.DOD.tf <- data$resid == "DOD"
  resid.WAT.tf <- data$resid == "WAT"
  resid.h2o.tf <- as.logical(resid.HOH.tf + resid.DOD.tf + resid.WAT.tf)

  ##----- oxygen atoms
  elety.OXY.tf <- data$elety == "O"

  ##----- protein atom indices
  prot.idc <- which( (type.ATOM.tf + !resid.h2o.tf)  == 2)

  ##----- het-atom indices
  het.idc <- which( (type.HETAM.tf + !resid.h2o.tf) == 2)

  ##----- water atom indices
  h2o.idc <- which( (type.HETAM.tf + resid.h2o.tf + elety.OXY.tf) == 3)


  ##----- return the information
  list(prot.idc = prot.idc,
       het.idc = het.idc,
       h2o.idc = h2o.idc
  )

}


## RescaleValues docs ------------------------------------------------------------
#' @title Rescale Values
#' @description Rescales provided vector of values to a user defined range.
#' @details Rescale the values to a new user defined range.
#'
#' @param data A vector of numerical values to be rescaled
#' @param newMin A numerical value indicating the new minimum value; default: 0
#' @param newMax A numerical value indicating the new maximum value; default: 1
#'
#' @return vector of rescaled numerical values.
#'
#' @examples
#'   RescaleValues(0:10, newMin = 0, newMax = 1)
#'   # [1] 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
#'
#' @export
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family utilities
#'
RescaleValues <- function(data, newMin = 0, newMax = 1) {

  ##----- make sure newMin is less than newMax
  if ( newMin >= newMax ) {
    mess <- paste("newMin value (", newMin,
                  ") needs to be less than newMax (", newMax, ") value",
                  sep="")
    stop(mess)
  }

  ##----- determine the range
  MinMax <- range(data, na.rm = TRUE, finite = TRUE)

  ##----- rescale the data
  data.rescaled <- newMin + (data - MinMax[1]) *
    ( (newMax - newMin) / (MinMax[2] - MinMax[1]) )

  ##----- return to user
  return(data.rescaled)

}


## residue indices to coordinate indices docs ----------------------------------
#' @title Residue Indices to Coordinate Indices
#' @description Return the coordinate indices for the provided residue indices.
#' @details Using the residue indices of the atoms
#'   [base::data.frame()] (_e.g._, `pdb$atom`) determine the
#'   coodinate indices of the residue atoms (_e.g._, `pdb$xyz`).
#'
#' @param res.idc Indicies of residues to convert to coordinate indices
#'
#' @return Vector of coordinate indicies to be applied to `pdb$xyz`
#'
#' @examples
#'   res.idc <- c(5:10)
#'   res2xyz(res.idc)
#'   # [1] 13 15 15 16 18 18 19 21 21 22 24 24 25 27 27 28 30
#'
#' @export
#'
#' @family utilities
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
res2xyz <- function(res.idc) {

  ##----- the internal function...
  z.idc <- res.idc * 3

  xyz.idc.values <- sort(c((z.idc - 2),
                           (z.idc - 1),
                           z.idc) )

  ##----- return needed indices
  return(xyz.idc.values)

}


## AtomType to AtomClass docs --------------------------------------------------
#' @title Convert Residue-AtomType to AtomType Class
#' @description Converts the residue-AtomType to AtomType Class.
#' @details See examples...
#'
#' @param resAT residue and AtomType; _e.g._, `"LYS NZ"`
#'
#' @return A string with the AtomType's class:
#'
#'   - Nitrogen
#'   - Nitrogen (+)
#'   - Oxygen
#'   - Oxygen (-)
#'   - Carbon
#'   - Sulfur
#'
#' @examples
#'   resAtomType2AtomClass(resAT="LYS NZ")
#'   # [1] "Nitrogen (+)"
#'   resAtomType2AtomClass(resAT="GLU N")
#'   # [1] "Nitrogen"
#'   resAtomType2AtomClass(resAT="VAL O")
#'   # [1] "Oxygen"
#'   resAtomType2AtomClass(resAT="ASP OD2")
#'   # [1] "Oxygen (-)"
#'   resAtomType2AtomClass(resAT="GLN CA")
#'   # [1] "Carbon"
#'   resAtomType2AtomClass(resAT="CYS SG")
#'   # [1] "Sulfur"
#'
#' @export
#'
#' @family utilities
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
resAtomType2AtomClass <- function(resAT) {

  ##----- make input UPPER CASE
  resAT <- toupper(resAT)

  ##----- determine element
  ##--- split reside-AtomType
  residue.AtomType <- unlist(strsplit(x=resAT, split=" "))
  ##--- get the element
  AtomType.element <- substr(residue.AtomType[2], 1, 1)

  ##----- nitrogen atom classes
  if ( AtomType.element == "N") {
    ##--- nitrogen neutral
    if ( any(resAT %in% names.resATs.nitro.neut)  ||
         (residue.AtomType == "UNK") ) {
      ATclass <- "Nitrogen"
    }

    ##--- nitrogen positive
    if ( any(resAT %in% names.resATs.nitro.pos) ) {
      ATclass <- "Nitrogen (+)"
    }

    ##--- return AtomType Class
    return(ATclass)
  }

  ##----- oxygen atom classes
  if ( AtomType.element == "O") {
    ##--- oxygen neutral
    if ( any(resAT %in% names.resATs.oxy.neut) ||
         (residue.AtomType == "UNK") ) {
      ATclass <- "Oxygen"
    }

    ##--- oxygen negative
    if ( any(resAT %in% names.resATs.oxy.neg) ) {
      ATclass <- "Oxygen (-)"
    }

    ##--- return AtomType Class
    return(ATclass)
  }

  ##----- carbon atoms
  if ( AtomType.element == "C" ) {
    ATclass <- "Carbon"

    ##--- return AtomType Class
    return(ATclass)
  }

  ##----- carbon atoms
  if ( AtomType.element == "S" ) {
    ATclass <- "Sulfur"

    ##--- return AtomType Class
    return(ATclass)
  }

}


## ReturnPDBfullPath docs ----------------------------------------------------
#' @title Return PDB Full Path
#' @description Determine the full path of the PDB files and return the complete
#'   path of each file within the provided directory.
#' @details The complete path of the PDB file(s) in the user provided
#'   `prefix` is returned.
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param prefix The directory with the PDB files of interest; _e.g._,
#'   ProteinSystem_Aligned
#'
#' @return collection of string values with the complete (normalized) path for
#'   each PDB file within the provided directory/folder.
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family utilities
#'
ReturnPDBfullPath <- function(prefix){

  ##----- get list of PDB files within prefix
  prefix <- normalizePath(prefix)
  pdb.files <- list.files(prefix, pattern = "pdb")

  ##----- directory and PDB files
  pdb.location <- paste(prefix, pdb.files, sep = "/")

  ##----- return full path and file name
  return(pdb.location)

}


## TimeSpan docs ---------------------------------------------------------------
#' @title Time Span
#' @description Calculate the duration of a set of calculations.
#' @details Using the time a set of calculations started, the duration of the
#'   calculations is returned.
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param time.start The start time determined using the [base::Sys.time()]
#'
#' @return character string of the calculation duration
#'
#' @examples
#'   \dontrun{
#'   time.start <- Sys.time() - 25   ## subtract 25 seconds from time.start
#'   TimeSpan(time.start)
#'   # [1] "00:00:25"
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family utilities
#'
TimeSpan <- function(time.start) {

  ##----- convert the start time
  time.start <- as.POSIXct(time.start)

  ##----- calculate the duration
  time.diff <- difftime(Sys.time(), time.start, units = "secs")

  ##---- return the duration for the calculation in HH:MM:SS
  format(.POSIXct(time.diff, tz = "GMT"), "%H:%M:%S")

}

## UniqueAtomHashes docs -------------------------------------------------------
#' @title Create Unique Atom Hashes
#' @description Constructs unique atom hashes from the provided
#' @details Using atom specific identifiers from a PDB-like formatted
#'   data.frame, unique atom hashes are constructed. The identifiers are
#'   separated by a user-defined separator, the default separator is an
#'   underscores ("_"), and the constructed hashes are returned as a vector.
#'
#'   Select a separator to allow easy splitting of the  the unique atom hashes
#'   using the [base::strsplit()] function to access the individual components.
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param atoms.oi A data.frame containing the common PDB information in columns
#' @param cols.oi A vector of column names to be used in the construction of the
#'   unique atom hashes
#' @param separator A single character string to separate the atom specific
#'   identifiers. Acceptable separators include: _ (default), -, +, ., :, |, " "
#'   (space), and "" (no separator).
#'
#' @return a vector of strings containing the unique atom hashes
#'
#' @examples
#'   \dontrun{
#'   atoms.oi <- thrombin.1hai$atom[1:10, ]
#'   cols.oi <- c("elety", "resid", "chain", "resno")
#'   UniqueAtomHashes(atoms.oi, cols.oi, separator = "_")
#'   # [1] "N_THR_L_1"   "CA_THR_L_1"  "C_THR_L_1"   "O_THR_L_1"   "CB_THR_L_1"
#'   #     "OG1_THR_L_1" "CG2_THR_L_1" "N_PHE_L_1"   "CA_PHE_L_1"  "C_PHE_L_1"
#'
#'   UniqueAtomHashes(atoms.oi, cols.oi, separator = "!")
#'   # The provided separator "!" is not acceptable. The default separator "_" is being used.
#'   #  [1] "N_THR_L_1"   "CA_THR_L_1"  "C_THR_L_1"   "O_THR_L_1"   "CB_THR_L_1"
#'   #      "OG1_THR_L_1" "CG2_THR_L_1" "N_PHE_L_1"   "CA_PHE_L_1"  "C_PHE_L_1"
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family utilities
#'
UniqueAtomHashes <- function(atoms.oi, cols.oi, separator = "_") {

  ##----- check separator
  separator.ok <- c("_", "-", "+", ".", ":", "|", " ", "")
  ##----- separator check
  if ( any(separator.ok == separator) != TRUE ) {
    mess <- paste("The provided separator \"", separator,
                  "\" is not acceptable. ",
                  "The default separator \"_\" is being used.", sep = "")
    separator <- "_"
    message(mess)
  }

  ##----- create the unique atom hashes
  uniq.atom.hashes <- as.vector(apply(atoms.oi[, cols.oi],
                                      1, paste, sep = "",
                                      collapse = separator) )

  ##----- remove blank spaces if not the separator
  if ( separator != " ") {
    uniq.atom.hashes <- gsub(x = uniq.atom.hashes,
                             pattern = " ",
                             replacement = "")
  }

  ##----- return the results
  return(uniq.atom.hashes)

}


## write basic PDB docs ---------------------------------------------
#' @title Write Basic PDB File
#' @description Writes standard PDB file.
#' @details Using the [bio3d::write.pdb()] function this function writes a PDB
#'   file from a [base::data.frame()] containing the typical PDB file
#'   information. This function is called from the [FreeSASA.diff()] function
#'   within the [HydrophilicityEvaluation()] function.
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param file Filename with ".pdb" extension.
#' @param atoms.oi The atoms [base::data.frame()].
#'
#' @return Writes a PDB file for the [FreeSASA.diff()] function.
#'
#' @examples
#'   \dontrun{
#'     write.basic.pdb(file = "just_some_PDB.pdb", atoms.oi)
#'   }
#'
#' @import bio3d
#'
#' @family utilities
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
write.basic.pdb <- function(file, atoms.oi) {

  bio3d::write.pdb(pdb = NULL,
                   file = file,
                   type = atoms.oi$type,
                   eleno = atoms.oi$eleno,
                   elety = atoms.oi$elety,
                   resid = atoms.oi$resid,
                   chain = atoms.oi$chain,
                   resno = atoms.oi$resno,
                   xyz = as.vector(unlist(t(atoms.oi[, c("x", "y", "z")]))),
                   o = atoms.oi$o,
                   b = atoms.oi$b,
                   elesy = atoms.oi$elesy
  )

}


## write conserved waters PDB docs ---------------------------------------------
#' @title Write Conserved Waters to PDB File
#' @description Writes conserved water information to a PDB file.
#' @details Using the [bio3d::write.pdb()] function this function writes a PDB
#'   file for the conserved water oxygen atoms with the percentage of structures
#'   with a water participting in the cluster (written to the occupancy column)
#'   and the calculated B-value -- using the rmsf of the waters in the cluster
#'   -- for the waters participating in the cluster (written to the B-value
#'   column). This function is called from the [ConservedWaters()] function.
#'
#'   All water molecules will include the water's oxygen atom (`elety`), be
#'   assigned the residue name (`resid`) `HOH`, and the chain (`chain`) `A`
#'   while the atom number (`eleno`) and residue number (`resno`) both start at
#'   `1`.
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param file Filename with ".pdb" extension.
#' @param h2o.clusters.summary The conserved water clusters summary.
#'
#' @return Writes a PDB file with the X, Y, and Z coordinates, percent conserved
#'   within the analyzed structures, and the calculated B-value for the oxygen
#'   atoms of the clustered waters.
#'
#' @examples
#'   \dontrun{
#'     write.conservedWaters.pdb(file = "system_conservedWaters.pdb",
#'                               h2o.clusters.summary)
#'   }
#'
#' @import bio3d
#'
#' @family utilities
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
write.conservedWaters.pdb <- function(file,
                                      h2o.clusters.summary) {

  num.clusters <- nrow(h2o.clusters.summary)
  bio3d::write.pdb(pdb = NULL,
                   file = file,
                   type = rep("HETATM", num.clusters),
                   eleno = 1:num.clusters,
                   elety = rep("O", num.clusters),
                   resid = rep("HOH", num.clusters),
                   chain = rep("A", num.clusters),
                   resno = 1:num.clusters,
                   xyz = as.vector(unlist(t(h2o.clusters.summary[, c("x", "y", "z")]))),
                   o = (h2o.clusters.summary$pct.conserved / 100),
                   b = h2o.clusters.summary$b.calc,
                   elesy = rep("O", num.clusters)
  )

}
