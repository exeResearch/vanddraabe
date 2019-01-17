## CleanProteinStructures.R
##
## oct-06-2016 (exe) created
## oct-24-2016 (exe) added D elesy (deuterated water hydrogen) to
##                   RemoveHydrogenAtoms
## oct-25-2016 (exe) identify duplicate PDB IDs and use filename (without
##                   extension) for PDB ID within analysis in the function
##                   CleanProteinStructures
## nov-26-2016 (exe) updated documentation
## feb-03-2017 (exe) change model atom occupancy value from 0.00 to 0.01
## feb-08-2017 (exe) added binning of occupancy and B-values
## feb-08-2017 (exe) added number of modeled atoms removed
## feb-08-2017 (exe) added number of waters removed and retained
## feb-08-2017 (exe) added B-values out of range (RemoveOoR.b)
## feb-08-2017 (exe) added occupancy out of range (RemoveOoR.o)
## feb-09-2017 (exe) added binning of B-values, normalized B-values,
##                   occupancy, and mobility
## apr-13-2017 (exe) updated documentation
## apr-14-2017 (exe) writes data to Excel workbook
## jul-25-2017 (exe) updated documentation
## jul-31-2017 (exe) updated RetainWatersWithinX() documentation
## aug-08-2017 (exe) added min.num.h2o parameter and check to CleanProteinStructures()
## aug-08-2017 (exe) CleanProteinStructures() now also returns the PDBids.retained
## aug-10-2017 (exe) added @importFrom stats ...
## jan-17-2019 (exe) expanded o.bins (-1,2), b.bins (-100,200), and mob.bins (-6,6) to
##                   accomodate fringe values
##
## Please direct all questions to Emilio Xavier Esposito, PhD
## exeResearch LLC, East Lansing, Michigan 48823 USA
## http://www.exeResearch.com
## emilio AT exeResearch DOT com
## emilio DOT esposito AT gmail DOT com
##
## Copyright (c) 2019, Emilio Xavier Esposito
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


## RemoveHydrogenAtoms docs ----------------------------------------------------
#' @title Remove Hydrogen and Deuterium Atoms
#' @description Removes hydrogen atoms from a RCSB/PDB structure.
#' @details Removes hydrogen and deuterium atoms from a PDB formatted
#'   [base::data.frame()] with PDB formatted information.
#'
#' @param atoms.chains.oi The data.frame containing the PDB file information;
#'   aka the PDB structure
#'
#' @return data.frame of the PDB structure _without_ hydrogen or deuterium
#'   atoms
#'
#' @examples
#'   \donttest{
#'   PDB.5rxn.noHydrogens <- RemoveHydrogenAtoms(PDB.5rxn$atom)
#'   }
#'
#' @export
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Clean Protein Structure"
#'
RemoveHydrogenAtoms <- function(atoms.chains.oi) {

  ##----- identify hydrogen and deuterium atoms
  heavy.atoms.tf <- ( (atoms.chains.oi$elesy != "H") &
                        (atoms.chains.oi$elesy != "D") )

  ##----- remove hydrogen and deuterium atoms
  atoms.chains.oi <- atoms.chains.oi[heavy.atoms.tf, ]

  ##----- return data.frame withOUT hydrogen and deuterium atoms
  return(atoms.chains.oi)

}


## RemoveModeledAtoms docs -----------------------------------------------------
#' @title Remove Modeled Atoms
#' @description Removes modeled atoms from a RCSB/PDB structure.
#' @details Sometimes atoms are not well resolved within the electron density
#'   maps and the scientists resolving/determining the structures "model back
#'   into" the resulting structure the atoms based on historical data. This is
#'   most common for residues where a portion of the residue is missing and
#'   based on the structure the missing atoms are replaces. These modeled atoms
#'   have an occupancy value of 0.01 or less and are identified and removed.
#'
#'   The reported occupancy value of 0.01 is used as the cutoff because several
#'   PDB structures have comments in the `REMARK 3` section stating,
#'   `"...MISSING ABOVE 1SIGMA WERE GIVEN A 0.01 OCCUPANCY..."`` or
#'   `"...WITH NO DENSITIES ARE GIVEN OCCUPANCY VALUES OF 0.01..."`.
#'
#' @param atoms.chains.oi The data.frame containing the PDB file information;
#'   aka the PDB structure
#'
#' @return data.frame of the PDB structure _without_ the modeled atoms
#'
#' @examples
#'   \donttest{
#'   PDB.1ecd.noModeledAtoms <- RemoveModeledAtoms(PDB.1ecd$atom)
#'   }
#'
#' @export
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Clean Protein Structure"
#'
RemoveModeledAtoms <- function(atoms.chains.oi) {

  ##----- identify modeled heavy atoms
  modeled.heavy.tf <- atoms.chains.oi[atoms.chains.oi$elesy != "H", "o"] <= 0.01
  if ( sum(modeled.heavy.tf) > 0 ) {
    ##-- retain only non-modeled atoms
    atoms.chains.oi <- atoms.chains.oi[!modeled.heavy.tf, ]
  }

  ##----- return data.frame withOUT modeled atoms
  return(atoms.chains.oi)

}


## Remove Occupancy Out of Range atoms docs ------------------------------------
#' @title Remove Occupancy Out of Range Atoms
#' @description Removes atoms with occupancy values out of accepted range.
#' @details Accepted occupancy values range from 0 to 1 with values for
#'   modeled atoms being 0.0 or 0.01 and highly conserved or represented atoms
#'   throughout the lattice having values greater than 0.9 and commonly
#'   possessing values of 1.0. This function identifies occupancy values
#'   less than 0 and greater than 1 and removes them from the structure.
#'
#' @param atoms.chains.oi The data.frame containing the PDB file information;
#'   aka the PDB structure
#'
#' @return data.frame of the PDB structure _without_ the offending atoms
#'
#' @examples
#'   \donttest{
#'   nrow(PDB.4dfr$atom)
#'   PDB.4dfr.OoR.o <- RemoveOoR.o(PDB.4dfr$atom)
#'   nrow(PDB.4dfr.OoR.o)
#'   }
#'
#' @export
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Clean Protein Structure"
#'
RemoveOoR.o <- function(atoms.chains.oi) {

  ##----- get occupancy values
  o.values <- atoms.chains.oi$o

  ##----- identify atoms with occupancy values In Range (0 - 1)
  o.values.ir.tf <- ((o.values >= 0) & (o.values <= 1))

  ##----- remove atoms with out of range values
  atoms.chains.oi <- atoms.chains.oi[o.values.ir.tf, ]

  ##----- return data.frame withOUT out of range atoms based on occupancy
  return(atoms.chains.oi)

}


## Remove B-value Out of Range atoms docs ------------------------------------
#' @title Remove B-value Out of Range Atoms
#' @description Removes atoms with B-values out of accepted range.
#' @details Accepted B-value values range from 0 to 100 with values. Atoms are
#'   considered stationary -- possessing low thermal energy -- when possessing
#'   values between 20 and 40 while larger values between 60 and 100 indicate a
#'   large amount of position variability within the lattice. This function
#'   identifies occupancy values less than 0 and greater than 100 and removes
#'   them from the structure.
#'
#' @param atoms.chains.oi The data.frame containing the PDB file information;
#'   aka the PDB structure
#'
#' @return data.frame of the PDB structure _without_ the offending atoms
#'
#' @examples
#'   \donttest{
#'   nrow(PDB.4ape$atom)
#'   PDB.4ape.OoR.b <- RemoveOoR.b(PDB.4ape$atom)
#'   nrow(PDB.4ape.OoR.b)
#'   }
#'
#' @export
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Clean Protein Structure"
#'
RemoveOoR.b <- function(atoms.chains.oi) {

  ##----- get occupancy values
  b.values <- atoms.chains.oi$b

  ##----- identify atoms with B-value values In Range (0 - 100)
  b.values.ir.tf <- ((b.values >= 0) & (b.values <= 100))

  ##----- remove atoms with out of range values
  atoms.chains.oi <- atoms.chains.oi[b.values.ir.tf, ]

  ##----- return data.frame withOUT out of range atoms based on B-value
  return(atoms.chains.oi)

}


## RetainWatersWithinX docs ----------------------------------------------------
#' @title Retain Waters Within X Angstroms of Protein
#' @description Retains water oxygen atoms within a user defined distance
#' @details Retain water oxygen atoms within a user defined distance. This
#'   function is a coarse grain method of removing waters beyond a predefined
#'   distance to reduce the computational load associated with the
#'   [stats::dist()] function for a collection of protein structure.
#'
#' @param atoms.dist Atomic distances calculated with the [stats::dist()]
#'   function
#' @param prot.het.h2o.idc List of protein, HET-atom, and water atom indices
#' @param cutoff.prot.h2o.dist User defined maximum numerical distance, in
#'   Angstroms, between the protein and water oxygen atoms to be retained.
#'
#' @return numerical vector of water oxygen atom indicies to retain
#'
#' @examples
#'   \donttest{
#'   ##--- determine the protein, hetatom, and  water indices
#'   prot.het.h2o.idc <- ProtHetWatIndices(data=PDB.1hah.aoi.clean)
#'
#'   ##--- calculate the distances
#'   atoms.dist <- as.matrix(dist(PDB.1hah.aoi.clean[, c("x","y","z")],
#'                                method="euclidean",
#'                                diag=TRUE, upper=TRUE))
#'   diag(atoms.dist) <- NA
#'
#'   water.idc.within.6 <- RetainWatersWithinX(atoms.dist,
#'                                             prot.het.h2o.idc,
#'                                             cutoff.prot.h2o.dist=6.0)
#'   # - 204 of the 204 water oxygen atoms are within 6 Angstroms of the protein
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Clean Protein Structure"
#'
RetainWatersWithinX <- function(atoms.dist,
                                prot.het.h2o.idc,
                                cutoff.prot.h2o.dist) {

  ##----- "unpack" the indicies
  ##--- protein indices
  prot.idc <- prot.het.h2o.idc$prot.idc
  ##--- water indices
  h2o.idc <- prot.het.h2o.idc$h2o.idc

  ##----- retain waters within user defined distance (Angstroms) of the protein
  prot.h2o.dists <- atoms.dist[prot.idc, h2o.idc]
  prot.h2o.dists.tf <- prot.h2o.dists <= (cutoff.prot.h2o.dist + 0.001)
  prot.h2o.nearby <- which(colSums(prot.h2o.dists.tf) >= 1)
  prot.h2o.nearby.idc <- h2o.idc[prot.h2o.nearby]
  h2o.mess <- paste(" -", length(prot.h2o.nearby.idc), "of the",
                    length(h2o.idc), "water oxygen atoms are within",
                    cutoff.prot.h2o.dist, "Angstroms of the protein",
                    sep=" ")
  message(h2o.mess)
  ##--- the water indices for waters within user defined
  ##    Angstroms of the protein
  h2o.idc <- h2o.idc[as.vector(prot.h2o.nearby)]

  ##----- return the water indices for waters within user defined
  ##      Angstroms of the protein
  return(h2o.idc)

}


## CleanProteinStructures docs -------------------------------------------------
#' @title Clean Protein Structures
#' @description Removes hydrogen and modeled atoms from a RCSB/PDB structure
#'   along with waters beyond a user defined distance from protein atoms.
#' @details PDB files obtained from the PDB conform to a specific set of
#'   formatting standards but this does not mean the data within the PDB files
#'   is always correct. This function _cleans_ the PDB file and summaries the
#'   atom evaluations.
#'
#'   This function does the following (in this order):
#'   * Reads in the PDB file
#'   * Adds/updates the element symbol (`elesy`) using the atom type (`elety`)
#'     via the [bio3d::atom2ele()] function
#'   * Removes hydrogen atoms via [RemoveHydrogenAtoms()] (user option)
#'   * Removes atoms with occupancy values determined to be out of range (OoR)
#'     via [RemoveOoR.o()]
#'   * Removes atoms with B-values determined to be out of range (OoR)
#'     via [RemoveOoR.b()]
#'   * Bins (counts) the occupancy values
#'   * Bins (counts) the B-values
#'   * Bins (counts) the normalized B-values
#'   * Bins (counts) the mobility values
#'   * Removes modeled atoms via [RemoveModeledAtoms()] (user option)
#'   * Removes water oxygen atoms greater than user defined value
#'     `cutoff.prot.h2o.dist` from the protein via [RetainWatersWithinX()]
#'     (user option)
#'   * Writes cleaned protein structure to a PDB file
#'
#' @param prefix The directory with the PDB files to be cleaned
#' @param CleanHydrogenAtoms A logical indication if hydrogen atoms should be
#'   removed; default: `TRUE`
#' @param CleanModeledAtoms A logical indication if modeled atoms should be
#'   removed; default: `TRUE`
#' @param cutoff.prot.h2o.dist A numerical value setting the maximum distance
#'   between a protein atom (heteroatoms are ignored) and water oxygen atoms.
#'   The oxygen atoms equal to or less than this distance are retained;
#'   default: `6.0` Angstroms
#' @param min.num.h2o Minimum number of water oxygen atoms within a protein
#'   structure for it to be included in the conserved water analysis;
#'   default: 20
#' @param cleanDir A character string for the "cleaned" PDB structures to be
#'   written. The provided character string are appended with "_CLEANED";
#'   default: `"ProteinSystem"`
#' @param filename The filename prefix for the returned results. Default is
#'   "ProteinSystem"
#'
#' @return
#'   The following data is returned:
#'   * **cleaning.summary**: summary indicating
#'     - if hydrogen atoms were removed `TRUE/FALSE`
#'     - number of out of range atoms for B-values and occupancy values
#'     - number of modeled (and thus removed)
#'     - number of atoms _**NOT**_ modeled (and thus retained)
#'     - number of water oxygen atoms beyond the user defined cutoff
#'     - the number of water oxygen atoms within the user defined cutoff.
#'   * **Bvalue.counts**: binned B-value values with binwidths = 5 (0 to 100)
#'   * **normBvalue.counts**: binned normalized B-value values with binwidths
#'     = 0.1 (-4 to 6)
#'   * **occupancy.counts**: binned occupancy values with binwidths = 0.1 (0 to
#'     1)
#'   * **mobility.counts**: binned mobility values with binwidths = 0.1 (0 to 6)
#'   * **Excel workbook**: containing the cleaning.summary, Bvalue.counts,
#'     normBvalue.counts, occupancy.counts, and mobility.counts
#'     data as individual tabs
#'   * **PDBids.retained**: a vector of PDBids
#'   * **call**: parameters provided by the user
#'
#' @export
#'
#' @import bio3d
#' @importFrom stats dist
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Clean Protein Structure"
#'
CleanProteinStructures <- function(prefix="./alignTesting",
                                   CleanHydrogenAtoms=TRUE,
                                   CleanModeledAtoms=TRUE,
                                   cutoff.prot.h2o.dist=6.0,
                                   min.num.h2o = 20,
                                   cleanDir="ProteinSystem",
                                   filename="ProteinSystem") {


  ##----- the provided request
  the.call <- match.call()

  ##----- check the user provided cutoff values
  if ( is.numeric(cutoff.prot.h2o.dist) ) {
    RetainAtoms <- TRUE
  } else {
    RetainAtoms <- FALSE
  }

  ## CREATE now.date.time AND WORKSHEET NAMES ----------------------------------
  ##----- create the date and time portion of the filenames
  current.time <- Sys.time()
  now.date.time <- FileTimeStamp(current.time)

  ##----- construct Excel sheetnames
  PDBs.cleaning.summary <- paste("PDB_cleanSumm", now.date.time, sep="_")
  PDBs.Bvalue.counts    <- paste("BvalueBins", now.date.time, sep="_")
  PDBs.nBvalue.counts   <- paste("nBvalueBins", now.date.time, sep="_")
  PDBs.occ.counts       <- paste("OccBins", now.date.time, sep="_")
  PDBs.mob.counts       <- paste("MobilityBins", now.date.time, sep="_")


  ## GET THE PDB LOCATION INFORMATION ------------------------------------------
  ##----- get list of PDB files within prefix
  pdb.location <- ReturnPDBfullPath(prefix)

  ##----- determine PDB IDs
  pdb.ids <- ExtractPDBids(pdb.location)
  ##--- look for duplicate names and use filenames if duplicates occur
  duplicates.tf <- sum(duplicated(pdb.ids)) > 0
  if ( duplicates.tf == TRUE ) {
    pdb.files <- basename(pdb.location)
    pdb.ids <- gsub(pattern = ".pdb",
                    replacement = "",
                    x = pdb.files,
                    fixed = TRUE)
  }

  ##----- determine number of structures
  num.pdb.structures <- length(pdb.ids)


  ## CONSTRUCT THE DIRECTORY FOR CLEANED STRUCTURES ---------------------------
  os.type <- .Platform$OS.type
  dirName <- dirname(cleanDir)
  baseName <- basename(cleanDir)
  dir.clean <- paste0(dirName, "/", baseName, "_CLEANED")
  if ( os.type == "windows" ) {
    dir.clean <- gsub(pattern = "/", replacement = "\\",
                      x = dir.clean, fixed = TRUE)
  }
  dir.create(dir.clean, showWarnings = FALSE)


  ## SETUP BINNING FOR OCCUPANCY -----------------------------------------------
  o.bins <- c(-1, seq(from = 0, to = 1, by = 0.1), 2)
  # o.bins.colNames <- paste("o_", round(o.bins, digits = 1), sep = "")
  o.bins.results <- data.frame(
    matrix(rep(0, num.pdb.structures * length(o.bins) ),
           nrow=num.pdb.structures), stringsAsFactors = FALSE)
  colnames(o.bins.results) <- round(o.bins, digits = 1)


  ## SETUP BINNING FOR B-VALUES -----------------------------------------------
  b.bins <- c(-100, seq(from = 0, to = 100, by = 5), 200)
  # b.bins.colNames <- paste("B_", b.bins, sep = "")
  b.bins.results <- data.frame(
    matrix(rep(0, num.pdb.structures * length(b.bins) ),
           nrow=num.pdb.structures), stringsAsFactors = FALSE)
  colnames(b.bins.results) <- b.bins


  ## SETUP BINNING FOR MOBILITY ------------------------------------------------
  mob.bins <- c(-6, seq(from = 0, to = 6, by = 0.1))
  # mob.bins.colNames <- paste("mob_", round(mob.bins, digits = 1), sep = "")
  mob.bins.results <- data.frame(
    matrix(rep(0, num.pdb.structures * length(mob.bins) ),
           nrow=num.pdb.structures), stringsAsFactors = FALSE)
  colnames(mob.bins.results) <- round(mob.bins, digits = 1)


  ## SETUP BINNING FOR NORMALIZED B-VALUES ------------------------------------
  normB.bins <- seq(from = -7, to = 7, by = 0.1)
  # normB.bins.colNames <- paste("normB_", round(normB.bins, digits = 1),
  #                              sep = "")
  normB.bins.results <- data.frame(
    matrix(rep(0, num.pdb.structures * length(normB.bins) ),
           nrow=num.pdb.structures), stringsAsFactors = FALSE)
  colnames(normB.bins.results) <- round(normB.bins, digits = 1)


  ## CREATE DATA.FRAME FOR REPORTING MODIFICATIONS ----------------------------
  false.vector <- rep(FALSE, num.pdb.structures)
  zero.vector  <- rep(0, num.pdb.structures)
  df.cleaning <- data.frame(row.names = pdb.ids,
                            removedHydrogens = false.vector,
                            num.o.OoR = zero.vector,
                            num.b.OoR = zero.vector,
                            num.Modeled = zero.vector,
                            num.notModeled = zero.vector,
                            num.WatersDistantRemoved = zero.vector,
                            num.WatersRetained = zero.vector,
                            stringsAsFactors = FALSE)


  ## CLEAN THE STRUCTURES ------------------------------------------------------
  for (curr.pdb.idx in 1:num.pdb.structures) {

    ## READ IN PDB -------------------------------------------------------------
    curr.pdb.id <- pdb.ids[curr.pdb.idx]
    message(paste("Cleaning ", curr.pdb.id, "...", sep = ""))
    pdb <- bio3d::read.pdb2(file=pdb.location[curr.pdb.idx])
    atoms.oi <- pdb$atom
    num.atoms <- nrow(atoms.oi)


    ## ADD/UPDATE ELEMENTY SYMBOL (elesy) USING ATOM TYPE (elety) --------------
    suppressWarnings( atoms.oi$elesy <- bio3d::atom2ele(atoms.oi$elety) )


    ## REMOVE HYDROGEN ATOMS ---------------------------------------------------
    if ( CleanHydrogenAtoms == TRUE ) {
      atoms.oi <- RemoveHydrogenAtoms(atoms.chains.oi = atoms.oi)
      num.atoms.noHydrogen <- nrow(atoms.oi)

      ##--- information to the user
      if ( num.atoms > num.atoms.noHydrogen) {
        message(" - Removed hydrogen atoms")
        num.atoms <- nrow(atoms.oi)
        df.cleaning[curr.pdb.idx, "removedHydrogens"] <- TRUE
      }
    }

    ## CHECK FOR WATERS. IF NOT, NEXT STRUCTURE -----
    atoms.oi.resid <- atoms.oi$resid
    ##--- check to see if structure has at least the user defined number of
    ##    waters. based on Kuhn et al. hydrophilicity article it should be
    ##    20 waters...
    has.h2o.tf <- HasXWaters(atoms.oi.resid, min.num.h2o = min.num.h2o)
    if (has.h2o.tf$has.h2o.tf == FALSE) {
      mess <- paste(" <<<|||>>> NOTE: ", curr.pdb.id, " contains ",
                    has.h2o.tf$num.water, " waters! This is less than the ",
                    "minimum of ", min.num.h2o,
                    " required waters for the hydrophilicity evaluation. ",
                    "Moving to the next structure <<<|||>>>.", sep = "")
      message(mess)
      next()
    }

    ## REMOVE ATOMS WITH OCCUPANCY VALUES OUT OF THE 0-1 RANGE -----------------
    num.atoms <- nrow(atoms.oi)
    atoms.oi <- RemoveOoR.o(atoms.chains.oi = atoms.oi)
    df.cleaning[curr.pdb.idx, "num.o.OoR"] <- num.atoms - nrow(atoms.oi)


    ## REMOVE ATOMS WITH B-VALUES OUT OF THE 0-100 RANGE ----------------
    num.atoms <- nrow(atoms.oi)
    atoms.oi <- RemoveOoR.b(atoms.chains.oi = atoms.oi)
    df.cleaning[curr.pdb.idx, "num.b.OoR"] <- num.atoms - nrow(atoms.oi)


    ## BIN OCCUPANCY VALUES ----------------------------------------------------
    o.values <- atoms.oi$o
    o.counts <- data.frame( table( findInterval(o.values, vec = o.bins) ) )
    o.bins.results[curr.pdb.idx,
                   as.numeric( as.character(o.counts$Var1) )] <- o.counts$Freq


    ## BIN B-VALUES -----------------------------------------------------
    b.values <- atoms.oi$b
    b.counts <- data.frame( table( findInterval(b.values, vec = b.bins) ) )
    b.bins.results[curr.pdb.idx,
                   as.numeric( as.character(b.counts$Var1) )] <- b.counts$Freq


    ## BIN NORMALIZED B-VALUES ------------------------------------------
    nBvalues <- NormalizedBvalue(Bvalues = b.values)
    normB.counts <- data.frame( table( findInterval(nBvalues,
                                                    vec = normB.bins) ) )
    normB.bins.results[curr.pdb.idx,
                       as.numeric( as.character(normB.counts$Var1) )] <-
      normB.counts$Freq


    ## BIN MOBILITY VALUES -----------------------------------------------------
    mobility <- Mobility(Bvalues = b.values, occupancy = o.values)
    mob.counts <- data.frame( table( findInterval(mobility, vec = mob.bins) ) )
    mob.bins.results[curr.pdb.idx,
                     as.numeric( as.character(mob.counts$Var1) )] <-
      mob.counts$Freq


    ## REMOVE MODELED ATOMS ----------------------------------------------------
    if ( CleanModeledAtoms == TRUE ){
      num.atoms <- nrow(atoms.oi)
      atoms.oi <- RemoveModeledAtoms(atoms.chains.oi = atoms.oi)
      num.atoms.notModeled <- nrow(atoms.oi)
      df.cleaning[curr.pdb.idx, c("num.Modeled", "num.notModeled")] <-
        c(num.atoms - num.atoms.notModeled, num.atoms.notModeled)

      ##--- information to the user
      if ( num.atoms > num.atoms.notModeled) {
        message(" - Removed modeled atoms")
      }
    }


    ## REMOVE DISTANT WATER MOLECULES/OXYGEN ATOMS -----------------------------
    num.atoms <- nrow(atoms.oi)
    ##----- determine the protein, hetatom, and  water indices
    prot.het.h2o.idc <- ProtHetWatIndices(data=atoms.oi)
    prot.idc <- prot.het.h2o.idc$prot.idc
    het.idc <-  prot.het.h2o.idc$het.idc
    h2o.idc.orig <-  prot.het.h2o.idc$h2o.idc
    num.h2o.idc.orig <- length(h2o.idc.orig)
    df.cleaning[curr.pdb.idx, "num.WatersRetained"] <- num.h2o.idc.orig

    if ( RetainAtoms == TRUE ) {
      ##--- calculate the distances
      atoms.dist <- as.matrix(dist(atoms.oi[, c("x","y","z")],
                                   method = "euclidean",
                                   diag = TRUE, upper = TRUE))
      diag(atoms.dist) <- NA

      ##--- determine waters to retain
      h2o.idc <- RetainWatersWithinX(atoms.dist,
                                     prot.het.h2o.idc,
                                     cutoff.prot.h2o.dist)

      ##--- protein, hetatm, and retained waters
      atoms.oi <- atoms.oi[c(prot.idc, het.idc, h2o.idc), ]

      ##--- update df.cleaning
      df.cleaning[curr.pdb.idx,
                  c("num.WatersDistantRemoved", "num.WatersRetained")] <-
        c(num.h2o.idc.orig - length(h2o.idc), length(h2o.idc))

    }
    num.atoms <- nrow(atoms.oi)

    ## WRITE CLEANED PDB TO FILE -----------------------------------------------
    filename.cleaned <- paste0(curr.pdb.id, "_cleaned.pdb")
    path.filename.cleaned <- paste(dir.clean,
                                   filename.cleaned,
                                   sep = "/")
    mess <- paste(" - Wrote", filename.cleaned,
                  "to", dir.clean, "\n",
                  sep = " ")
    message(mess)
    write.basic.pdb(file = path.filename.cleaned, atoms.oi)

  }


  ## ADD PDBids TO b, normB, o, and mobility .bins.results ---------------------
  pdb.ids <- rownames(df.cleaning)
  rownames(b.bins.results) <- rownames(o.bins.results) <-
    rownames(mob.bins.results) <- rownames(normB.bins.results) <- pdb.ids

  b.bins.results <- data.frame(PDBid = pdb.ids, b.bins.results,
                               stringsAsFactors = FALSE)
  normB.bins.results <- data.frame(PDBid = pdb.ids, normB.bins.results,
                                   stringsAsFactors = FALSE)
  o.bins.results <- data.frame(PDBid = pdb.ids, o.bins.results,
                               stringsAsFactors = FALSE)
  mob.bins.results <- data.frame(PDBid = pdb.ids, mob.bins.results,
                                 stringsAsFactors = FALSE)


  ## WRITE RESULTS TO EXCEL WORKBOOK -------------------------------------------
  ##----- construct workbook name
  filename.xlsx <- paste(filename, "_DATA_RESULTS.xlsx", sep="")
  ##--- open existing excel workbook if available
  if ( file.exists(filename.xlsx) == TRUE ) {
    results.wb <- openxlsx::loadWorkbook(file=filename.xlsx)
  } else {  ##--- construct the workbook
    results.wb <- openxlsx::createWorkbook()
  }

  ##--- cleaned structure summary information
  ox.df.cleaning <- cbind(PDBids=rownames(df.cleaning),
                          df.cleaning)
  results.wb <- oxPDBcleanedSummarySheet(wb.name = results.wb,
                                         sheet.name = PDBs.cleaning.summary,
                                         df = ox.df.cleaning)
  ##--- B-values counts data
  results.wb <- oxPlainDataSheet(wb.name = results.wb,
                                 sheet.name = PDBs.Bvalue.counts,
                                 df = b.bins.results)
  ##--- normalized B-values counts data
  results.wb <- oxPlainDataSheet(wb.name = results.wb,
                                 sheet.name = PDBs.nBvalue.counts,
                                 df = normB.bins.results)
  ##--- occupancy values counts data
  results.wb <- oxPlainDataSheet(wb.name = results.wb,
                                 sheet.name = PDBs.occ.counts,
                                 df = o.bins.results)
  ##--- mobility values counts data
  results.wb <- oxPlainDataSheet(wb.name = results.wb,
                                 sheet.name = PDBs.mob.counts,
                                 df = mob.bins.results)

  ##--- write the workbook
  openxlsx::saveWorkbook(results.wb, filename.xlsx, overwrite=TRUE)

  message("----- Results written to Excel workbook _____\n")


  ## DETERMINE RETAINED (CLEAN) PDBids -----------------------------------------
  PDBids.retained.tf <- df.cleaning$num.WatersRetained > 0
  PDBids.retained <- rownames(df.cleaning)[PDBids.retained.tf]


  ## RETURN THE RESULTS --------------------------------------------------------
  list(cleaning.summary = df.cleaning,
       Bvalue.counts = b.bins.results,
       normBvalue.counts = normB.bins.results,
       occupancy.counts = o.bins.results,
       mobility.counts = mob.bins.results,
       PDBids.retained = PDBids.retained,
       call = the.call)

}
