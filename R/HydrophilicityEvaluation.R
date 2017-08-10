## HydrophilicityEvaluation.R
##
## sep-16-2016 (exe) created
## feb-23-2017 (exe) added getResidueData function
## feb-24-2017 (exe) added getProtAtomsNearWater function
## feb-27-2017 (exe) added check for existing directory(ies)
## feb-27-2017 (exe) added check for more than zero structures
## feb-28-2017 (exe) remove hydrogen atoms within HydrophilicityEvaluation
## mar-02-2017 (exe) added atom count information to PDB.info output
## mar-02-2017 (exe) removed warnings() from the hydrogen atom check for removal
## mar-06-2017 (exe) added calcAtomClassHydrophilicity()
## mar-06-2017 (exe) added calcAtomHydrationEstimate()
## jul-25-2017 (exe) updated documentation
## jul-31-2017 (exe) corrected \dontrun in HydrophilicityEvaluation()'s example
## jul-31-2017 (exe) updated HydrophilicityEvaluation() and
##                   calcAtomHydrationEstimate() documentation
## aug-08-2017 (exe) updated HydrophilicityEvaluation()'s check for waters message
## aug-09-2017 (exe) added @importFrom stats ...
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


## HydrophilicityEvaluation docs -----------------------------------------------
#' @title Hydrophilicity Evaluation
#' @description Calculate the hydrophilicity values for a set of protein
#'   structures.
#' @details The hydrophilicity values of individual atomtypes is determined
#'   using a collection of protein structures. For each water oxygen atom within
#'   at the most 4 Angstroms of a solvent accessible (exposed) protein atom,
#'   these occurrences are recorded. The number of solvent accessible atom types
#'   interacting with a water molecule are divided by the number of solvent
#'   accessible atom types. In general the more diverse data available, the
#'   better the informatics based hydrophilicity values should correlate with
#'   various experimental values.
#'
#'   _**NOTE**_: Hydrogen atoms are removed for instances when the protein
#'   structures have not be cleaned with [CleanProteinStructures()].
#'
#' @param prefix The directory containing the protein structures; _e.g._,
#'   "alignTesting/"
#' @param h2o.prot.dist.max Maximum distance between the water oxygen atoms and
#'   the protein for consideration in the determination for hydrophilicity
#'   values; default: 6.0
#' @param bound.h2o.dist.max Maximum distance between the water oxygen atoms and
#'   the protein for inclusion in the calculation of hydrophilicity values;
#'   default: 4.0
#' @param min.num.h2o Minimum number of water oxygen atoms within a protein
#'   structure for it to be included in the calculation of hydrophilicity
#'   values; default: 20
#' @param probeRadius Water molecule probe radius; default: 1.4
#' @param dataset Name of the dataset to be used; _e.g._,"top56"
#'
#' @return This function returns:
#'   * **PDB.info**: a summary of the data for each protein structure analyzed
#'     - _PDBid_: PDB id
#'     - _time_: duration for hydrophilicity evaluation
#'     - _num.res_: number of protein residues
#'     - _num.res.buried_: number of protein residues with _**NO**_ solvent
#'       exposure
#'     - _num.res.SurExp_: number of protein residues with solvent accessible
#'       surface area
#'     - _pct.res.SurExp_: percentage of protein residues with solvent
#'     - _SASA.total_: total protein solvent accessible surface area;
#'       Angstroms^2^
#'     - _SASA.lost_: total protein solvent accessible surface area lost due to
#'       bound waters; Angstroms^2^
#'     - _pct.SASA.exposed_: percentage protein solvent accessible surface area
#'       \eqn{(SASA.total - SASA.lost) / SASA.total}
#'     - _num.prot.atom_: number of protein atoms
#'     - _num.atom.buried_: number of protein atoms with _**NO**_ solvent
#'       exposure
#'     - _num.atom.SurExp_: number of protein atoms with solvent accessible
#'       surface area
#'     - _pct.atom.SurExp_: percentage protein atoms with solvent accessible
#'       surface area \eqn{(SASA.total - SASA.lost) / SASA.total}
#'     - _num.h2o_: number of waters in the system
#'     - _num.h2o.lte.prot.max_: number of waters within `h2o.prot.dist.max`
#'       cutoff
#'     - _num.SurBound.h2o_: number of surface bound waters; water within
#'       `bound.h2o.dist.max` cutoff
#'     - _num.bb.h2o.inter_: number of backbone - water interactions
#'     - _num.sc.h2o.inter_: number of sidechain - water interactions
#'     - _num.res.h2o.inter_: number of interactions between residues and water
#'     - _num.h2o.res.inter_: number of interactions between water and residue
#'       (residues are a unit)
#'     - _num.h2o.resAtom.inter_: number of water-atom interactions
#'   * **SASA.results**: `data.frame` of protein atoms within the
#'     `h2o.prot.dist.max` of each water oxygen atom
#'   * **df.AtomTypes.all**: total number of AtomTypes for each structure
#'   * **df.AtomTypes.buried**: number of buried AtomTypes for each structure
#'   * **df.AtomTypes.SurExp**: number of surface exposed AtomTypes for each
#'       structure
#'   * **df.AtomTypes.h2o.nearby**: number of surface exposed AtomTypes within
#'       h2o.prot.dist.max (default 6 Ang) of an individual water
#'   * **df.AtomTypes.h2o.bound**: number of surface exposed AtomTypes within
#'       bound.h2o.dist.max (default 4 Ang) of an indvidual water
#'   * **df.AtomTypes.h2o.inter**: number of surface exposed AtomTypes with the
#'       shortest distance to an individual water
#'   * **df.residue.hydro**:
#'   * **HydrophilicityTable**: hydrophilicity table based on provided protein
#'       structures
#'   * **AtomTypeClasses.hydratFract**:
#'   * **no.h2o**: proteins (PDB IDs) without the _minimum_ number of user
#'       defined waters `min.num.h2o`
#'   * **call**: parameters provided by the user
#'   * **duration**: duration of complete [HydrophilicityEvaluation()]
#'       calculation
#'
#' @export
#'
#' @import bio3d
#' @importFrom stats dist
#'
#' @examples
#' \dontrun{
#'  HydrophilicityEvaluation <- function(prefix = "alignTesting/",
#'                                       h2o.prot.dist.max = 6.0,
#'                                       bound.h2o.dist.max = 4.0,
#'                                       min.num.h2o = 20,
#'                                       probeRadius = 1.4,
#'                                       dataset = "top56")
#' }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Hydrophilicity Evaluation" "Bound Water Environment"
#'
#' @references
#'   Leslie A Kuhn, Craig A Swanson, Michael E Pique, John A Tainer,
#'   and Elizabeth D Getzof. Atomic and Residue Hydrophilicity in the Context of
#'   Folded Protein Structures. _PROTEINS: Structure, Function, and
#'   Genetics_, 1995, **23** (_4_), pp 536-547.
#'   [DOI: 10.1002/prot.340230408](http://doi.org/10.1002/prot.340230408)
#'   [PMID: 8749849](http://www.ncbi.nlm.nih.gov/pubmed/8749849)
#'
HydrophilicityEvaluation <- function(prefix = "alignTesting/",
                                     h2o.prot.dist.max = 6.0,
                                     bound.h2o.dist.max = 4.0,
                                     min.num.h2o = 20,
                                     probeRadius = 1.4,
                                     dataset = "top56") {

  ##----- the starting time
  time.start.hydrophilicity.evaluation <- Sys.time()


  ## SETUP AND FILES -----------------------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##----- what the user entered
  the.call <- match.call()


  ##----- check if the directory(ies) exist
  prefix <- prefix[!duplicated(prefix)]
  num.dirs <- length(prefix)
  dirExists.tf <- dir.exists(prefix)
  if ( num.dirs == sum(!dirExists.tf) ) {
    stop("None of the provided directory(ies) in the prefix parameter exist. \n
         Exiting Hydrophilicity Evaluation.")
  }
  if ( num.dirs != sum(dirExists.tf) ) {
    prefix.missing <- prefix[!dirExists.tf]
    prefix.missing.list <- paste("  - ", prefix.missing, "\n", collapse = NULL)
    message("The following directory(ies) do NOT exist:")
    message(prefix.missing.list)

    ##--- update prefix with existing directory(ies)
    prefix <- prefix[dirExists.tf]
  }


  ##----- get list of PDB files within prefix
  pdb.location <- unlist(lapply(prefix, ReturnPDBfullPath))


  ##----- determine PDB IDs
  # pdb.ids <- ExtractPDBids(pdb.location)
  pdb.files <- basename(pdb.location)
  pdb.ids <- gsub(pattern = "_cleaned.pdb",
                  replacement = "",
                  x = pdb.files,
                  fixed = TRUE)


  ##----- determine number of structures
  num.pdb.structures <- length(pdb.ids)
  ##--- stop evaluation if no structures present
  if ( num.pdb.structures == 0 ) {
    stop("No protein structures provided. \n
         Exiting Hydrophilicity Evaluation.")
  }



  ## CREATE EMPTY DATA.FRAMES FOR DATA/RESULTS ---------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## data set information ------------------------------------------------------
  df.dataset <- data.frame(PDBid = pdb.ids,             ## column #1
                           time = NA,                   ## column #2
                           num.res = NA,                ## column #3
                           num.res.buried = NA,         ## column #4
                           num.res.SurExp = NA,         ## column #5
                           pct.res.SurExp = NA,         ## column #6
                           SASA.total = NA,             ## column #7
                           SASA.lost = NA,              ## column #8
                           pct.SASA.exposed = NA,       ## column #9
                           num.prot.atom = NA,          ## column #10
                           num.atom.buried = NA,        ## column #11
                           num.atom.SurExp = NA,        ## column #12
                           pct.atom.SurExp = NA,        ## column #13
                           num.h2o = NA,                ## column #14 num h2o in system
                           num.h2o.lte.prot.max = NA,   ## column #15 | num h2o within h2o.prot.dist.max
                           num.SurBound.h2o = NA,       ## column #16 | num  (surface bound waters) h2o within bound.h2o.dist.max
                           num.bb.h2o.inter = NA,       ## column #17
                           num.sc.h2o.inter = NA,       ## column #18
                           num.res.h2o.inter = NA,      ## column #19 | num interactions between residues and h2o
                           num.h2o.res.inter = NA,      ## column #20 | num interactions between h2o and residue (residues are a unit)
                           num.h2o.resAtom.inter = NA,  ## column #21 | num interactions between h2o and individual residue atoms
                           stringsAsFactors = FALSE)
  ##--- update column names to reflect user provided parameters
  colnames(df.dataset)[15] <- paste("num.h2o.lte.", h2o.prot.dist.max,
                                    "ang.prot",
                                    sep="")
  colnames(df.dataset)[16] <- paste("num.SurBound.h2o.lte.", bound.h2o.dist.max,
                                    "ang.prot",
                                    sep="")


  ## per residue hydrophilicity ------------------------------------------------
  df.residue.hydro <- as.data.frame(matrix(data = 0, nrow = 20, ncol = 7))
  rownames(df.residue.hydro) <- names.residues
  colnames(df.residue.hydro) <- c("num.res.tot",
                                  "num.res.buried",
                                  "num.res.SurExp",
                                  "num.bb.h2o",
                                  "num.sc.h2o",
                                  "num.res.h2o",
                                  "num.SurBound.h2o")
  df.residue.hydro.FINAL <- df.residue.hydro.ZEROs <- df.residue.hydro


  ## AtomTypes counts ----------------------------------------------------------
  num.names.res.AtomTypes <- length(names.res.AtomTypes)

  df.AtomTypes.all <-             ## total number of AtomTypes
    df.AtomTypes.buried <-        ## number of buried AtomTypes
    df.AtomTypes.SurExp <-        ## number of surface exposed AtomTypes
    df.AtomTypes.h2o.nearby <-    ## number of surface exposed AtomTypes within h2o.prot.dist.max (default 6 Ang) of an individual water
    df.AtomTypes.h2o.bound <-     ## number of surface exposed AtomTypes within bound.h2o.dist.max (default 4 Ang) of an indvidual water
    df.AtomTypes.h2o.inter <-     ## number of surface exposed AtomTypes with the shortest distance to an individual water
    as.data.frame(matrix(data = 0,
                         nrow = num.pdb.structures,
                         ncol = num.names.res.AtomTypes))

  ##--- rename rows
  rownames(df.AtomTypes.all) <- rownames(df.AtomTypes.buried) <-
    rownames(df.AtomTypes.SurExp) <- rownames(df.AtomTypes.h2o.nearby) <-
    rownames(df.AtomTypes.h2o.bound) <- rownames(df.AtomTypes.h2o.inter) <- pdb.ids
  ##--- rename columns
  colnames(df.AtomTypes.all) <- colnames(df.AtomTypes.buried) <-
    colnames(df.AtomTypes.SurExp) <- colnames(df.AtomTypes.h2o.nearby) <-
    colnames(df.AtomTypes.h2o.bound) <- colnames(df.AtomTypes.h2o.inter) <- names.res.AtomTypes


  ## ResTypes counts -----------------------------------------------------------
  ## ResTypes are the 20 naturally occurring amino acids
  num.names.residues <- length(names.residues)

  df.ResTypes.all <-             ## total number of each residue type (ResTypes)
    df.ResTypes.buried <-        ## number of buried residue types
    df.ResTypes.SurExp <-        ## number of surface exposed residue types
    df.ResTypes.sc.h2o <-        ## number of surface exposed side chains interacting with a water for each residue type
    df.ResTypes.bb.h2o <-        ## number of surface exposed backbones interacting with a water for each residue type
    df.ResTypes.res.h2o <-       ## number of surface exposed residues interacting with a water for each residue type
    df.ResTypes.SurBound.h2o <-  ## number of surface bound waters within bound.h2o.dist.max of each residue type
    as.data.frame(matrix(data = 0,
                         nrow = num.pdb.structures,
                         ncol = num.names.residues))

  ##--- rename rows
  rownames(df.ResTypes.all) <- rownames(df.ResTypes.buried) <-
    rownames(df.ResTypes.SurExp) <- rownames(df.ResTypes.sc.h2o) <-
    rownames(df.ResTypes.bb.h2o) <- rownames(df.ResTypes.res.h2o) <-
    rownames(df.ResTypes.SurBound.h2o) <- pdb.ids
  ##--- rename columns
  colnames(df.ResTypes.all) <- colnames(df.ResTypes.buried) <-
    colnames(df.ResTypes.SurExp) <- colnames(df.ResTypes.sc.h2o) <-
    colnames(df.ResTypes.bb.h2o) <- colnames(df.ResTypes.res.h2o) <-
    colnames(df.ResTypes.SurBound.h2o) <- names.residues


  ## create empty data.frames for results/data ---------------------------------
  df.results <- as.data.frame(matrix(data = NA,
                                     nrow = num.pdb.structures * 7500,
                                     ncol = 31),
                              stringsAsFactors = FALSE)


  ## START EVALUATING THE PROTEIN STRUCTURES -----------------------------------
  entry.start <- 1
  for (curr.pdb.idx in 1:num.pdb.structures) {

    time.start.curr.pdb.idx <- Sys.time()

    ##----- read in PDB
    curr.pdb.id <- pdb.ids[curr.pdb.idx]
    pdb <- bio3d::read.pdb2(file = pdb.location[curr.pdb.idx])
    atoms.oi <- pdb$atom
    ##--- check for OXT oxygens and correct
    atomType.oxt.tf <- atoms.oi$elety == "OXT"
    read.mess <- paste("Reading structure ", curr.pdb.id, " (index: ",
                       curr.pdb.idx, ")...", sep = "")
    if ( any(atomType.oxt.tf) ) {
      atoms.oi$elety[atoms.oi$elety == "OXT"] <- "O"
      read.mess <- paste(read.mess,
                         " FYI: Converting ", sum(atomType.oxt.tf),
                         " OXT C-terminus oxygen atom types to O",
                         sep = "")
    }
    message(read.mess)


    ##----- remove hydrogen atoms
    if ( any(atoms.oi$elesy == "H") | any(atoms.oi$elesy == "D") ) {
      atoms.oi <- RemoveHydrogenAtoms(atoms.chains.oi = atoms.oi)
    }
    # if ( any(atoms.oi$elesy == "D") ) {
    #   atoms.oi <- RemoveHydrogenAtoms(atoms.chains.oi = atoms.oi)
    # }


    ##----- check for waters
    atoms.oi.resid <- atoms.oi$resid
    ##--- check to see if structure has at least 20 waters (minimum number
    ##    based on Kuhn et al. hydrophilicity article)
    has.h2o.tf <- HasXWaters(atoms.oi.resid, min.num.h2o = min.num.h2o)
    df.dataset[curr.pdb.idx, c(14)] <- has.h2o.tf$num.water
    if (has.h2o.tf$has.h2o.tf == FALSE) {
      mess <- paste(" <<<|||>>> NOTE: ", curr.pdb.id, " contains ",
                    has.h2o.tf$num.water, " waters! This is less than the ",
                    "minimum of ", min.num.h2o,
                    " required waters for the hydrophilicity evaluation. ",
                    "Moving to the next structure <<<|||>>>.", sep = "")
      message(mess)
      next()
    }


    ##----- find and convert various charged form names to standard names
    atoms.oi$resid <- aaStandardizeNames(atoms.oi$resid)


    ##----- add unique protein atom hashes
    res.atom.ids <- UniqueAtomHashes(atoms.oi = atoms.oi,
                                     cols.oi = c("resid", "elety"),
                                     separator = " ")

    res.chain.resno.ids <- UniqueAtomHashes(atoms.oi,
                                      cols.oi = c("resid", "chain", "resno"),
                                      separator = "_")

    uniq.atom.ids <- UniqueAtomHashes(atoms.oi = atoms.oi,
                                      cols.oi = c("resid", "resno", "chain",
                                                "elety", "eleno"),
                                      separator = "_")

    PDBid.uniq.atom.ids <- paste(curr.pdb.id, uniq.atom.ids, sep = "_")

    ##--- add in the hashes
    atoms.oi <- cbind(PDBid = curr.pdb.id,
                      PDBid.uniq.atom.ids = PDBid.uniq.atom.ids,
                      res.chain.resno.ids = res.chain.resno.ids,
                      uniq.atom.ids = uniq.atom.ids,
                      res.atom.ids = res.atom.ids,
                      atoms.oi,
                      stringsAsFactors = FALSE)


    ##----- retain only protein and water atoms
    ##--- determine protein, hetero, and water indices
    prot.het.h2o.idc <- ProtHetWatIndices(data = atoms.oi)
    ##--- remove the hetero atoms but retain the waters!
    atoms.oi <- atoms.oi[c(prot.het.h2o.idc$prot.idc,
                           prot.het.h2o.idc$h2o.idc), ]
    ##--- determine the protein and water indices to reflect the removal
    prot.het.h2o.idc <- ProtHetWatIndices(data = atoms.oi)
    prot.idc <- prot.het.h2o.idc$prot.idc
    h2o.idc  <- prot.het.h2o.idc$h2o.idc


    ## GET PROTEIN ATOMS WITHIN X ANGSTROMS OF WATER OXYGEN ATOMS --------------
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    atoms.dist <- as.matrix(dist(atoms.oi[, c("x","y","z")],
                                 method = "euclidean",
                                 diag = TRUE, upper = TRUE))
    diag(atoms.dist) <- NA

    ##--- subset of water (rows) -- protein (columns) distances
    h2o.prot.dists <- atoms.dist[h2o.idc, prot.idc]

    ##--- check for waters within defined distances
    h2o.prot.dists.tf <- h2o.prot.dists <= h2o.prot.dist.max
    bound.h2o.dists.tf <- h2o.prot.dists <= bound.h2o.dist.max

    ##--- update df.dataset with number of waters within criteria
    num.h2o.prot.dists.lte.max <- sum(rowSums(h2o.prot.dists.tf) > 0)
    num.bound.h2o.dist.lte.max <- sum(rowSums(bound.h2o.dists.tf) > 0)
    df.dataset[curr.pdb.idx, c(15)] <- num.h2o.prot.dists.lte.max
    df.dataset[curr.pdb.idx, c(16)] <- num.bound.h2o.dist.lte.max
    if ( num.bound.h2o.dist.lte.max == 0 ) {
      mess <- paste("   NOTE: ", curr.pdb.id, " contains no bound waters ",
                    "within ", bound.h2o.dist.max, " Angstroms of the protein.",
                    " Moving to the next structure.", sep="")
      message(mess)
      next()
    }


    ##----- water indices within water-protein maximum distance
    h2o.near.prot.idc <- as.vector(which(rowSums(h2o.prot.dists.tf) >= 1))


    ## CALCULATE THE SOLVENT ACCESSIBLE SURFACE AREA ---------------------------
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    SASA.diff <- FreeSASA.diff(atoms.oi, probeRadius=probeRadius)


    ##----- merge atoms.oi with SASA.data
    atoms.oi <- merge(x = atoms.oi, y = SASA.diff,
                      by.x = "uniq.atom.ids", by.y = "uniq.atom.ids",
                      all = TRUE, sort = FALSE)


    ## NUMBER OF RESIDUES AND SOLVENT ACCESSIBLE/EXPOSED RESIDUES --------------
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # ##--- number of residues
    # type.ATOM.tf <- atoms.oi$type == "ATOM"
    # atoms.oi.prot <- atoms.oi[type.ATOM.tf, ]
    # num.res <- length(unique(atoms.oi.prot$res.chain.resno.ids))
    # ##--- surface exposed (solvent accessible) residues
    # SurExp.res.atoms.tf <- atoms.oi.prot$SASA.prot > 0
    # SurExp.res <- atoms.oi.prot$res.chain.resno.ids[SurExp.res.atoms.tf]
    # Buried.res <- atoms.oi.prot$res.chain.resno.ids[!SurExp.res.atoms.tf]
    # Buried.res <- Buried.res[!(Buried.res %in% SurExp.res)]
    # ##--- number surface exposed (solvent accessible) and buried
    # num.res.SurExp <- length(unique(SurExp.res))
    # num.res.buried <- num.res - num.res.SurExp
    # ##--- percent surface exposed (solvent accessible) residues
    # pct.res.SurExp <- num.res.SurExp / num.res
    #
    # ##--- SASA values
    # SASA.total <- sum(atoms.oi$SASA.prot, na.rm=TRUE)
    # ##--- SASA lost
    # SASA.lost <- sum(atoms.oi$SASA.lost, na.rm=TRUE)
    # ##--- percent SASA exposed
    # pct.SASA.exposed <- (SASA.total - SASA.lost) / SASA.total
    ##--- update the df.dataset data.frame
    # df.dataset[curr.pdb.idx, c(3:9)] <-
    #   c(num.res, num.res.buried, num.res.SurExp, pct.res.SurExp,
    #     SASA.total, SASA.lost, pct.SASA.exposed)
    #


    ## NUMBER OF RESIDUES AND SOLVENT ACCESSIBLE/EXPOSED RESIDUES --------------
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    type.ATOM.tf <- atoms.oi$type == "ATOM"
    atoms.oi.prot <- atoms.oi[type.ATOM.tf, ]
    ##--- surface exposed (solvent accessible) atoms
    SurExp.res.atoms.tf <- atoms.oi.prot$SASA.prot > 0
    residueSASA.data <- getResidueData(atoms.oi.prot, SurExp.res.atoms.tf)
    ##--- update the df.dataset data.frame
    df.dataset[curr.pdb.idx, c(3:9)] <- as.vector(unlist(residueSASA.data))


    ## NUMBER OF PROTEIN ATOMS AND SOLVENT ACCESSIBLE/EXPOSED ATOMS --------------
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ##----- atoms surface exposed and buried
    ##--- number surface exposed (solvent accessible) and buried atoms
    num.prot.atoms <- sum(type.ATOM.tf)
    num.atom.SurExp <- sum(SurExp.res.atoms.tf)
    num.atom.buried <- sum(!SurExp.res.atoms.tf)
    ##--- percent surface exposed (solvent accessible) atoms
    pct.atom.SurExp <- num.atom.SurExp / num.prot.atoms
    ##--- update the df.dataset data.frame
    df.dataset[curr.pdb.idx, c(10:13)] <- c(num.prot.atoms, num.atom.SurExp,
                                            num.atom.buried, pct.atom.SurExp)


    ## get protein atoms near each water oxygen atom within -------------------
    ## user defined cutoff/distance
    # df.nearby.prot.atoms <- data.frame(stringsAsFactors = FALSE)
    # ##--- loop over all waters within the user defined distance from the protein
    # for (h2o.oi in h2o.near.prot.idc) {
    #
    #   prot.atoms.oi.idc <- as.vector(which(h2o.prot.dists.tf[h2o.oi, ]))
    #   nearby.prot.atoms <- atoms.oi[prot.atoms.oi.idc, ]
    #
    #   ##--- get the distances and identify the minimum distance
    #   distances <- h2o.prot.dists[h2o.oi, prot.atoms.oi.idc]
    #   dist.is.min <- SASA.and.minDist <- rep(FALSE, length(distances))
    #   dist.is.min[which.min(distances)] <- TRUE
    #
    #   ##--- determine minimum distance with solvent accessibility
    #   SASA.minDist.idc <- 1:nrow(nearby.prot.atoms)
    #   ##- those with SASA
    #   has.SASA.idc <- SASA.minDist.idc[nearby.prot.atoms$SASA.prot > 0.0]
    #   SASA.and.minDist[has.SASA.idc[which.min(distances[has.SASA.idc])]] <- TRUE
    #
    #   ##--- add the nearby water information
    #   nearby.prot.atoms <- cbind(nearby.prot.atoms,
    #                              distances=distances,
    #                              dist.is.min=dist.is.min,
    #                              SASA.and.minDist=SASA.and.minDist,
    #                              h2o.atom.ids=atoms.oi[h2o.idc[h2o.oi], "uniq.atom.ids"],
    #                              h2o.x=atoms.oi[h2o.idc[h2o.oi], "x"],
    #                              h2o.y=atoms.oi[h2o.idc[h2o.oi], "y"],
    #                              h2o.z=atoms.oi[h2o.idc[h2o.oi], "z"],
    #                              stringsAsFactors=FALSE)
    #
    #   ##--- add current information to existing
    #   df.nearby.prot.atoms <- rbind(df.nearby.prot.atoms,
    #                                 nearby.prot.atoms,
    #                                 stringsAsFactors=FALSE)
    #
    # }
    list.nearby.prot.atoms <- lapply(X = h2o.near.prot.idc,
                                     FUN = getProtAtomsNearWater,
                                     h2o.idc = h2o.idc,
                                     atoms.oi = atoms.oi,
                                     h2o.prot.dists = h2o.prot.dists,
                                     h2o.prot.dists.tf = h2o.prot.dists.tf)

    df.nearby.prot.atoms <- do.call(rbind.data.frame, list.nearby.prot.atoms)


    ##----- add nearby.prot.atoms information to the df.results
    entry.end <- nrow(df.nearby.prot.atoms) + entry.start - 1
    df.results[entry.start:entry.end, ] <- df.nearby.prot.atoms
    entry.start <- entry.end + 1


    ##----- retain solvent exposed atoms and those within bound.h2o.dist.max (default: 4.0)
    SASA.and.boundDist.tf <- ((df.nearby.prot.atoms$SASA.prot > 0.0) &
                                (df.nearby.prot.atoms$distances <= bound.h2o.dist.max))
    df.SurExp.prot.SurBound.h2o <- df.nearby.prot.atoms[SASA.and.boundDist.tf, ]


    ## number of each ATOM type ------------------------------------------------
    ##--- number of atoms of each atom type
    all.AtomTypes <- atoms.oi.prot$res.atom.ids
    df.AtomTypes.all[curr.pdb.idx, ] <- getAtomTypeCounts(all.AtomTypes)

    ##--- number of buried atoms of each atom type
    buried.AtomTypes <- atoms.oi.prot$res.atom.ids[!SurExp.res.atoms.tf]
    df.AtomTypes.buried[curr.pdb.idx, ] <- getAtomTypeCounts(buried.AtomTypes)

    ##--- number of surface exposed (solvent accessible) atoms of each atom type
    SurExp.AtomTypes <- atoms.oi.prot$res.atom.ids[SurExp.res.atoms.tf]
    df.AtomTypes.SurExp[curr.pdb.idx, ] <- getAtomTypeCounts(SurExp.AtomTypes)

    ##--- number of surface exposed (solvent accessible) atom types NEARBY waters (h2o.prot.dist.max; default: 6.0)
    h2o.nearby.AtomTypes.tf <- df.nearby.prot.atoms$SASA.prot > 0.0
    h2o.nearby.AtomTypes <- df.nearby.prot.atoms$res.atom.ids[h2o.nearby.AtomTypes.tf]
    df.AtomTypes.h2o.nearby[curr.pdb.idx, ] <- getAtomTypeCounts(h2o.nearby.AtomTypes)

    ##--- number of surface exposed (solvent accessible) atom types BOUND TO waters (bound.h2o.dist.max; default: 4.0)
    h2o.bound.AtomTypes.tf <- ( (df.nearby.prot.atoms$SASA.prot > 0.0) &
                                  (df.nearby.prot.atoms$distances <= 4.0) )
    h2o.bound.AtomTypes <- df.nearby.prot.atoms$res.atom.ids[h2o.bound.AtomTypes.tf]
    df.AtomTypes.h2o.bound[curr.pdb.idx, ] <- getAtomTypeCounts(h2o.bound.AtomTypes)

    ##--- number of surface exposed (solvent accessible) atom types INTERACTING WITH waters
    ##      (shortest distance between protein atoms and unique/individual water)
    h2o.inter.AtomTypes.tf <- df.nearby.prot.atoms$SASA.and.minDist
    h2o.inter.AtomTypes <- df.nearby.prot.atoms$res.atom.ids[h2o.inter.AtomTypes.tf]
    df.AtomTypes.h2o.inter[curr.pdb.idx, ] <- getAtomTypeCounts(h2o.inter.AtomTypes)


    ## number of each RESIDUE type ---------------------------------------------
    ##--- number of residues of each residue type
    all.ResTypes.uniq <- unique(atoms.oi.prot$res.chain.resno.ids)
    all.ResTypes <- substring(all.ResTypes.uniq, 1, 3)
    df.ResTypes.all[curr.pdb.idx, ] <- getResTypeCounts(all.ResTypes)

    ##--- number of surface exposed (solvent accessible) residues of each residue type
    SurExp.ResTypes.uniq <- unique(atoms.oi.prot$res.chain.resno.ids[SurExp.res.atoms.tf])
    SurExp.ResTypes <- substring(SurExp.ResTypes.uniq, 1, 3)
    df.ResTypes.SurExp[curr.pdb.idx, ] <- getResTypeCounts(SurExp.ResTypes)

    ##--- number of buried residues of each residue type
    buried.ResTypes.uniq <- all.ResTypes.uniq[!all.ResTypes.uniq %in% SurExp.ResTypes.uniq]
    buried.ResTypes <- substring(buried.ResTypes.uniq, 1, 3)
    df.ResTypes.buried[curr.pdb.idx, ] <- getResTypeCounts(buried.ResTypes)

    ##--- number of surface exposed (solvent accessible) side chain residue
    ##    types INTERACTING WITH waters (shortest distance between protein atoms
    ##    and unique/individual water)
    sc.h2o.atoms.tf <- (df.nearby.prot.atoms$elety %in% names.sidechain.atoms) +
      (df.nearby.prot.atoms$SASA.and.minDist) == 2
    sc.h2o.ResTypes <- df.nearby.prot.atoms$resid[sc.h2o.atoms.tf]
    df.ResTypes.sc.h2o[curr.pdb.idx, ] <- getResTypeCounts(sc.h2o.ResTypes)

    ##--- number of surface exposed (solvent accessible) backbone (mainchain) residue
    ##    types INTERACTING WITH waters (shortest distance between protein atoms
    ##    and unique/individual water)
    bb.h2o.atoms.tf <- (df.nearby.prot.atoms$elety %in% names.backbone.atoms) +
      (df.nearby.prot.atoms$SASA.and.minDist) == 2
    bb.h2o.ResTypes <- df.nearby.prot.atoms$resid[bb.h2o.atoms.tf]
    df.ResTypes.bb.h2o[curr.pdb.idx, ] <- getResTypeCounts(bb.h2o.ResTypes)

    ##--- number of surface exposed (solvent accessible) residues
    ##    types INTERACTING WITH waters (shortest distance between protein atoms
    ##    and unique/individual water)
    res.h2o.ResTypes <- df.nearby.prot.atoms$resid[df.nearby.prot.atoms$SASA.and.minDist]
    df.ResTypes.res.h2o[curr.pdb.idx, ] <- getResTypeCounts(res.h2o.ResTypes)

    ##--- number of surface bound waters within bound.h2o.dist.max of each residue type
    SurBound.h2o.tf <- ( (df.nearby.prot.atoms$distances <= bound.h2o.dist.max) &
                           (df.nearby.prot.atoms$SASA.and.minDist) )
    SurBound.h2o.ResTypes <- df.nearby.prot.atoms$resid[SurBound.h2o.tf]
    df.ResTypes.SurBound.h2o[curr.pdb.idx, ] <- getResTypeCounts(SurBound.h2o.ResTypes)

    ##--- update df.dataset with surface bound water, backbone, side chain, and residue counts
    df.dataset[curr.pdb.idx, c(16)] <- sum(df.ResTypes.SurBound.h2o[curr.pdb.idx, ])
    df.dataset[curr.pdb.idx, c(17:19)] <- c(sum(df.ResTypes.bb.h2o[curr.pdb.idx, ]),
                                            sum(df.ResTypes.sc.h2o[curr.pdb.idx, ]),
                                            sum(df.ResTypes.res.h2o[curr.pdb.idx, ])
    )


    ## NUMBER OF WATER-RESIDUE & WATER-ATOM INTERACTIONS -----------------------
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ##----- water-residue interactions
    h2o.res.inter.ids <- paste0(df.SurExp.prot.SurBound.h2o$res.chain.resno.ids,
                                "_",
                                df.SurExp.prot.SurBound.h2o$h2o.atom.ids)
    df.dataset[curr.pdb.idx, c(20)] <- length(unique(h2o.res.inter.ids))
    # h2o.res.inter.ids <- paste0(df.nearby.prot.atoms$res.chain.resno.ids,
    #                             "_",
    #                             df.nearby.prot.atoms$h2o.atom.ids)
    # df.dataset[curr.pdb.idx, "num.h2o.res.inter"] <- length(unique(h2o.res.inter.ids))


    ##----- water-atom interactions
    h2o.resAtom.inter.ids <- paste0(df.SurExp.prot.SurBound.h2o$uniq.atom.ids,
                                    "_",
                                    df.SurExp.prot.SurBound.h2o$h2o.atom.ids)
    df.dataset[curr.pdb.idx, c(21)] <- length(unique(h2o.resAtom.inter.ids))


    ##----- calculate duration of analysis
    df.dataset[curr.pdb.idx, "time"] <- difftime(Sys.time(),
                                                 time.start.curr.pdb.idx,
                                                 units = "secs")

  }


  ## number of solvent exposed residue types (amino acids) -------------------
  df.residue.hydro$num.res.tot <- colSums(df.ResTypes.all)

  ##--- number of each type of residue with solvent accessibility/exposure
  df.residue.hydro$num.res.SurExp <- colSums(df.ResTypes.SurExp)

  ##--- number of each type of residue with NO solvent accessibility/exposure
  df.residue.hydro$num.res.buried <- colSums(df.ResTypes.buried)

  ##--- number of each type of residue with backbone (mainchain) solvent
  df.residue.hydro$num.bb.h2o <- colSums(df.ResTypes.bb.h2o)

  ##--- number of each type of residue with sidechain solvent
  df.residue.hydro$num.sc.h2o <- colSums(df.ResTypes.sc.h2o)

  ##--- number of each type of residue with solvent accessibility/exposure and
  ##    water contact/inateraction
  df.residue.hydro$num.res.h2o <- colSums(df.ResTypes.res.h2o)

  ##--- number of surface bound waters
  df.residue.hydro$num.SurBound.h2o <- colSums(df.ResTypes.SurBound.h2o)



  ## CLEANUP df.results (REMOVE ROWS OF NA) ------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##------ remove NA rows from df.results
  na.rows.tf <- rowSums(is.na(df.results)) == 31
  df.results <- df.results[!na.rows.tf, ]

  ##------ add column names
  colnames(df.results) <- colnames(df.nearby.prot.atoms)


  ## STRUCTURES WITHOUT WATERS -------------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##----- determine PDBids with and without waters
  pdb.ids.w.h2o <- unique(df.results$PDBid)
  pdb.ids.wo.h2o <- pdb.ids[!(pdb.ids %in% pdb.ids.w.h2o)]


  ## WATERS CLOSEST (WITHIN bound.h2o.dist.max ANGSTROMS) TO THE SOLVENT EXPOSED PROTEIN --------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##----- add unique water names (PDBid_HOH_res#_chain#_elesy_atom#)
  h2o.atom.ids <- paste(df.results$PDBid, df.results$h2o.atom.ids, sep = "_")
  df.results$h2o.atom.ids <- h2o.atom.ids
  ##--- new column names
  AtomicHydro.colNames <- c("PDBid", "resid", "elety", "distances", "SASA.prot",
                            "SASA.hetatm", "SASA.lost", "PDBid.uniq.atom.ids",
                            "res.chain.resno.ids", "res.atom.ids")
  ##--- waters closest to the surface exposed (solvent accessible) protein
  surface.h2o.tf <- ( (df.results$SASA.prot > 0) &
                        (df.results$distances <= bound.h2o.dist.max) )
  df.AtomicHydro.surface <- df.results[surface.h2o.tf, AtomicHydro.colNames]


  ## ATOMIC HYDROPHILICITY -----------------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##----- calculate number of interacting waters and bound waters
  AtomTypes.count     <- colSums(df.AtomTypes.all)
  AtomTypes.SurExp    <- colSums(df.AtomTypes.SurExp)
  # AtomTypes.h2o.nearby <- colSums(df.AtomTypes.h2o.nearby)
  # AtomTypes.h2o.bound <- colSums(df.AtomTypes.h2o.bound)
  AtomTypes.h2o.inter <- colSums(df.AtomTypes.h2o.inter)

  ##----- calculate atomic surface exposure probability
  AtomTypes.prob.SurExp <- AtomTypes.SurExp / AtomTypes.count
  AtomTypes.prob.SurExp.3digits <- round(AtomTypes.prob.SurExp, digits = 3)

  ##----- atomic hydrophilicity values
  AtomicHydrophilicity.values <- AtomTypes.h2o.inter / AtomTypes.SurExp
  AtomicHydrophilicity.values <- round(AtomicHydrophilicity.values, digits = 3)
  AtomicHydrophilicity.values <- as.vector(unlist(AtomicHydrophilicity.values))

  ##----- calculate AtomType classes hydrophilicity
  df.AtomHydroTEMP <- data.frame("residueAtomType" = names(AtomTypes.h2o.inter),
                                 "surfaceOccurrences" = AtomTypes.SurExp,
                                 "hydratOccurrences" = AtomTypes.h2o.inter,
                                 "probSurfaceOccurrence" = AtomTypes.prob.SurExp
                                 )
  AT.hydratFract <- calcAtomClassHydrophilicity(df.AtomHydroTEMP)

  ##----- calculate estimation of hydration for an atom with unknown surface
  ##      exposure
  AT.hydratFract.estimates <- calcAtomHydrationEstimate(df.AtomHydroTEMP,
                                                        AT.hydratFract)
  AT.hydratFract.estimates <- round(AT.hydratFract.estimates, digits = 3)

  ##----- construct the atomic hydrophilicity data.frame
  df.AtomicHydrophilicityValues <- data.frame(names(AtomTypes.h2o.inter),
                                              as.integer(AtomTypes.SurExp),
                                              AtomTypes.prob.SurExp.3digits,
                                              as.integer(AtomTypes.h2o.inter),
                                              AtomicHydrophilicity.values,
                                              AT.hydratFract.estimates,
                                              stringsAsFactors = FALSE)

  colnames(df.AtomicHydrophilicityValues) <- c("residueAtomType",
                                               "surfaceOccurrences",
                                               "prob.surfaceOccurrence",
                                               "hydratOccurrences",
                                               "hydratFraction",
                                               "est.hydratFraction")


  ##----- add new values to the original hydrophilicity table
  dataset.name <- paste0(".", dataset, "_", probeRadius)

  HydrophilicityTable.new <- merge(x = HydrophilicityTable,
                                   y = df.AtomicHydrophilicityValues,
                                   by.x = "residueAtomType",
                                   by.y = "residueAtomType",
                                   all = TRUE, sort = FALSE)

  colnames(HydrophilicityTable.new) <- gsub(pattern = ".x",
                                            replacement = ".Kuhn",
                                            x = colnames(HydrophilicityTable.new),
                                            fixed = TRUE)

  colnames(HydrophilicityTable.new) <- gsub(pattern = ".y",
                                            replacement = dataset.name,
                                            x = colnames(HydrophilicityTable.new),
                                            fixed = TRUE)


  ## FINISH ------------------------------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hydrophilicity.eval.duration <- difftime(Sys.time(),
                                           time.start.hydrophilicity.evaluation,
                                           units = "mins")


  ##----- return the data
  list(PDB.info = df.dataset,
       SASA.results = df.results,
       df.AtomTypes.all = df.AtomTypes.all,   ## total number of AtomTypes
       df.AtomTypes.buried = df.AtomTypes.buried,   ## number of buried AtomTypes
       df.AtomTypes.SurExp = df.AtomTypes.SurExp,     ## number of surface exposed AtomTypes =
       df.AtomTypes.h2o.nearby = df.AtomTypes.h2o.nearby,   ## number of surface exposed AtomTypes within h2o.prot.dist.max (default 6 Ang) of an individual water
       df.AtomTypes.h2o.bound = df.AtomTypes.h2o.bound,     ## number of surface exposed AtomTypes within bound.h2o.dist.max (default 4 Ang) of an indvidual water
       df.AtomTypes.h2o.inter = df.AtomTypes.h2o.inter,     ## number of surface exposed AtomTypes with the shortest distance to an individual water
       df.residue.hydro = df.residue.hydro,
       HydrophilicityTable = HydrophilicityTable.new,
       AtomTypeClasses.hydratFract = AT.hydratFract,        ## AtomType classes hydration fraction (hydrophilicity)
       no.h2o = pdb.ids.wo.h2o,
       call = the.call,
       duration = hydrophilicity.eval.duration)

}


## getResidueData docs ---------------------------------------------------------
#' @title Number of Residues and Solvent Accessible/Exposed Residues
#' @description Calculate the number of residues and solvent exposed residues.
#' @details This function is called within [HydrophilicityEvaluation()] to
#'   provide general solvent accessibility data for the protein structure of
#'   interest.
#'
#' @param atoms.oi.prot The protein `data.frame` with the `SASA` and `SASA lost`
#'   values for each protein atom.
#' @param SurExp.res.atoms.tf `TRUE`/`FALSE` vector indicating if an atom is
#'   solvent exposed/accessible
#'
#' @return
#'   This function returns:
#'   * **num.res**: number of residues within the structure
#'   * **num.res.buried**: number of residues with _**NO**_ solvent accessible
#'       surface area
#'   * **num.res.SurExp**: number of residues with solvent accessible surface
#'       area
#'   * **pct.res.SurExp**: percentage of residues with solvent accessible
#'       surface area
#'   * **SASA.total**: total protein solvent accessible surface area;
#'       Angstroms^2^
#'   * **SASA.lost**: total protein solvent accessible surface area lost due to
#'       bound waters; Angstroms^2^
#'   * **pct.SASA.exposed**: percentage protein solvent accessible surface area
#'       \eqn{(SASA.total - SASA.lost) / SASA.total}
#'
#'   These values are returned in `df.residue.hydro` of the results of
#'   [HydrophilicityEvaluation()]
#'
#' @export
#'
#' @examples
#'   \dontrun{
#'   getResidueData(atoms.oi.prot = PDB.1hai.aoi.clean.SASA.prot,
#'         SurExp.res.atoms.tf = PDB.1hai.SurExp.res.atoms.tf)
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Hydrophilicity Evaluation" "Bound Water Environment"
#'
getResidueData <- function(atoms.oi.prot, SurExp.res.atoms.tf) {

  ##----- number of residues
  num.res <- length(unique(atoms.oi.prot$res.chain.resno.ids))

  ##----- surface exposed (solvent accessible) residues
  SurExp.res <- atoms.oi.prot$res.chain.resno.ids[SurExp.res.atoms.tf]
  Buried.res <- atoms.oi.prot$res.chain.resno.ids[!SurExp.res.atoms.tf]
  Buried.res <- Buried.res[!(Buried.res %in% SurExp.res)]

  ##----- residue surface exposed and buried
  ##--- number surface exposed (solvent accessible) and buried residues
  num.res.SurExp <- length(unique(SurExp.res))
  num.res.buried <- num.res - num.res.SurExp
  ##--- percent surface exposed (solvent accessible) residues
  pct.res.SurExp <- num.res.SurExp / num.res

  ##----- calculate the SASA values, SASA lost, and SASA percent exposed
  ##--- SASA values
  SASA.total <- sum(atoms.oi.prot$SASA.prot, na.rm = TRUE)
  ##--- SASA lost
  SASA.lost <- sum(atoms.oi.prot$SASA.lost, na.rm = TRUE)
  ##--- percent SASA exposed
  pct.SASA.exposed <- (SASA.total - SASA.lost) / SASA.total

  ##----- return the data
  list(num.res = num.res,
       num.res.buried = num.res.buried,
       num.res.SurExp = num.res.SurExp,
       pct.res.SurExp = pct.res.SurExp,
       SASA.total = SASA.total,
       SASA.lost = SASA.lost,
       pct.SASA.exposed = pct.SASA.exposed)

}


## getProtAtomsNearWater docs --------------------------------------------------
#' @title Number of Solvent Accessible/Exposed Protein Atoms Near a Water
#' @description Calculate the number of solvent exposed protein atoms near a
#'   water.
#' @details This function is called within [HydrophilicityEvaluation()] to
#'   determine protein atoms near each water oxygen.
#'
#'   This function is designed to work with the [base::lapply()] function and
#'   thus each `h2o.oi` is independently evaluated
#'
#' @param h2o.oi The index of the water of interest
#' @param h2o.idc The indices of the waters within the protein structure
#' @param atoms.oi The protein `data.frame` with the `SASA` and `SASA lost`
#'   values for each atom within the protein.
#' @param h2o.prot.dists Distance `matrix` for all water-protein through space
#'   distances
#' @param h2o.prot.dists.tf The `TRUE`/`FALSE` matrix indicating if the protein-
#'   water distances are less than or equal to the user defined cutoff value
#'   denoted by the `h2o.prot.dist.max` parameter for
#'   [HydrophilicityEvaluation()]. From [HydrophilicityEvaluation()]: the
#'   maximum distance between the water oxygen atoms and the protein for
#'   consideration in the determination for hydrophilicity values; default: 6.0
#'
#' @return This function returns a `data.frame` with:
#'   * **nearby.prot.atoms**: protein atoms within the user specified distance
#'     of a water's oxygen atom
#'   * **distances**: The distance -- in Angstroms -- from the water to the
#'       closest solvent accessible protein atom so long as the distance is
#'       equal to or less than the user provided value; see `h2o.prot.dists.tf`
#'       above
#'   * **dist.is.min**: ; see `h2o.prot.dists.tf`
#'       above
#'   * **SASA.and.minDist**: `TRUE`/`FALSE` indicating if the protein atom is
#'       _**BOTH**_ solvent accessible and at least the user defined number
#'       of Angstroms from a water's oxygen atom; see `h2o.prot.dists.tf`
#'       above
#'   * **h2o.atom.ids**: Unique water atom ID
#'   * **h2o.x**: Atom coordinate `X` for the water's oxygen atom
#'   * **h2o.y**: Atom coordinate `Y` for the water's oxygen atom
#'   * **h2o.z**: Atom coordinate `Z` for the water's oxygen atom
#'
#'   These values are returned in `df.nearby.prot.atoms` of the results of
#'   [HydrophilicityEvaluation()]
#'
#' @export
#'
#' @examples
#'   \dontrun{
#'   getProtAtomsNearWater(h2o.oi = PDB.1hai.h2o.oi,
#'                         h2o.idc = PDB.1hai.clean.h2o.idc,
#'                         atoms.oi = PDB.1hai.aoi.clean.SASA,
#'                         h2o.prot.dists = PDB.1hai.h2o.prot.dists,
#'                         h2o.prot.dists.tf = PDB.1hai.h2o.prot.dists.tf)
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Hydrophilicity Evaluation" "Bound Water Environment"
#'
getProtAtomsNearWater <- function(h2o.oi,
                                  h2o.idc,
                                  atoms.oi,
                                  h2o.prot.dists,
                                  h2o.prot.dists.tf) {

  prot.atoms.oi.idc <- as.vector(which(h2o.prot.dists.tf[h2o.oi, ]))
  nearby.prot.atoms <- atoms.oi[prot.atoms.oi.idc, ]

  ##--- get the distances and identify the minimum distance
  distances <- h2o.prot.dists[h2o.oi, prot.atoms.oi.idc]
  dist.is.min <- SASA.and.minDist <- rep(FALSE, length(distances))
  dist.is.min[which.min(distances)] <- TRUE

  ##--- determine minimum distance with solvent accessibility
  SASA.minDist.idc <- 1:nrow(nearby.prot.atoms)
  ##- those with SASA
  has.SASA.idc <- SASA.minDist.idc[nearby.prot.atoms$SASA.prot > 0.0]
  SASA.and.minDist[has.SASA.idc[which.min(distances[has.SASA.idc])]] <- TRUE

  ##--- add the nearby water information
  data.frame(nearby.prot.atoms,
             distances = distances,
             dist.is.min = dist.is.min,
             SASA.and.minDist = SASA.and.minDist,
             h2o.atom.ids = atoms.oi[h2o.idc[h2o.oi], "uniq.atom.ids"],
             h2o.x = atoms.oi[h2o.idc[h2o.oi], "x"],
             h2o.y = atoms.oi[h2o.idc[h2o.oi], "y"],
             h2o.z = atoms.oi[h2o.idc[h2o.oi], "z"],
             stringsAsFactors = FALSE
  )

}


## calcAtomClassHydrophilicity docs --------------------------------------------
#' @title Atom Class Hydration Fraction
#' @description Calculates the mean hydration value for atoms within a class.
#' @details This function is called within [HydrophilicityEvaluation()] to
#'   calculate the hydration fraction for the five atom classes listed in the
#'   *Value* section.
#'
#'   \deqn{(num surface exposures)/(num atom occurrences)}
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param df.AtomHydroTEMP The newly calculated (determined) atomic
#'   hydrophilicity values
#'
#' @return
#'   This function returns:
#'   * **hydratFraction.oxy.neut**: Neutral oxygen atoms; enter
#'     `names.resATs.oxy.neut` to see list of residue-atomtypes
#'   * **hydratFraction.oxy.neg**: Negative oxygen atoms; enter
#'     `names.resATs.oxy.neg` to see list of residue-atomtypes
#'   * **hydratFraction.nitro.neut**: Neutral nitrogen atoms; enter
#'     `names.resATs.nitro.neut` to see list of residue-atomtypes
#'   * **hydratFraction.nitro.pos**: Positive nitrogen atoms; enter
#'     `names.resATs.nitro.pos` to see list of residue-atomtypes
#'   * **hydratFraction.carb.sulf**: Carbon and sulfur atoms; enter
#'     `names.resATs.carb.sulf` to see list of residue-atomtypes
#'
#'   These values are returned in `HydrophilicityValues.AtomTypeClasses` of the
#'   results of [HydrophilicityEvaluation()]
#'
#' @examples
#'  \dontrun{
#'   calcAtomClassHydrophilicity(df.AtomHydroTEMP)
#'  }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Hydrophilicity Evaluation" "Bound Water Environment"
#'
#' @references
#'   Leslie A Kuhn, Craig A Swanson, Michael E Pique, John A Tainer,
#'   and Elizabeth D Getzof. Atomic and Residue Hydrophilicity in the Context of
#'   Folded Protein Structures. _PROTEINS: Structure, Function, and
#'   Genetics_, 1995, **23** (_4_), pp 536-547.
#'   [DOI: 10.1002/prot.340230408](http://doi.org/10.1002/prot.340230408)
#'   [PMID: 8749849](http://www.ncbi.nlm.nih.gov/pubmed/8749849)
#'
calcAtomClassHydrophilicity <- function(df.AtomHydroTEMP) {

  ## get data.frame OF INTEREST FUNCTION ---------------------------------------
  get.data.frame <- function(df, rows.of.interest) {

    ##----- determine rows of interest
    rows.tf <- df$residueAtomType %in% rows.of.interest

    ##----- return data.frame
    df[rows.tf, ]

  }


  ## CALCULATE THE ATOMTYPE CLASS HYDROPHILICITY VALUES ------------------------
  ##----- neutral oxygen
  df.oxy.neut <- get.data.frame(df = df.AtomHydroTEMP,
                                rows.of.interest = names.resATs.oxy.neut)

  hydratFract.oxy.neut <- mean(df.oxy.neut$hydratOccurrences /
    df.oxy.neut$surfaceOccurrences)

  ##----- negative oxygen
  df.oxy.neg <- get.data.frame(df = df.AtomHydroTEMP,
                               rows.of.interest = names.resATs.oxy.neg)

  hydratFract.oxy.neg <- mean(df.oxy.neg$hydratOccurrences /
    df.oxy.neg$surfaceOccurrences)

  ##----- neutral nitrogen
  df.nitro.neut <- get.data.frame(df = df.AtomHydroTEMP,
                                  rows.of.interest = names.resATs.nitro.neut)

  hydratFract.nitro.neut <- mean(df.nitro.neut$hydratOccurrences /
    df.nitro.neut$surfaceOccurrences)

  ##----- positive nitrogen
  df.nitro.pos <- get.data.frame(df = df.AtomHydroTEMP,
                                 rows.of.interest = names.resATs.nitro.pos)

  hydratFract.nitro.pos <- mean(df.nitro.pos$hydratOccurrences /
    df.nitro.pos$surfaceOccurrences)

  ##----- carbon and sulfur
  df.carb.sulf <- get.data.frame(df = df.AtomHydroTEMP,
                                 rows.of.interest = names.resATs.carb.sulf)

  hydratFract.carb.sulf <- mean(df.carb.sulf$hydratOccurrences /
    df.carb.sulf$surfaceOccurrences)


  ## RETURN RESULTS ------------------------------------------------------------
  list(hydratFraction.oxy.neut   = hydratFract.oxy.neut,
       hydratFraction.oxy.neg    = hydratFract.oxy.neg,
       hydratFraction.nitro.neut = hydratFract.nitro.neut,
       hydratFraction.nitro.pos  = hydratFract.nitro.pos,
       hydratFraction.carb.sulf  = hydratFract.carb.sulf)

}



## calcAtomHydrationEstimate docs ----------------------------------------------
#' @title Estimated Atomic Hydration Fraction
#' @description Calculates the estimated atomic hydration fraction for an atom
#'   with unknown surface exposure.
#' @details This function is called within [HydrophilicityEvaluation()] to
#'   calculate the estimated hydration of an atom with unknown surface exposure.
#'
#'   \deqn{(num surface exposures)/(num atom occurrences) * atom class
#'   hydration fraction}
#'
#'   _**NOTE**_: This is a non-public function.
#'
#' @param df.AtomHydroTEMP The newly calculated (determined) atomic
#'   hydrophilicity values
#' @param AT.hydratFract The `AtomTypeClasses.hydratFract` variable calculated
#'   with the [HydrophilicityEvaluation()] function; the mean hydration fraction
#'   for the AtomTypes
#'
#' @return
#'   This function returns the hydration estimate values in a string to the variable
#'   `AT.hydratFract.estimates` and are included in the `HydrophilicityTable`
#'   results of [HydrophilicityEvaluation()].
#'
#' @examples
#'  \dontrun{
#'   calcAtomHydrationEstimate(df.AtomHydroTEMP, AT.hydratFract)
#'  }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Hydrophilicity Evaluation" "Bound Water Environment"
#'
#' @references
#'   Leslie A Kuhn, Craig A Swanson, Michael E Pique, John A Tainer,
#'   and Elizabeth D Getzof. Atomic and Residue Hydrophilicity in the Context of
#'   Folded Protein Structures. _PROTEINS: Structure, Function, and
#'   Genetics_, 1995, **23** (_4_), pp 536-547.
#'   [DOI: 10.1002/prot.340230408](http://doi.org/10.1002/prot.340230408)
#'   [PMID: 8749849](http://www.ncbi.nlm.nih.gov/pubmed/8749849)
#'
calcAtomHydrationEstimate <- function(df.AtomHydroTEMP, AT.hydratFract) {


  ## THE SETUP -----------------------------------------------------------------
  estimates <- rep(NA, nrow(df.AtomHydroTEMP))
  residueAtomType <- df.AtomHydroTEMP$residueAtomType
  probSurfaceOccurrence <- df.AtomHydroTEMP$probSurfaceOccurrence


  ## CALCULATE THE ESTIMATES ---------------------------------------------------
  ##----- neutral oxygen atoms
  resATs.oi.tf <- residueAtomType %in% names.resATs.oxy.neut
  estimates[resATs.oi.tf] <- probSurfaceOccurrence[resATs.oi.tf] *
    AT.hydratFract$hydratFraction.oxy.neut

  ##----- negative oxygen atoms
  resATs.oi.tf <- residueAtomType %in% names.resATs.oxy.neg
  estimates[resATs.oi.tf] <- probSurfaceOccurrence[resATs.oi.tf] *
    AT.hydratFract$hydratFraction.oxy.neg

  ##----- neutral nitrogen atoms
  resATs.oi.tf <- residueAtomType %in% names.resATs.nitro.neut
  estimates[resATs.oi.tf] <- probSurfaceOccurrence[resATs.oi.tf] *
    AT.hydratFract$hydratFraction.nitro.neut

  ##----- positive nitrogen atoms
  resATs.oi.tf <- residueAtomType %in% names.resATs.nitro.pos
  estimates[resATs.oi.tf] <- probSurfaceOccurrence[resATs.oi.tf] *
    AT.hydratFract$hydratFraction.nitro.pos

  ##----- carbon and sulfur atoms
  resATs.oi.tf <- residueAtomType %in% names.resATs.carb.sulf
  estimates[resATs.oi.tf] <- probSurfaceOccurrence[resATs.oi.tf] *
    AT.hydratFract$hydratFraction.carb.sulf


  ## RETURN RESULTS ------------------------------------------------------------
  return(estimates)

}
