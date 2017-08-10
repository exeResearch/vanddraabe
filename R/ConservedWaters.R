## ConservedWaters.R
##
## dec-28-2015 (exe) created
## may-05-2016 (exe) ConservedWaters: copy PDB files to used and rejected
##                   directories
## may-05-2016 (exe) ConservedWaters: rewrote PDB quality evaluations
## jul-30-2016 (exe) ConservedWaters: remove modeled heavy atoms
## sep-16-2016 (exe) split original file to create ConservedWaters.R
## oct-07-2016 (exe) added chains of interest check for valid chain IDs
## oct-07-2016 (exe) the following functions replace blocks of code
##                   - getRCSBdata function
##                   - DetermineChainsOfInterest function
##                   - RetainChainsOfInterest function
##                   - RemoveHydrogenAtoms function
##                   - RemoveModeledAtoms function
##                   - RetainWatersWithinX function
## oct-07-2016 (exe) creates ALL and PASSED conserved waters PyMOL script files
## jan-23-2017 (exe) update user provided parameter checks
## feb-07-2017 (exe) updated documentation
## apr-18-2017 (exe) added check for same PDBids to indicate MDS trajectory
##                   being use for structures
## apr-25-2017 (exe) added ConservedWaters.MDS()
## jul-25-2017 (exe) updated documentation
## jul-28-2017 (exe) updated NormalizedBvalue and Mobility calls;
##                   changed Bvalue parameter to Bvalues to match code
## jul-31-2017 (exe) updated the ConservedWaters(), ConservedWaters.MDS() and
##                   ConservedWaterStats() documentation
## aug-02-2017 (exe) added prefix to filenames.used variable
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




## conserved waters docs -------------------------------------------------------
#' @title Conserved Crystallographic Waters
#' @description Identifies conserved crystallographic waters from a collection
#'   of PDBs.
#' @details Only atoms within (less than or equal to) 5.10 Angstroms of the
#'   protein structures are included.
#'
#' @param prefix Directory of aligned structures; string.
#' @param cluster Oxygen atoms within 2.4 Angstroms or less of each other are
#'   considered a cluster; numeric. Default value is 2.4 Angstroms.
#' @param mobility A normalization method to identify the amount of variance an
#'   atom has within a structure; numeric. Calculated mobility values equal to
#'   or greater than the provided value will be removed from analysis. Default
#'   value is 2.0. See [Mobility()] for more information.
#' @param nBvalue The number of standard deviations from the mean for the water
#'   oxygens' B-values within the structure of interest; numeric. Calculated
#'   mobility values equal to or greater than the provided value will be removed
#'   from analysis. Default value is 1.0. See [NormalizedBvalue()] for
#'   more information.
#' @param chain The chain to examine. The user can define "first" and the first
#'   chain alphabetically will be selected; this is the default. Defining "all"
#'   will result in all chains being explored. Alternatively the user can define
#'   individual the chains to include in the analysis; for example, `c("A",
#'   "B", "C")`. When defining chains, the chain designation _**must
#'   be characters**_.
#' @param prot.h2o.dist.min The minimum distance (in Angstroms) between the
#'   protein and waters to be considered for the conserved water clusters. Water
#'   oxygen atoms greater than this distance are removed from the analysis.
#'   Default value is 5.10 Angstroms.
#' @param cluster.method Method of clustering the waters; default is "complete".
#'   Any other method accepted by the [stats::hclust()] or
#'   [fastcluster::hclust()] functions are appropriate. The original method used
#'   by Sanschagrin and Kuhn is the complete linkage clustering method and is
#'   the default. Other options include "ward.D" (equivilant to the only Ward
#'   option in `R` versions 3.0.3 and earlier), "ward.D2" (implements Ward's
#'   1963 criteria; see Murtagh and Legendre 2014), or "single" (related to the
#'   minimal spanning tree method and adopts a "friend of friends" clustering
#'   method). Please see [fastcluster::hclust()] for additional and complete
#'   information regarding clustering explanations.
#' @param PDBinfo The PDB information for all structures within the analysis.
#'   This information is obtained using the [getRCSBdata()] function.
#' @param filename The filename prefix for the returned results. Default is
#'   "ProteinSystem"
#'
#' @return
#'   This function returns:
#'   * **h2o.cluster.all**: Clusters constructed from all waters present in
#'     the aligned PDB structures.
#'   * **h2o.cluster.passed**: Clusters constructed from waters that passed
#'     the [Mobility()] and [NormalizedBvalue()] evaluations.
#'   * **h2o.cluster.summary**: Summary of water clusters
#'   * **Excel workbook**: containing the Cluster Statistics, Cluster Summaries
#'     for **all** and **passed** waters, Occurrence Summaries for **all** and
#'     **passed** waters, and the Initial Water Data data as individual tabs
#'   * **call**: The user provided parameters for the function
#'
#' @export
#'
#' @import bio3d
#' @importFrom stats aggregate dist
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @references
#'   Paul C Sanschagrin and Leslie A Kuhn. Cluster analysis of
#'   consensus water sites in thrombin and trypsin shows conservation between
#'   serine proteases and contributions to ligand specificity. _Protein
#'   Science_, 1998, **7** (_10_), pp 2054-2064.
#'   [DOI: 10.1002/pro.5560071002](http://doi.org/10.1002/pro.5560071002)
#'   [PMID: 9792092](http://www.ncbi.nlm.nih.gov/pubmed/9792092)
#'   [WatCH webpage](http://www.kuhnlab.bmb.msu.edu/software/watch/index.html)
#'
#'   Hitesh Patel, Bjorn A. Gruning, Stefan Gunther, and Irmgard Merfort.
#'   PyWATER: a PyMOL plug-in to find conserved water molecules in proteins by
#'   clustering. _Bioinformatics_, 2014, **30** (_20_), pp 2978-2980.
#'   [DOI: 10.1093/bioinformatics/btu424](http://doi.org/10.1093/bioinformatics/btu424)
#'   [PMID: 24990608](http://www.ncbi.nlm.nih.gov/pubmed/24990608)
#'   [PyWATER on GitHub](https://github.com/hiteshpatel379/PyWATER/blob/master/README.rst)
#'
ConservedWaters <- function(prefix="", # "thrombin_fitlsq_1.0ang/",
                            cluster=2.4, mobility=2.0, nBvalue=1.0,
                            chain="first",
                            prot.h2o.dist.min=5.10,
                            cluster.method="complete",
                            PDBinfo,
                            filename="ProteinSystem") {


  ## THE SETUP -----------------------------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##----- the provided call
  the.call <- match.call()

  ##----- create the date and time portion of the filenames
  current.time <- Sys.time()
  now.date.time <- FileTimeStamp(current.time)

  ##----- rename user defined variables
  cutoff.nBvalue <- nBvalue
  cutoff.mobility <- mobility
  cutoff.cluster <- cluster + 0.001
  cutoff.prot.h2o.dist <- prot.h2o.dist.min + 0.001
  chains.explore <- toupper(chain)


  ## CHECK THE USER PROVIDED VALUES --------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##----- nBvalue, mobility, cluster, and prot.h2o.dist.min
  if ( !is.numeric(cutoff.nBvalue) ||
       is.nan(cutoff.nBvalue) ||
       (cutoff.nBvalue <= 0.0) ) { cutoff.nBvalue <- NULL }

  if ( !is.numeric(cutoff.mobility) ||
       is.nan(cutoff.mobility) ||
       (cutoff.mobility <= 0.0) ) { cutoff.mobility <- NULL }

  if ( !is.numeric(cutoff.cluster) ||
       is.nan(cutoff.cluster) ||
       (cutoff.cluster <= 0.0) ) {
    cutoff.cluster <- 2.401
    mess <- paste("Cluster cutoff must be positive and greater than 0.0. Being",
                  "set to the default of 2.4", sep=" ")
    message(mess)
  }

  if ( !is.numeric(cutoff.prot.h2o.dist) ||
       is.nan(cutoff.prot.h2o.dist) ||
       (cutoff.prot.h2o.dist <= 0.0) ) { cutoff.prot.h2o.dist <- NA }

  ##----- check the user provided chain(s) to explore
  chains.valid <- c("FIRST", "ALL", LETTERS)
  chains.explore.vs.valid <- match(chains.explore, chains.valid)
  chains.explore.vs.valid.NA <- anyNA(chains.explore.vs.valid)
  if ( chains.explore.vs.valid.NA == TRUE ) {
    chains.invalid.idx <- match(NA, chains.explore.vs.valid)
    chains.invalid <- chains.explore[chains.invalid.idx]
    chains.stop <- paste("The provided chains of: ", toString(chains.explore),
                         " contains the invalid chain identifier: ",
                         toString(chains.invalid), ". Please only use the ",
                         "following chain identifiers: ",
                         toString(chains.valid), sep = "")
    stop(chains.stop)
  }

  ##----- check the clustering method
  cluster.method <- check.cluster.method(cluster.method)

  ## READ IN THE PDB STRUCTURES ------------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##----- get list of PDB files within prefix
  pdb.location <- ReturnPDBfullPath(prefix)

  ##----- determine PDB IDs
  pdb.ids <- ExtractPDBids(pdb.location)
  ##--- determine if the pdb.ids are all the same
  num.pdb.structures <- length(pdb.ids)
  num.pdb.ids.unique <- length(unique(pdb.ids))
  if ( num.pdb.structures != num.pdb.ids.unique ) {
    mess <- paste("There are", num.pdb.structures, "PDB structures yet there",
                  "is only", num.pdb.ids.unique, "unique PDB structure names.",
                  "Thus, the base filename is used to identify the individual",
                  "structures.",sep = " ")
    message(mess)
    pdb.ids <- basename(pdb.location)

    ##--- check for flsq in name
    num.flsq <- sum(grepl(pattern = "flsq", x = pdb.ids))
    if ( num.flsq == num.pdb.structures ) {
      pdb.ids <- gsub(pattern = ".pdb_flsq.pdb", replacement = "", x = pdb.ids)
    } else {
      pdb.ids <- gsub(pattern = ".pdb", replacement = "", x = pdb.ids)
    }
  }


  ##----- determine if enough structures are provided
  if ( num.pdb.structures <= 3 ) {
    stop("Water cluster analysis requires more than 3 PDB (structure) file.")
  }

  ##----- which chains to explore?
  chains.oi <- DetermineChainsOfInterest(chains.explore)

  ##----- read in aligned PDBs and create waters data.frame for clustering
  message("----- Reading in the aligned structures _____")
  h2o.df <- NULL
  modeled.heavy.atoms.pdb.ids <- no.h2o.pdb.ids <- got.h2o.pdb.ids <- NULL
  num.h2o.per.pdb <- rep(NA, num.pdb.structures)
  for (curr.pdb in 1:num.pdb.structures) {

    ##--- read in PDB
    curr.pdb.id <- pdb.ids[curr.pdb]
    read.mess <- paste("Reading structure ", curr.pdb.id, "...", sep = "")
    message(read.mess)
    pdb <- bio3d::read.pdb(file = pdb.location[curr.pdb])
    atoms.oi <- pdb$atom

    ##----- determine the chain(s) of interest
    atoms.chains.oi <- RetainChainsOfInterest(atoms.oi,
                                              chains.explore,
                                              chains.oi)

    ##----- remove hydrogen atoms (yeah, they shouldn't be there.
    ##      but sometimes...)
    atoms.chains.oi <- RemoveHydrogenAtoms(atoms.chains.oi)

    ##----- remove modeled heavy atoms
    ##--- identify the protein, het, and waters atoms
    prot.het.h2o.idc.pre.length <-
      unlist(lapply(ProtHetWatIndices(data = atoms.chains.oi), length))
    atoms.num.pre <- nrow(atoms.chains.oi)
    ##--- remove modeled atoms
    atoms.chains.oi <- RemoveModeledAtoms(atoms.chains.oi)
    ##--- compare post and initial structures
    atoms.num <- nrow(atoms.chains.oi)
    if ( atoms.num < atoms.num.pre ) {
      prot.het.h2o.idc.post.length <-
        unlist(lapply(ProtHetWatIndices(data = atoms.chains.oi), length))

      modeled.heavy.atoms.pdb.ids <- append(modeled.heavy.atoms.pdb.ids,
                                            curr.pdb.id)

      ##-- number of modeled atom(s) for each type of atom
      prot.het.h2o.lengths <- prot.het.h2o.idc.pre.length -
        prot.het.h2o.idc.post.length
      ##-- make the message
      modeled.mess <- paste("There are",
                            prot.het.h2o.lengths[1],
                            "modeled protein atoms;",
                            prot.het.h2o.lengths[2],
                            "modeled hetro (non-water) atoms;",
                            prot.het.h2o.lengths[3],
                            "modeled water atoms.\n", sep = " ")
      message(modeled.mess)
    }

    ##--- identify the protein, het, and waters atoms (again...)
    prot.het.h2o.idc <- ProtHetWatIndices(data = atoms.chains.oi)
    prot.idc <- prot.het.h2o.idc$prot.idc
    het.idc <- prot.het.h2o.idc$het.idc
    h2o.idc <- prot.het.h2o.idc$h2o.idc

    ##--- calculate the pairwise distances
    atoms.dist <- as.matrix(dist(atoms.chains.oi[, c("x","y","z")],
                                 method = "euclidean",
                                 diag = TRUE, upper = TRUE))
    diag(atoms.dist) <- NA

    ##--- retain waters within user defined distance (Angstroms) of the protein
    if ( (is.finite(cutoff.prot.h2o.dist)) & (cutoff.prot.h2o.dist > 0.0) ) {
      ##--- the water indices for waters within user defined
      ##    Angstroms of the protein
      h2o.idc <- RetainWatersWithinX(atoms.dist,
                                     prot.het.h2o.idc,
                                     cutoff.prot.h2o.dist)
    }

    ##--- update the structure and distance matrix
    atom.idc <- sort(c(prot.idc, het.idc, h2o.idc))
    atoms.dist <- atoms.dist[atom.idc, atom.idc]
    atoms.chains.oi <- atoms.chains.oi[atom.idc, ]

    ##--- update the indices
    ##--- identify the protein, het, and waters atoms (again...)
    prot.het.h2o.idc <- ProtHetWatIndices(data = atoms.chains.oi)
    prot.idc <- prot.het.h2o.idc$prot.idc
    het.idc <- prot.het.h2o.idc$het.idc
    h2o.idc <- prot.het.h2o.idc$h2o.idc
    num.atoms <- nrow(atoms.chains.oi)

    ##--- number of atoms in each group
    # num.prot.atoms <- length(prot.idc)
    # num.het.atoms <- length(het.idc)
    num.h2o.per.pdb[curr.pdb] <- curr.num.h2o <- length(h2o.idc)

    ##--- check if the file has waters
    h2o.message <- paste(" - Structure", curr.pdb, "of", num.pdb.structures,
                         "-->>", curr.pdb.id, "has", curr.num.h2o, "waters.",
                         sep = " ")
    message(h2o.message)
    if ( curr.num.h2o > 0 ) {
      got.h2o.pdb.ids <- append(got.h2o.pdb.ids, curr.pdb.id)
    } else {
      no.h2o.pdb.ids <- append(no.h2o.pdb.ids, curr.pdb.id)
      next()
    }

    ##--- get residue and atom names
    names.residues <- atoms.chains.oi$resid
    ##- atom names for the het atoms are the element symbol
    names.atoms <- rep(NA, num.atoms)
    names.atoms[prot.idc] <- atoms.chains.oi$elety[prot.idc]
    names.atoms[het.idc]  <- atoms.chains.oi$elesy[het.idc]
    names.atoms[h2o.idc]  <- atoms.chains.oi$elety[h2o.idc]
    names.res.atoms <- paste(names.residues, names.atoms, sep = " ")

    ##--- calculate mobility and normalized b-value for each group
    ##-- protein atoms
    prot.mobility <- Mobility(Bvalues = atoms.chains.oi$b[prot.idc],
                              occupancy = atoms.chains.oi$o[prot.idc])
    prot.nBvalue <- NormalizedBvalue(Bvalues = atoms.chains.oi$b[prot.idc])
    ##-- het atoms
    het.mobility <- Mobility(Bvalues = atoms.chains.oi$b[het.idc],
                             occupancy = atoms.chains.oi$o[het.idc])
    het.nBvalue <- NormalizedBvalue(Bvalues = atoms.chains.oi$b[het.idc])
    ##-- water atoms
    h2o.mobility <- Mobility(Bvalues = atoms.chains.oi$b[h2o.idc],
                             occupancy = atoms.chains.oi$o[h2o.idc])
    h2o.nBvalue <- NormalizedBvalue(Bvalues = atoms.chains.oi$b[h2o.idc])

    ##--- add helper columns and the mobility and normalized b-value values
    ##    to the structure
    pdb.id <- pdb.ids[curr.pdb]
    ##- find current PDBid row in PDBinfo data.frame
    PDBinfo.curr.pdb.idx <- which(tolower(PDBinfo$structureId) %in% pdb.id)
    atoms.chains.oi <- cbind(atoms.chains.oi,
                            PDBid = pdb.id,
                            resolution = PDBinfo[PDBinfo.curr.pdb.idx, "resolution"],
                            rObserved = PDBinfo[PDBinfo.curr.pdb.idx, "rObserved"],
                            rFree = PDBinfo[PDBinfo.curr.pdb.idx, "rFree"],
                            mobility = NA,
                            nBvalue = NA,
                            stringsAsFactors = FALSE)
    ##-- mobility
    atoms.chains.oi$mobility[prot.idc] <- prot.mobility
    atoms.chains.oi$mobility[het.idc] <- het.mobility
    atoms.chains.oi$mobility[h2o.idc] <- h2o.mobility
    ##-- norm B-value
    atoms.chains.oi$nBvalue[prot.idc] <- prot.nBvalue
    atoms.chains.oi$nBvalue[het.idc] <- het.nBvalue
    atoms.chains.oi$nBvalue[h2o.idc] <- h2o.nBvalue

    ##--- perform bound water environment analysis
    bwe.colNames <- c("adn", "ahp.sum", "ahp.mu", "ahp.sd",
                      "hbonds",
                      "o.sum", "o.mu", "o.sd",
                      "b.exp.sum", "b.exp.mu", "b.exp.sd",
                      "mobility.sum", "mobility.mu", "mobility.sd",
                      "nBvalue.sum", "nBvalue.mu", "nBvalue.sd")
    ##-- protein atoms
    prot.bwe.results <- apply(atoms.dist[h2o.idc, ], 1,
                              FUN = BoundWaterEnvironment,
                              set.oi.idc = prot.idc,
                              # h2o.idc = h2o.idc,
                              names.atoms = names.atoms,
                              names.res.atoms = names.res.atoms,
                              structure = atoms.chains.oi,
                              radius = 3.6)
    df.prot.bwe.results <- data.frame(matrix(unlist(prot.bwe.results),
                                             nrow = length(h2o.idc),
                                             byrow = TRUE),
                                      stringsAsFactors = FALSE)
    colnames(df.prot.bwe.results) <- paste0("prot.", bwe.colNames)
    ##-- water atoms
    h2o.bwe.results <- apply(atoms.dist[h2o.idc, ], 1,
                              FUN = BoundWaterEnvironment,
                              set.oi.idc = h2o.idc,
                              # h2o.idc = h2o.idc,
                              names.atoms = names.atoms,
                              names.res.atoms = names.res.atoms,
                              structure = atoms.chains.oi,
                              radius = 3.6)
    df.h2o.bwe.results <- data.frame(matrix(unlist(h2o.bwe.results),
                                            nrow = length(h2o.idc),
                                            byrow = TRUE),
                                      stringsAsFactors = FALSE)
    colnames(df.h2o.bwe.results) <- paste0("h2o.", bwe.colNames)
    ##-- het atoms
    het.bwe.results <- apply(atoms.dist[h2o.idc, ], 1,
                             FUN = BoundWaterEnvironment,
                             set.oi.idc = het.idc,
                             # h2o.idc = h2o.idc,
                             names.atoms = names.atoms,
                             names.res.atoms = names.res.atoms,
                             structure = atoms.chains.oi,
                             radius = 3.6)
    df.het.bwe.results <- data.frame(matrix(unlist(het.bwe.results),
                                            nrow = length(h2o.idc),
                                            byrow = TRUE),
                                     stringsAsFactors = FALSE)
    colnames(df.het.bwe.results) <- paste0("het.", bwe.colNames)
    ##-- all atoms
    all.bwe.results <- apply(atoms.dist[h2o.idc, ], 1,
                             FUN = BoundWaterEnvironment,
                             set.oi.idc = c(prot.idc, h2o.idc, het.idc),
                             # h2o.idc = h2o.idc,
                             names.atoms = names.atoms,
                             names.res.atoms = names.res.atoms,
                             structure = atoms.chains.oi,
                             radius = 3.6)
    df.all.bwe.results <- data.frame(matrix(unlist(all.bwe.results),
                                            nrow = length(h2o.idc),
                                            byrow = TRUE),
                                     stringsAsFactors = FALSE)
    colnames(df.all.bwe.results) <- paste0("all.", bwe.colNames)
    ##-- put all the BoundWaterEnvironment results together
    df.bwe.results <- cbind(df.prot.bwe.results,
                            df.h2o.bwe.results,
                            df.het.bwe.results,
                            df.all.bwe.results,
                            stringsAsFactors = FALSE)
    ##-- add the BWE results to the PDB information
    bwe.placeholder <- data.frame(matrix(data = NA,
                                         nrow = num.atoms,
                                         ncol = ncol(df.bwe.results)),
                                         stringsAsFactors = FALSE)
    colnames(bwe.placeholder) <- colnames(df.bwe.results)
    atoms.chains.oi <- cbind(atoms.chains.oi, bwe.placeholder)
    atoms.chains.oi[h2o.idc, colnames(df.bwe.results)] <- df.bwe.results

    ##--- construct data.frame of ONLY water residues
    h2o.res <- atoms.chains.oi[h2o.idc, ]

    ##--- construct and assign unique rownames for each water residue
    h2o.rownames <- paste(pdb.id,
                          h2o.res$resid,
                          h2o.res$chain,
                          h2o.res$resno,
                          sep = "_")
    rownames(h2o.res) <- h2o.rownames

    ##--- append extracted water residues to 'water residue data.frame'
    h2o.df <- rbind(h2o.df, h2o.res)
  }


  ## PDB STRUCTURES SUMMARY ----------------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##----- summary of structures imported
  num.h2o.per.pdb.imported <- num.h2o.per.pdb[num.h2o.per.pdb > 0]
  num.h2o <- sum(num.h2o.per.pdb.imported)
  pdb.id.imported <- unique(h2o.df$PDBid)
  num.pdbs.got.h2o <- length(got.h2o.pdb.ids)
  num.pdbs.no.h2o <- length(no.h2o.pdb.ids)
  num.pdbs.modeled.heavy <- length(modeled.heavy.atoms.pdb.ids)
  ##--- construct messages
  message.import <- paste(num.pdb.structures,
                          "structures were read for water cluster",
                          "analysis containing", num.h2o,
                          "water molecules.", sep = " ")
  message.got.h2o <- paste("  -", num.pdbs.got.h2o,
                           "Strucutures have water molecules.",
                           sep = " ")
  message.no.h2o <- paste("  -", num.pdbs.no.h2o,
                          "Strucuture(s) do NOT have water molecules.",
                          sep = " ")
  ##--- print the messages
  message("\n----- Summary about imported structures _____")
  message(message.import)
  message(message.got.h2o)
  message(message.no.h2o)
  if ( num.pdbs.no.h2o > 0 ) {
    message.no.h2o.pdb.ids <- paste(no.h2o.pdb.ids, collapse = ", ")
    message.no.h2o.2 <- paste("This structure(s) does NOT have waters: ",
                              message.no.h2o.pdb.ids, sep = "")
    message(message.no.h2o.2)
  }
  if ( num.pdbs.got.h2o <= 1 ) {
    stop("Water cluster analysis requires more than 1 PDB (structure)
         file to have water molecules.")
  }

  # if ( num.pdbs.modeled.heavy > 0 ) {
  #   message.modeled.heavy <- paste(modeled.heavy.atoms.pdb.ids, collapse=", ")
  #   message.modeled.2 <- paste("This structure(s) has heavy atoms with an occupancy value of 0: ",
  #                              message.modeled.heavy, sep="")
  #   message(message.modeled.2)
  # }

  ##----- filter waters using mobility and normalized b-value
  mobility.keep.tf <- nBvalue.keep.tf <- rep(TRUE, num.h2o)
  ##--- mobility
  if ( is.null(cutoff.mobility) == FALSE ) {
    mobility.keep.tf <- h2o.df$mobility < cutoff.mobility
    h2o.keep.tf <- mobility.keep.tf
  }
  ##--- normalized b-value
  if ( is.null(cutoff.nBvalue) == FALSE ) {
    nBvalue.keep.tf <- h2o.df$nBvalue < cutoff.nBvalue
    h2o.keep.tf <- nBvalue.keep.tf
  }
  ##--- based on the user-provided cutoffs, determine the waters for
  ##    conservation analysis
  if ( (is.null(cutoff.mobility) == FALSE) &
       (is.null(cutoff.nBvalue) == FALSE) ) {
    h2o.keep.tf <- (mobility.keep.tf + nBvalue.keep.tf) == 2
  }

  ##----- add mobility, normalized b-value, and keep/passed-cutoff info to
  ##      data.frame
  h2o.df <- cbind(h2o.df,
                  mobility.keep = mobility.keep.tf,
                  nBvalue.keep = nBvalue.keep.tf,
                  passed.cutoffs = h2o.keep.tf,
                  stringsAsFactors = FALSE)

  ##----- keep PDB structures where 50% of the waters PASS the mobility
  ##      and normalized b-value checks
  message(paste("\n----- Checking for high quality structures using B-value",
          "Normalization and Mobility _____", sep = " ") )
  num.h2o.per.pdb.passed <- as.vector(unlist(aggregate(h2o.df$passed.cutoffs,
                                                       list(h2o.df$PDBid),
                                                       sum)[2]))
  pct.h2o.per.pdb.eval <- num.h2o.per.pdb.passed / num.h2o.per.pdb.imported
  ##--- identify structures that PASS the 50% rule
  pdb.keep.tf <- pct.h2o.per.pdb.eval >= 0.50
  pdb.keep <- pdb.id.imported[pdb.keep.tf]
  ##--- keep structures that PASS the 50% rule
  h2o.df$passed.cutoffs[!(h2o.df$PDBid %in% pdb.keep)] <- FALSE
  h2o.df.passed <- h2o.df[h2o.df$passed.cutoffs == TRUE, ]
  ##--- report the structures that did NOT pass the 50% rule
  if ( length(pdb.keep) != length(pdb.id.imported) ) {
  message.50pct <- paste("The following structures were removed from analysis",
                         "because 50% or greater of their waters the did NOT",
                         "pass the mobility or normalized b-value evaluations.",
                         sep = " ")
  message(message.50pct)
  message( paste(pdb.id.imported[!pdb.keep.tf], collapse = ", ") )
  } else {
    message.50pct <- paste("Water molecules from all imported structures will",
                           "be used in the water conservation analysis.",
                           sep = " ")
    message(message.50pct)
  }


  ## CLUSTER WATER MOLECULES ---------------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##----- clustering
  message("\n----- Clustering ALL waters from the provided structures _____")
  h2o.cluster.all <- ClusterWaters(data = h2o.df,
                                   cutoff.cluster,
                                   cluster.method = cluster.method)
  message(paste("----- Clustering waters that PASSED the B-value Normalization",
                "and Mobility _____"), sep = " ")
  h2o.cluster.passed <- ClusterWaters(data = h2o.df.passed,
                                      cutoff.cluster,
                                      cluster.method = cluster.method)

  ##----- merge the ALL and PASSED h2o.clusters.raw data.frames
  h2o.clusters.raw.all <- h2o.cluster.all$h2o.clusters.raw
  h2o.clusters.raw.all <- cbind(h2o.clusters.raw.all,
                                unique.name = rownames(h2o.clusters.raw.all),
                                stringsAsFactors = FALSE)
  h2o.clusters.raw.passed <- h2o.cluster.passed$h2o.clusters.raw
  h2o.clusters.raw.passed <- cbind(h2o.clusters.raw.passed,
                                   unique.name = rownames(h2o.clusters.raw.passed),
                                   stringsAsFactors = FALSE)
  h2o.raw.all.passed <- merge(x = h2o.clusters.raw.all,
                              y = h2o.clusters.raw.passed,
                              by.x = "unique.name",
                              by.y = "unique.name",
                              all.x = TRUE,
                              sort = FALSE)
  rownames(h2o.raw.all.passed) <- h2o.raw.all.passed$unique.name

  h2o.raw.all.passed <- h2o.raw.all.passed[order(rownames(h2o.raw.all.passed)), ]
  h2o.clusters.raw.all <- h2o.clusters.raw.all[order(rownames(h2o.clusters.raw.all)), ]

  h2o.raw.all.passed <- cbind(h2o.clusters.raw.all[, colnames(h2o.clusters.raw.all) != "unique.name"],
                              cluster.passed = h2o.raw.all.passed$cluster.y,
                              pct.conserved.passed = h2o.raw.all.passed$pct.conserved.y,
                              stringsAsFactors = FALSE)


  h2o.raw.all.passed <- h2o.raw.all.passed[order(h2o.raw.all.passed$cluster), ]





  ## CONSTRUCT CLUSTERED WATERS SUMMARY ----------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##----- construct the structure and water summary for the clusters
  message("----- Constructing the summary table _____\n")
  h2o.cluster.all.stats <- ConservedWaterStats(h2o.cluster = h2o.cluster.all,
                                               num.h2o.inital = num.h2o,
                                               num.pdbs.got.h2o = num.pdbs.got.h2o)
  h2o.cluster.passed.stats <- ConservedWaterStats(h2o.cluster = h2o.cluster.passed,
                                                  num.h2o.inital = num.h2o,
                                                  num.pdbs.got.h2o = sum(pdb.keep.tf))
  h2o.cluster.stats <- cbind(h2o.cluster.all.stats,
                             h2o.cluster.passed.stats,
                             stringsAsFactors = FALSE)
  message("----- Summary table written _____\n")


  ## WRITE CLUSTERED WATERS PDB FILES ------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##----- write out PDB file with waters (o=% of structures with this water,
  ##      b=mean of b-values)
  ##--- all the provided waters
  message("----- Writing the conserved waters to PDB files _____")
  h2o.all.filename <- paste(filename,
                            "_ConservedWaters_ALL_",
                            now.date.time, ".pdb", sep = "")
  write.conservedWaters.pdb(file = h2o.all.filename,
                            h2o.clusters.summary = h2o.cluster.all$h2o.clusters.summary)

  ##--- the passed waters
  h2o.passed.filename <- paste(filename,
                               "_ConservedWaters_PASSED_",
                               now.date.time, ".pdb", sep = "")
  write.conservedWaters.pdb(file = h2o.passed.filename,
                            h2o.clusters.summary = h2o.cluster.passed$h2o.clusters.summary)

  message("----- PDB files written _____\n")


  ## WRITE RESULTS TO EXCEL WORKBOOK -------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##----- create worksheet names
  CW.ClusterStatistics <- paste("ClusterStats", now.date.time, sep="_")
  CW.all.h2o.summ      <- paste("all_ClustSumm", now.date.time, sep="_")
  CW.all.occ.summ      <- paste("all_OccurSumm", now.date.time, sep="_")
  CW.pass.h2o.summ     <- paste("pass_ClustSumm", now.date.time, sep="_")
  CW.pass.occ.summ     <- paste("pass_OccurSumm", now.date.time, sep="_")
  CW.init.h2o.data     <- paste("InitWaterData", now.date.time, sep="_")

  ##----- construct workbook name
  filename.xlsx <- paste(filename, "_DATA_RESULTS.xlsx", sep="")
  ##--- open existing excel workbook if available
  if ( file.exists(filename.xlsx) == TRUE ) {
    results.wb <- openxlsx::loadWorkbook(file=filename.xlsx)
  } else {  ##--- construct the workbook
    results.wb <- openxlsx::createWorkbook()
  }

  ##----- write the results to an excel workbook
  message("----- Writing results to Excel workbook _____")

  ##--- water cluster summary statistics
  results.wb <- oxClusterStatsSheet(wb.name = results.wb,
                                    sheet.name = CW.ClusterStatistics,
                                    df = h2o.cluster.stats)

  ##--- all clustered waters
  ##- cluster summary
  results.wb <- oxClusterSummarySheet(wb.name = results.wb,
                                      sheet.name = CW.all.h2o.summ,
                                      df = h2o.cluster.all$h2o.clusters.summary)

  ##- water occurrence summary
  results.wb <- oxWaterOccurrenceSheet(wb.name = results.wb,
                                       sheet.name = CW.all.occ.summ,
                                       df = h2o.cluster.all$h2o.occurrence)

  ##-- passed clustered waters
  ##- cluster summary
  results.wb <- oxClusterSummarySheet(wb.name = results.wb,
                                      sheet.name = CW.pass.h2o.summ,
                                      df = h2o.cluster.passed$h2o.clusters.summary)

  ##- water occurrence summary
  results.wb <- oxWaterOccurrenceSheet(wb.name = results.wb,
                                       sheet.name = CW.pass.occ.summ,
                                       df = h2o.cluster.passed$h2o.occurrence)

  ##-- the initial coordinates
  results.wb <- oxInitWaterDataSheet(wb.name = results.wb,
                                     sheet.name = CW.init.h2o.data,
                                     df = h2o.raw.all.passed)

  ##--- write the workbook
  openxlsx::saveWorkbook(results.wb, filename.xlsx, overwrite=TRUE)

  message("----- Results written to Excel workbook _____\n")


  ## RESULTS TO USER -----------------------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##----- filenames and date-time information
  filenames.used <- list(prefix=prefix,
                         conserved.all=h2o.all.filename,
                         conserved.passed=h2o.passed.filename,
                         xlsx=filename.xlsx,
                         when=now.date.time)

  ##----- return the results
  message("----- Done! _____")
  list(h2o.cluster.stats = h2o.cluster.stats,
       h2o.cluster.all = h2o.cluster.all,
       h2o.cluster.passed = h2o.cluster.passed,
       h2o.clusters.summary = h2o.raw.all.passed,  ## merged cluster number and conservation percentage
       filenames.used = filenames.used,
       MDS = FALSE,
       call = the.call
  )
}




## conserved MDS waters docs ---------------------------------------------------
#' @title Conserved Molecular Dynamics Simulation Waters
#' @description Identifies conserved molecular dynamics simulation (MDS) waters
#'   from a collection of PDBs.
#' @details Only atoms within (less than or equal to) 5.10 Angstroms of the
#'   protein structures are included.
#'
#' @param prefix Directory of aligned structures; string.
#' @param cluster Oxygen atoms within 2.4 Angstroms or less of each other are
#'   considered a cluster; numeric. Default value is 2.4 Angstroms.
#' @param chain The chain to examine. The user can define "first" and the first
#'   chain alphabetically will be selected; this is the default. Defining "all"
#'   will result in all chains being explored. Alternatively the user can define
#'   individual the chains to include in the analysis; for example, `c("A",
#'   "B", "C")`. When defining chains, the chain designation _**must
#'   be characters**_.
#' @param prot.h2o.dist.min The minimum distance (in Angstroms) between the
#'   protein and waters to be considered for the conserved water clusters. Water
#'   oxygen atoms greater than this distance are removed from the analysis.
#'   Default value is 5.10 Angstroms.
#' @param cluster.method Method of clustering the waters; default is "complete".
#'   Any other method accepted by the [stats::hclust()] or
#'   [fastcluster::hclust()] functions are appropriate. The original method used
#'   by Sanschagrin and Kuhn is the complete linkage clustering method and is
#'   the default. Other options include "ward.D" (equivilant to the only Ward
#'   option in `R` versions 3.0.3 and earlier), "ward.D2" (implements Ward's
#'   1963 criteria; see Murtagh and Legendre 2014), or "single" (related to the
#'   minimal spanning tree method and adopts a "friend of friends" clustering
#'   method). Please see [fastcluster::hclust()] for additional and complete
#'   information regarding clustering explanations.
#' @param filename The filename prefix for the returned results. Default is
#'   "ProteinSystem"
#'
#' @return
#'   This function returns:
#'   * **h2o.cluster.all**: Clusters constructed from all waters present in
#'     the aligned PDB structures.
#'   * **h2o.cluster.passed**: Clusters constructed from waters that passed
#'     the [Mobility()] and [NormalizedBvalue()] evaluations.
#'   * **h2o.cluster.summary**: Summary of water clusters
#'   * **Excel workbook**: containing the Cluster Statistics, Cluster Summaries
#'     for **all** and **passed** waters, Occurrence Summaries for **all** and
#'     **passed** waters, and the Initial Water Data data as individual tabs
#'   * **call**: The user provided parameters for the function
#'
#' @export
#'
#' @import bio3d
#' @importFrom stats dist
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @references
#'   Paul C Sanschagrin and Leslie A Kuhn. Cluster analysis of
#'   consensus water sites in thrombin and trypsin shows conservation between
#'   serine proteases and contributions to ligand specificity. _Protein
#'   Science_, 1998, **7** (_10_), pp 2054-2064.
#'   [DOI: 10.1002/pro.5560071002](http://doi.org/10.1002/pro.5560071002)
#'   [PMID: 9792092](http://www.ncbi.nlm.nih.gov/pubmed/9792092)
#'   [WatCH webpage](http://www.kuhnlab.bmb.msu.edu/software/watch/index.html)
#'
#'   Hitesh Patel, Bjorn A. Gruning, Stefan Gunther, and Irmgard Merfort.
#'   PyWATER: a PyMOL plug-in to find conserved water molecules in proteins by
#'   clustering. _Bioinformatics_, 2014, **30** (_20_), pp 2978-2980.
#'   [DOI: 10.1093/bioinformatics/btu424](http://doi.org/10.1093/bioinformatics/btu424)
#'   [PMID: 24990608](http://www.ncbi.nlm.nih.gov/pubmed/24990608)
#'   [PyWATER on GitHub](https://github.com/hiteshpatel379/PyWATER/blob/master/README.rst)
#'
ConservedWaters.MDS <- function(prefix="", # "thrombin_fitlsq_1.0ang/",
                                cluster=2.4,
                                chain="all",
                                prot.h2o.dist.min=5.10,
                                cluster.method="complete",
                                filename="ProteinSystem") {


  ## THE SETUP -----------------------------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##_ the provided call -----
  the.call <- match.call()

  ##_ create the date and time portion of the filenames -----
  current.time <- Sys.time()
  now.date.time <- FileTimeStamp(current.time)

  ##_ rename user defined variables -----
  cutoff.cluster <- cluster + 0.001
  cutoff.prot.h2o.dist <- prot.h2o.dist.min + 0.001
  chains.explore <- toupper(chain)


  ## CHECK THE USER PROVIDED VALUES --------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##_ cluster and prot.h2o.dist.min -----
  if ( !is.numeric(cutoff.cluster) ||
       is.nan(cutoff.cluster) ||
       (cutoff.cluster <= 0.0) ) {
    cutoff.cluster <- 2.401
    mess <- paste("Cluster cutoff must be positive and greater than 0.0. Being",
                  "set to the default of 2.4", sep=" ")
    message(mess)
  }

  if ( !is.numeric(cutoff.prot.h2o.dist) ||
       is.nan(cutoff.prot.h2o.dist) ||
       (cutoff.prot.h2o.dist <= 0.0) ) { cutoff.prot.h2o.dist <- NA }

  ##_ check the user provided chain(s) to explore -----
  chains.valid <- c("FIRST", "ALL", LETTERS)
  chains.explore.vs.valid <- match(chains.explore, chains.valid)
  chains.explore.vs.valid.NA <- anyNA(chains.explore.vs.valid)
  if ( chains.explore.vs.valid.NA == TRUE ) {
    chains.invalid.idx <- match(NA, chains.explore.vs.valid)
    chains.invalid <- chains.explore[chains.invalid.idx]
    chains.stop <- paste("The provided chains of: ", toString(chains.explore),
                         " contains the invalid chain identifier: ",
                         toString(chains.invalid), ". Please only use the ",
                         "following chain identifiers: ",
                         toString(chains.valid), sep = "")
    stop(chains.stop)
  }

  ##_ check the clustering method -----
  cluster.method <- tolower(cluster.method)
  cluster.methods.avail <- c("ward.d", "ward.d2", "single", "complete")
  if ( !any(cluster.method == cluster.methods.avail) ) {
    cluster.stop <- paste("Please select one of the following clustering",
                          "methods: \"complete\" (the default), \"ward.D\",",
                          "\"ward.D2\", or \"single\".", sep = " ")
    stop(cluster.stop)
  }
  if ( cluster.method == "ward.d" ) cluster.method <- "ward.D"
  if ( cluster.method == "ward.d2" ) cluster.method <- "ward.D2"


  ## READ IN THE PDB STRUCTURES ------------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##_ get list of PDB files within prefix -----
  pdb.location <- ReturnPDBfullPath(prefix)

  ##_ determine PDB IDs -----
  # pdb.ids <- basename.MDS.pdb(pdb.location)
  pdb.ids <- ExtractPDBids(pdb.location)
  ##__ determine if the pdb.ids are all the same -----
  num.pdb.structures <- length(pdb.ids)
  num.pdb.ids.unique <- length(unique(pdb.ids))
  if ( num.pdb.structures != num.pdb.ids.unique ) {
    mess <- paste("There are", num.pdb.structures, "PDB structures yet there",
                  "is only", num.pdb.ids.unique, "unique PDB structure names.",
                  "Thus, the base filename is used to identify the individual",
                  "structures.",sep = " ")
    message(mess)
    pdb.ids <- basename(pdb.location)

    ##__ check for flsq in name -----
    num.flsq <- sum(grepl(pattern = "flsq", x = pdb.ids))
    if ( num.flsq == num.pdb.structures ) {
      pdb.ids <- gsub(pattern = ".pdb_flsq.pdb", replacement = "", x = pdb.ids)
    } else {
      pdb.ids <- gsub(pattern = ".pdb", replacement = "", x = pdb.ids)
    }
  }


  ##_ determine if enough structures are provided -----
  if ( num.pdb.structures <= 3 ) {
    stop("Water cluster analysis requires more than 3 PDB (structure) file.")
  }

  ##_ which chains to explore? -----
  chains.oi <- DetermineChainsOfInterest(chains.explore)

  ##_ read in aligned PDBs and create waters data.frame for clustering -----
  message("----- Reading in the aligned structures _____")
  h2o.df <- NULL
  modeled.heavy.atoms.pdb.ids <- no.h2o.pdb.ids <- got.h2o.pdb.ids <- NULL
  num.h2o.per.pdb <- rep(NA, num.pdb.structures)
  for (curr.pdb in 1:num.pdb.structures) {

    ##__ read in PDB -----
    curr.pdb.id <- pdb.ids[curr.pdb]
    read.mess <- paste("Reading structure ", curr.pdb.id, "...", sep = "")
    message(read.mess)
    pdb <- bio3d::read.pdb(file = pdb.location[curr.pdb])
    atoms.oi <- pdb$atom

    ##__ determine the chain(s) of interest -----
    atoms.chains.oi <- RetainChainsOfInterest(atoms.oi,
                                              chains.explore,
                                              chains.oi)

    ##__ remove hydrogen atoms -----
    ## (yeah, they shouldn't be there. but sometimes...)
    atoms.chains.oi <- RemoveHydrogenAtoms(atoms.chains.oi)

    ##__ standardize residue names -----
    atoms.chains.oi$resid <- aaStandardizeNames(atoms.chains.oi$resid)

    ##__ rename ATOM to HETATM for waters -----
    h2o.residues.tf <- atoms.chains.oi$resid %in% names.waters
    atoms.chains.oi$type[h2o.residues.tf] <- "HETATM"

    ##__ rename ATOM to HETATM for possible ligands -----
    standard.AAs.tf <- atoms.chains.oi$resid %in% names.residues
    ligand.AAs.tf <- !h2o.residues.tf & !standard.AAs.tf
    atoms.chains.oi$type[ligand.AAs.tf] <- "HETATM"

    ##___ identify the protein, het, and waters atoms
    prot.het.h2o.idc.pre.length <-
      unlist(lapply(ProtHetWatIndices(data = atoms.chains.oi), length))
    atoms.num.pre <- nrow(atoms.chains.oi)

    ##___ identify the protein, het, and waters atoms (again...) -----
    prot.het.h2o.idc <- ProtHetWatIndices(data = atoms.chains.oi)
    prot.idc <- prot.het.h2o.idc$prot.idc
    het.idc <- prot.het.h2o.idc$het.idc
    h2o.idc <- prot.het.h2o.idc$h2o.idc

    ##___ calculate the pairwise distances -----
    atoms.dist <- as.matrix(dist(atoms.chains.oi[, c("x","y","z")],
                                 method = "euclidean",
                                 diag = TRUE, upper = TRUE))
    diag(atoms.dist) <- NA

    ##___ retain waters within user defined distance (Angstroms) of the protein -----
    if ( (is.finite(cutoff.prot.h2o.dist)) & (cutoff.prot.h2o.dist > 0.0) ) {
      ##--- the water indices for waters within user defined
      ##    Angstroms of the protein
      h2o.idc <- RetainWatersWithinX(atoms.dist,
                                     prot.het.h2o.idc,
                                     cutoff.prot.h2o.dist)
    }

    ##___ update the structure and distance matrix -----
    atom.idc <- sort(c(prot.idc, het.idc, h2o.idc))
    atoms.dist <- atoms.dist[atom.idc, atom.idc]
    atoms.chains.oi <- atoms.chains.oi[atom.idc, ]

    ##___ update the indices -----
    ##____ identify the protein, het, and waters atoms (again...) -----
    prot.het.h2o.idc <- ProtHetWatIndices(data = atoms.chains.oi)
    prot.idc <- prot.het.h2o.idc$prot.idc
    het.idc <- prot.het.h2o.idc$het.idc
    h2o.idc <- prot.het.h2o.idc$h2o.idc
    num.atoms <- nrow(atoms.chains.oi)

    ##___ number of atoms in each group -----
    # num.prot.atoms <- length(prot.idc)
    # num.het.atoms <- length(het.idc)
    num.h2o.per.pdb[curr.pdb] <- curr.num.h2o <- length(h2o.idc)

    ##___ check if the file has waters -----
    h2o.message <- paste(" - Structure", curr.pdb, "of", num.pdb.structures,
                         "-->>", curr.pdb.id, "has", curr.num.h2o, "waters.",
                         sep = " ")
    message(h2o.message)
    if ( curr.num.h2o > 0 ) {
      got.h2o.pdb.ids <- append(got.h2o.pdb.ids, curr.pdb.id)
    } else {
      no.h2o.pdb.ids <- append(no.h2o.pdb.ids, curr.pdb.id)
      next()
    }

    ##___ get residue and atom names -----
    names.residues <- atoms.chains.oi$resid
    ##- atom names for the het atoms are the element symbol
    names.atoms <- rep(NA, num.atoms)
    names.atoms[prot.idc] <- atoms.chains.oi$elety[prot.idc]
    names.atoms[het.idc]  <- atoms.chains.oi$elesy[het.idc]
    names.atoms[h2o.idc]  <- atoms.chains.oi$elety[h2o.idc]
    names.res.atoms <- paste(names.residues, names.atoms, sep = " ")

    ##___ add helper columns and the mobility and normalized b-value values -----
    ##    to the structure
    pdb.id <- pdb.ids[curr.pdb]
    ##- find current PDBid row in PDBinfo data.frame
    # PDBinfo.curr.pdb.idx <- which(tolower(PDBinfo$structureId) %in% pdb.id)
    atoms.chains.oi <- cbind(atoms.chains.oi,
                             PDBid = pdb.id,
                             # resolution = PDBinfo[PDBinfo.curr.pdb.idx, "resolution"],
                             # rObserved = PDBinfo[PDBinfo.curr.pdb.idx, "rObserved"],
                             # rFree = PDBinfo[PDBinfo.curr.pdb.idx, "rFree"],
                             # mobility = NA,
                             # nBvalue = NA,
                             stringsAsFactors = FALSE)
    # ##-- mobility
    # atoms.chains.oi$mobility[prot.idc] <- prot.mobility
    # atoms.chains.oi$mobility[het.idc] <- het.mobility
    # atoms.chains.oi$mobility[h2o.idc] <- h2o.mobility
    # ##-- norm B-value
    # atoms.chains.oi$nBvalue[prot.idc] <- prot.nBvalue
    # atoms.chains.oi$nBvalue[het.idc] <- het.nBvalue
    # atoms.chains.oi$nBvalue[h2o.idc] <- h2o.nBvalue



    ##___ perform bound water environment analysis -----
    bwe.colNames <- c("adn", "ahp.sum", "ahp.mu", "ahp.sd",
                      "hbonds")

    ##____ protein atoms -----
    prot.bwe.results <- apply(atoms.dist[h2o.idc, ], 1,
                              FUN = BoundWaterEnvironment.interact,
                              set.oi.idc = prot.idc,
                              names.atoms = names.atoms,
                              names.res.atoms = names.res.atoms,
                              radius = 3.6)
    df.prot.bwe.results <- data.frame(matrix(unlist(prot.bwe.results),
                                             nrow = length(h2o.idc), byrow = T),
                                      stringsAsFactors = FALSE)
    colnames(df.prot.bwe.results) <- paste0("prot.", bwe.colNames)
    ##____ water atoms -----
    h2o.bwe.results <- apply(atoms.dist[h2o.idc, ], 1,
                             FUN = BoundWaterEnvironment.interact,
                             set.oi.idc = h2o.idc,
                             names.atoms = names.atoms,
                             names.res.atoms = names.res.atoms,
                             radius = 3.6)
    df.h2o.bwe.results <- data.frame(matrix(unlist(h2o.bwe.results),
                                            nrow = length(h2o.idc), byrow = T),
                                     stringsAsFactors = FALSE)
    colnames(df.h2o.bwe.results) <- paste0("h2o.", bwe.colNames)
    ##____ het atoms -----
    het.bwe.results <- apply(atoms.dist[h2o.idc, ], 1,
                             FUN = BoundWaterEnvironment.interact,
                             set.oi.idc = het.idc,
                             names.atoms = names.atoms,
                             names.res.atoms = names.res.atoms,
                             radius = 3.6)
    df.het.bwe.results <- data.frame(matrix(unlist(het.bwe.results),
                                            nrow = length(h2o.idc), byrow = T),
                                     stringsAsFactors = FALSE)
    colnames(df.het.bwe.results) <- paste0("het.", bwe.colNames)
    ##____ all atoms -----
    all.bwe.results <- apply(atoms.dist[h2o.idc, ], 1,
                             FUN = BoundWaterEnvironment.interact,
                             set.oi.idc = c(prot.idc, h2o.idc, het.idc),
                             names.atoms = names.atoms,
                             names.res.atoms = names.res.atoms,
                             radius = 3.6)
    df.all.bwe.results <- data.frame(matrix(unlist(all.bwe.results),
                                            nrow = length(h2o.idc), byrow = T),
                                     stringsAsFactors = FALSE)
    colnames(df.all.bwe.results) <- paste0("all.", bwe.colNames)
    ##____ put all the BoundWaterEnvironment results together -----
    df.bwe.results <- cbind(df.prot.bwe.results,
                            df.h2o.bwe.results,
                            df.het.bwe.results,
                            df.all.bwe.results,
                            stringsAsFactors = FALSE)
    ##____ add the BWE results to the PDB information -----
    bwe.placeholder <- data.frame(matrix(data = NA,
                                         nrow = num.atoms,
                                         ncol = ncol(df.bwe.results)),
                                  stringsAsFactors = FALSE)
    colnames(bwe.placeholder) <- colnames(df.bwe.results)
    atoms.chains.oi <- cbind(atoms.chains.oi, bwe.placeholder)
    atoms.chains.oi[h2o.idc, colnames(df.bwe.results)] <- df.bwe.results

    ##___ construct data.frame of ONLY water residues -----
    h2o.res <- atoms.chains.oi[h2o.idc, ]

    ##___ construct and assign unique rownames for each water residue -----
    h2o.rownames <- paste(pdb.id,
                          h2o.res$resid,
                          h2o.res$chain,
                          h2o.res$resno,
                          sep = "_")
    rownames(h2o.res) <- h2o.rownames

    ##___ append extracted water residues to 'water residue data.frame' -----
    h2o.df <- rbind(h2o.df, h2o.res)
  }


  ## PDB STRUCTURES SUMMARY ----------------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##_ summary of structures imported -----
  num.h2o.per.pdb.imported <- num.h2o.per.pdb[num.h2o.per.pdb > 0]
  num.h2o <- sum(num.h2o.per.pdb.imported)
  pdb.id.imported <- unique(h2o.df$PDBid)
  num.pdbs.got.h2o <- length(got.h2o.pdb.ids)
  num.pdbs.no.h2o <- length(no.h2o.pdb.ids)
  num.pdbs.modeled.heavy <- length(modeled.heavy.atoms.pdb.ids)
  ##__ construct messages -----
  message.import <- paste(num.pdb.structures,
                          "structures were read for water cluster",
                          "analysis containing", num.h2o,
                          "water molecules.", sep = " ")
  message.got.h2o <- paste("  -", num.pdbs.got.h2o,
                           "Strucutures have water molecules.",
                           sep = " ")
  message.no.h2o <- paste("  -", num.pdbs.no.h2o,
                          "Strucuture(s) do NOT have water molecules.",
                          sep = " ")
  ##__ print the messages -----
  message("\n----- Summary about imported structures _____")
  message(message.import)
  message(message.got.h2o)
  message(message.no.h2o)
  if ( num.pdbs.no.h2o > 0 ) {
    message.no.h2o.pdb.ids <- paste(no.h2o.pdb.ids, collapse = ", ")
    message.no.h2o.2 <- paste("This structure(s) does NOT have waters: ",
                              message.no.h2o.pdb.ids, sep = "")
    message(message.no.h2o.2)
  }
  if ( num.pdbs.got.h2o <= 1 ) {
    stop("Water cluster analysis requires more than 1 PDB (structure)
              file to have water molecules.")
  }


  ## CLUSTER WATER MOLECULES ---------------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##_ clustering -----
  message("\n----- Clustering MDS waters from the provided structures _____")
  h2o.cluster.all <- ClusterWaters.MDS(data = h2o.df,
                                   cutoff.cluster,
                                   cluster.method = cluster.method)

  ##_ merge the ALL and PASSED h2o.clusters.raw data.frames -----
  h2o.clusters.raw.all <- h2o.cluster.all$h2o.clusters.raw
  h2o.clusters.raw.all <- cbind(h2o.clusters.raw.all,
                                unique.name = rownames(h2o.clusters.raw.all),
                                stringsAsFactors = FALSE)
  h2o.clusters.raw.all <- h2o.clusters.raw.all[order(rownames(h2o.clusters.raw.all)), ]


  ## CONSTRUCT CLUSTERED WATERS SUMMARY ----------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##_ construct the structure and water summary for the clusters -----
  message("----- Constructing the summary table _____\n")
  h2o.cluster.all.stats <- ConservedWaterStats(h2o.cluster = h2o.cluster.all,
                                               num.h2o.inital = num.h2o,
                                               num.pdbs.got.h2o = num.pdbs.got.h2o)
  h2o.cluster.stats <- cbind(h2o.cluster.all.stats,
                             # h2o.cluster.passed.stats,
                             stringsAsFactors = FALSE)
  message("----- Summary table complete _____\n")


  ## WRITE CLUSTERED WATERS PDB FILES ------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##_ write out PDB file with waters (o=% of structures with this water, -----
  ##      b=mean of b-values)
  ##--- all the provided waters
  message("----- Writing the conserved MDS waters to PDB file _____")
  h2o.all.filename <- paste(filename,
                            "_ConservedWaters_ALL_",
                            now.date.time, ".pdb", sep = "")
  write.conservedWaters.pdb(file = h2o.all.filename,
                            h2o.clusters.summary = h2o.cluster.all$h2o.clusters.summary)

  message("----- PDB file written _____\n")


  ## WRITE RESULTS TO EXCEL WORKBOOK -------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##_ create worksheet names -----
  CW.ClusterStatistics <- paste("ClusterStats", now.date.time, sep="_")
  CW.all.h2o.summ      <- paste("MDS_ClustSumm", now.date.time, sep="_")
  CW.all.occ.summ      <- paste("MDS_OccurSumm", now.date.time, sep="_")
  CW.init.h2o.data     <- paste("InitWaterData", now.date.time, sep="_")

  ##_ construct workbook name -----
  filename.xlsx <- paste(filename, "_DATA_RESULTS.xlsx", sep="")
  ##__ open existing excel workbook if available -----
  if ( file.exists(filename.xlsx) == TRUE ) {
    results.wb <- openxlsx::loadWorkbook(file=filename.xlsx)
  } else {  ##--- construct the workbook
    results.wb <- openxlsx::createWorkbook()
  }

  ##_ write the results to an excel workbook -----
  message("----- Writing results to Excel workbook _____")

  ##__ water cluster summary statistics -----
  results.wb <- oxClusterStatsSheet(wb.name = results.wb,
                                    sheet.name = CW.ClusterStatistics,
                                    df = h2o.cluster.stats)

  ##__ all clustered waters -----
  ##- cluster summary
  results.wb <- oxClusterSummarySheet(wb.name = results.wb,
                                      sheet.name = CW.all.h2o.summ,
                                      df = h2o.cluster.all$h2o.clusters.summary)

  ##- water occurrence summary
  results.wb <- oxWaterOccurrenceSheet(wb.name = results.wb,
                                       sheet.name = CW.all.occ.summ,
                                       df = h2o.cluster.all$h2o.occurrence)

  ##__ write the workbook -----
  openxlsx::saveWorkbook(results.wb, filename.xlsx, overwrite=TRUE)

  message("----- Results written to Excel workbook _____\n")


  ## RESULTS TO USER -----------------------------------------------------------
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##_ filenames and date-time information -----
  filenames.used <- list(prefix=prefix,
                         conserved.all=h2o.all.filename,
                         xlsx=filename.xlsx,
                         when=now.date.time)


  ##_ return the results -----
  message("----- Done! _____")
  list(h2o.cluster.stats = h2o.cluster.stats,
       h2o.cluster.all = h2o.cluster.all,
       # h2o.clusters.summary = h2o.raw.all.passed,  ## merged cluster number and conservation percentage
       filenames.used = filenames.used,
       MDS = TRUE,
       call = the.call
  )
}




## conserved water statistics docs ---------------------------------------------
#' @title Conserved Water Statistics
#' @description Calculates the Conserved Water Statistics for
#'   [ConservedWaters()]
#' @details Calculates the statistics for each conserved water analysis
#'   performed by [ConservedWaters()]. This summary information is useful for
#'   timings information and is written to the Excel workbook.
#'
#' @param h2o.cluster Conserved water cluster
#' @param num.h2o.inital Number of initial waters
#' @param num.pdbs.got.h2o Number of PDB structures with waters
#'
#' @return
#'   A table with the following information is returned:
#'   - Number of structures
#'   - Number of initial waters
#'   - Number of waters used in the calculations
#'   - Number of water clusters
#'   - Average water conservation
#'   - Number of conserved waters with
#'     + < 50\% conservation
#'     + 50 - 69\% conservation
#'     + 70 - 79\% conservation
#'     + 80 - 89\% conservation
#'     + 90 - 99\% conservation
#'     + 100\% conservation
#'   - Number of pairwise distances evaluated
#'   - Amount of memory used by the pairwise diatance matrix
#'   - Pairwise distance calculation time
#'   - Cluster centroid distance calculation time
#'
#' @family "vanddraabe utilities"
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
ConservedWaterStats <- function(h2o.cluster,
                                num.h2o.inital,
                                num.pdbs.got.h2o) {

  ##_ extract needed information from h2o cluster -----
  num.h2o.used <- as.integer(sum(h2o.cluster$h2o.clusters.summary$num.waters))
  num.h2o.clusters <- nrow(h2o.cluster$h2o.clusters.summary)
  conservation.mu <- mean(h2o.cluster$h2o.clusters.summary$num.waters)
  clustering.info <- h2o.cluster$clustering.info

  ##_ get the distance matrix size -----
  size.pairwise.dist <- unlist(strsplit(clustering.info$size.pairwise.dist, " "))
  size.pairwise.dist.num <- as.numeric(size.pairwise.dist[1])
  size.pairwise.dist.units <- size.pairwise.dist[2]

  ##_ calculate the percentages -----
  pct.conserved.passed <- h2o.cluster$h2o.clusters.summary$pct.conserved
  num.pct.lt0.5        <- sum(pct.conserved.passed  < 50)
  num.pct.gte0.5_lt0.7 <- sum(pct.conserved.passed >= 50 & pct.conserved.passed < 70)
  num.pct.gte0.7_lt0.8 <- sum(pct.conserved.passed >= 70 & pct.conserved.passed < 80)
  num.pct.gte0.8_lt0.9 <- sum(pct.conserved.passed >= 80 & pct.conserved.passed < 90)
  num.pct.gte0.9_lt1.0 <- sum(pct.conserved.passed >= 90 & pct.conserved.passed < 100)
  num.pct.gte.100 <- sum(pct.conserved.passed >= 100)

  ##_ combine the values -----
  h2o.cluster.summary.values <- c(num.pdbs.got.h2o,
                                  num.h2o.inital,
                                  num.h2o.used,
                                  num.h2o.clusters,
                                  conservation.mu,
                                  num.pct.lt0.5,
                                  num.pct.gte0.5_lt0.7,
                                  num.pct.gte0.7_lt0.8,
                                  num.pct.gte0.8_lt0.9,
                                  num.pct.gte0.9_lt1.0,
                                  num.pct.gte.100,
                                  as.numeric(clustering.info$num.pairwise.dist),
                                  size.pairwise.dist.num,
                                  as.numeric(clustering.info$time.pairwise.dist),
                                  as.numeric(clustering.info$time.cluster.dist) )

  ##_ combine the percentages -----
  h2o.cluster.summary.pct <- c(NA,
                               NA,
                               NA,
                               NA,
                               NA,
                               num.pct.lt0.5/num.h2o.clusters,
                               num.pct.gte0.5_lt0.7/num.h2o.clusters,
                               num.pct.gte0.7_lt0.8/num.h2o.clusters,
                               num.pct.gte0.8_lt0.9/num.h2o.clusters,
                               num.pct.gte0.9_lt1.0/num.h2o.clusters,
                               num.pct.gte.100/num.h2o.clusters,
                               NA,
                               NA,
                               NA,
                               NA)

  ##_ construct the data.frame of values and add rownames -----
  h2o.cluster.summary <- data.frame(values = h2o.cluster.summary.values,
                                    percentages = h2o.cluster.summary.pct * 100,
                                    stringsAsFactors = FALSE)
  rownames(h2o.cluster.summary) <- c("Number of structures",
                                     "Number of initial waters",
                                     "Number of used waters",
                                     "Number of water clusters",
                                     "Average conservation",
                                     "<50%",
                                     "50-69%",
                                     "70-79%",
                                     "80-89%",
                                     "90-99%",
                                     "100%",
                                     "Number of pairwise distances",
                                     paste("Memory size of pairwise distances (",
                                           size.pairwise.dist.units, ")", sep = ""),
                                     "Pairwise distance calc time (s)",
                                     "Cluster centroid distance calcs time (s)")

  ##_ return the summary statistics -----
  return(h2o.cluster.summary)

}



