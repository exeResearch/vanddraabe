## BoundWaterEnvironments.R
##
## apr-03-2016 (exe) created
## jan-20-2017 (exe) corrected formatting based on lintr
## feb-06-2017 (exe) added examples
## feb-07-2017 (exe) updated documentation
## mar-29-2017 (exe) update HydrophilicityTable column names in
##                   BoundWaterEnvironment()
## apr-25-2017 (exe) created sub-fxns for BoundWaterEnvironment()
## apr-25-2017 (exe) created calcNearbyHydrationFraction() fxn
## apr-26-2017 (exe) created calcNumHydrogenBonds() fxn
## apr-26-2017 (exe) created BoundWaterEnvironment.quality() fxn
## apr-26-2017 (exe) created BoundWaterEnvironment.interact() fxn
## jul-25-2017 (exe) updated documentation
## jul-28-2017 (exe) updated NormalizedBvalue and Mobility documentation;
##                   changed Bvalue parameter to Bvalues to match code
## jul-31-2017 (exe) added calcBvalue() example
## aug-10-2017 (exe) added @importFrom stats ... and utils ...
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


## Normalize B-values docs -----------------------------------------------------
#' @title B-value Normalization
#' @description Calculate the normalized B-value values of waters for a
#'   structure.
#' @details The normalized B-value values are the number of standard deviations
#'   from the mean for the water oxygens' B-values within the structure of
#'   interest.
#'
#'   The B-value normalization exclusion value is user defined within the main
#'   [ConservedWaters()] function but has a default value of 1.0.
#'
#' @param Bvalues B-value values
#'
#' @return Vector of normalized and unitless B-value values.
#'
#' @examples
#'   set.seed(13)
#'   Bvalues <- sample(thrombin.1hai$atom$b, 10)
#'   Bvalues
#'   # [1] 32.53 22.36 24.91 42.11 36.60
#'   #     54.66 37.71 14.93 27.65 17.84
#'   NormalizedBvalue(Bvalues)
#'   # [1]  0.1158 -0.7252 -0.5143  0.9080  0.4523
#'   #      1.9457  0.5441 -1.3396 -0.2878 -1.0990
#'
#' @export
#'
#' @importFrom stats sd
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Bound Water Environment"
#'
#' @references
#'   Oliviero Carugo. Correlation between occupancy and B value of
#'   water molecules in protein crystal structures. _Protein Engineering_,
#'   1999, **12** (_12_), pp 1021-1024.
#'   [DOI: 10.1093/protein/12.12.1021](http://doi.org/10.1093/protein/12.12.1021)
#'   [PMID: 10611392](http://www.ncbi.nlm.nih.gov/pubmed/10611392)
#'
NormalizedBvalue <- function (Bvalues) {

  Bvalues.mu <- mean(Bvalues)
  Bvalues.sd <- sd(Bvalues)
  Bvalues.norm <- (Bvalues - Bvalues.mu) / Bvalues.sd

  return(Bvalues.norm)

}


## Mobility docs ---------------------------------------------------------------
#' @title Water Molecule Mobility
#' @description Calculate the mobility values of waters for a structure.
#' @details The mobility of waters within a structure is normalization method to
#'   identify the amount of variance an atom has within a structure. In the case
#'   of waters, identified by an oxygen atom without hydrogen atoms, a
#'   water-oxygen atom with a mobility value of 0 is considered rigid and does
#'   not possess variance. The average mobility within a structure has value of
#'   1 while an atom's mobility value of x is considered x-times as mobile as an
#'   average atom.
#'
#'   \deqn{Mobility = \frac{\frac{B-value}{\mu_{B-value}}}{\frac{Occupancy}{\mu_{Occupancy}}}}{Mobility = (B-value/mu B-value)/(Occupancy/mu Occupancy)}
#'
#'   Mobility is calculated using the B-value and occupancy values; these
#'   values are a byproduct of solving the 3D molecular structure from electron
#'   density maps. The mobility values allows us to compare atomic mobility
#'   between molecular structures solved using different structural refinement
#'   methods. Atoms, in this instance water-oxygens, with a mobility value
#'   greater than 2.0 are removed from analysis.
#'
#'   The mobility exclusion value is user defined within the main
#'   [ConservedWaters()] function but has a default value of 2.0.
#'
#' @param Bvalues B-value values from the imported PDB file(s)
#' @param occupancy Occumpancy values from the imported PDB file(s)
#'
#' @return Vector of mobility values; unitless.
#'
#' @examples
#'   set.seed(13)
#'   sample.idc <- sample(1:nrow(thrombin.1hai$atom), 10)
#'   Bvalues <- thrombin.1hai$atom[sample.idc, "b"]
#'   Bvalues
#'   # [1] 32.53 22.36 24.91 42.11 36.60
#'   #     54.66 37.71 14.93 27.65 17.84
#'   occupancy <- thrombin.1hai$atom[sample.idc, "o"]
#'   occupancy
#'   # [1] 1.00 1.00 1.00 1.00 0.79
#'   #     1.00 1.00 1.00 1.00 1.00
#'   Mobility(Bvalues, occupancy)
#'   # [1] 1.0230 0.7032 0.7834 1.3243 1.4570
#'   #     1.7190 1.1859 0.4695 0.8696 0.5610
#'
#' @export
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Bound Water Environment"
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
Mobility <- function(Bvalues, occupancy) {

  Bvalues.mu <- mean(Bvalues)
  Occupancy.mu <- mean(occupancy)

  mobility.score <- (Bvalues / Bvalues.mu) / (occupancy / Occupancy.mu)

  return(mobility.score)

}



## calcBvalue docs -------------------------------------------------------------
#' @title Calculate B-value
#' @description Calculate the B-value for an atom.
#' @details The B-value (aka B-factor) is calcualted from the rmsf from a
#'   collection of atoms. The rmsf is calculated using [bio3d::rmsf()].
#'
#'   \deqn{B-value = rmsf^{2} * 8 * {pi}^{2}}{B-value = rmsf^2 * 8 * pi^2}
#'
#'   The calculated B-values are returned within the [BoundWaterEnvironment()]
#'   results and used to define the size of conserved waters for the depiction
#'   of MDS conserved waters.
#'
#' @param rmsfValue rmsf value calculated by [bio3d::rmsf()]
#'
#' @return B-value (aka B-factor) in Angstroms^2^
#'
#' @examples
#'   calcBvalue(rmsfValue=0.25)
#'   # [1] 4.935
#'   calcBvalue(rmsfValue=0.50)
#'   # [1] 19.74
#'   calcBvalue(rmsfValue=0.75)
#'   # [1] 44.41
#'   calcBvalue(rmsfValue=1.0)
#'   # [1] 78.96
#'   calcBvalue(rmsfValue=1.25)
#'   # [1] 123.4
#'
#' @export
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Bound Water Environment"
#'
#' @references
#'   Eaton E Lattman & Patrick J Loll. _Protein Crystallography: A Concise
#'   Guide_. Baltimore, Maryland, USA: The Johns Hopkins University Press, 2008.
#'   QP551.L345 2008. ISBN: 978-0-8018-8808-3
#'   [website](https://jhupbooks.press.jhu.edu/content/protein-crystallography)
#'
calcBvalue <- function(rmsfValue) {

  Bvalue <- rmsfValue * rmsfValue * 8 * pi * pi

  ##----- return B-value
  return(Bvalue)

}


## calculate nearby hydration fraction docs ------------------------------------
#' @title Calculate Nearby Atom Hydration Fraction
#' @description Calculate the mobility values of waters for a structure.
#' @details The summation, mean, and standard deviation of the hydrophilicity
#'   fraction for the protein atoms within the user specified distance for the
#'   [BoundWaterEnvironment()] function are calculated and returned.
#'
#' @param names.res.nearby.atoms string of residue-atom name for nearby atoms
#'
#' @return Hydrophilicity fraction sum, mean, and standard deviation.
#'
#' @export
#'
#' @importFrom stats sd
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Bound Water Environment"
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
calcNearbyHydrationFraction <- function(names.res.nearby.atoms) {

  if ( length(names.res.nearby.atoms) > 0 ) {
    ##----- get the hydration fraction from the lookup table
    nearby.hydro.idc <- match(names.res.nearby.atoms,
                              HydrophilicityTable$residueAtomType)
    nearby.hydro.values <- HydrophilicityTable$hydratFraction[nearby.hydro.idc]

    ##----- calculate the summation, mean, and standard deviation
    nearby.hydro.values.sum <- sum(nearby.hydro.values, na.rm = TRUE)
    nearby.hydro.values.mu  <- mean(nearby.hydro.values, na.rm = TRUE)
    nearby.hydro.values.sd  <- sd(nearby.hydro.values, na.rm = TRUE)
  } else {
    nearby.hydro.values.sum <-
      nearby.hydro.values.mu <-
      nearby.hydro.values.sd <- 0
  }

  ##----- return values
  list(ahp.sum = nearby.hydro.values.sum,
       ahp.mu = nearby.hydro.values.mu,
       ahp.sd = nearby.hydro.values.sd)

}


## calculate number of hydrogen bonds docs -------------------------------------
#' @title Calculate Number of Hydrogen Bonds
#' @description Calculate the number of hydrogen bonds.
#' @details The summation, mean, and standard deviation of the hydrophilicity
#'   fraction for the protein atoms within the user specified distance for the
#'   [BoundWaterEnvironment()] function are calculated and returned.
#'
#' @param distances between water atom of interest and the protein atoms, water
#'   oxygen atoms, or HETATMs
#' @param nearby.atoms.idc numeric vector of atom indices near water of interest
#' @param names.atoms names of atoms; _e.g._; c("CB", "CA", "N", "O", "CZ")
#' @param set.oi.idc numeric vector of indices for protein atoms, water oxygen
#'   atoms, or HETATMs
#'
#' @return Number of possible hydrogen bonds between the water of interest and
#'   the protein atoms within 3.5 Angstroms of the water.
#'
#' @examples
#'   \dontrun{
#'   distances <- PDB.1hai.h2o.prot.dists[3, ]
#'   nearby.atoms.idc <- Nearby(distances, set.idc = prot.idc, radius = 3.6)
#'   names.atoms <- PDB.1hai.aoi.clean$elety[prot.idc]
#'   calcNumHydrogenBonds(distances, nearby.atoms.idc, names.atoms,
#'     set.oi.idc = prot.idc)
#'   # [1] 4
#'   }
#'
#' @export
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Bound Water Environment"
#'
calcNumHydrogenBonds <- function(distances,
                                 nearby.atoms.idc,
                                 names.atoms,
                                 set.oi.idc) {

  ##----- determine the number of hydrogen bonds
  if ( length(nearby.atoms.idc) > 0 ){
    hbonds.idc <- Nearby(distances, set.idc = set.oi.idc, radius = 3.5)
    num.hbonds <- sum(!is.na(match(names.atoms[hbonds.idc], names.polar.atoms)))
  } else {
    num.hbonds <- 0
  }

  ##----- return value
  return(num.hbonds)

}


## BoundWaterEnvironment.interact docs -----------------------------------------
#' @title Bound Water Environment (interactions)
#' @description Various enviroment counts for bound waters.
#' @details For the heavy atoms near each water molecule (oxygen atom) the bound
#'   water environment is calculated. These values are defined in the **Return**
#'   section. The default radius distance is 3.6 Angstroms. While it is possible
#'   to define the radius to a value other than 3.6 this value is hardcoded into
#'   the [ConservedWaters()] function. This might change in future versions.
#'
#'   _**NOTE**_: This function is designed to work with [ConservedWaters()] via
#'   the [base::apply()] function processing rows (the `MARGIN = 1` option). For
#'   this reason it is **NOT** a public function. The [Nearby()] is specifically
#'   designed to work with this function.
#'
#' @param distances Matrix of atomic pairwise distances
#' @param set.oi.idc Indices of protein atoms; can also HETATMs if those are of
#'   interest
#' @param names.atoms Atom names from the PDB file in the PDB atomic naming
#'   convention.
#' @param names.res.atoms Atom names of the form "RES AT"; created by combining
#'   the residue and atom name while separating the two by a space. These do not
#'   need to be unique because these names will be used to lookup hydrophilicity
#'   values.
#' @param radius Distance in Angstroms between the atoms of interest; default:
#'   3.6 Angstroms
#'
#' @return
#'   A list of the bound water environment values for nearby heavy atoms.
#'   * **adn**: num of nearby heavy atoms
#'   * **ahp.sum**: sum of hydrodrophilicy values
#'   * **ahp.mu**: mean of hydrodrophilicy values
#'   * **ahp.sd**: standard deviation of hydrodrophilicy values
#'   * **hbonds**: number of possible hydrogen bonds
#'
#' @examples
#'   \dontrun{
#'   distances <- PDB.1hai.h2o.prot.dists[3, ]
#'   set.oi.idc <- prot.idc
#'   names.atoms <- PDB.1hai.aoi.clean$elety[prot.idc]
#'   names.res.atoms <- paste(PDB.1hai.aoi.clean$resid[prot.idc], names.atoms, sep =" ")
#'   BoundWaterEnvironment.interact(distances,
#'                                  set.oi.idc,
#'                                  names.atoms,
#'                                  names.res.atoms,
#'                                  radius = 3.6)
#'   # $adn
#'   # [1] 9
#'   #
#'   # $ahp.sum
#'   # [1] 2.001
#'   #
#'   # $ahp.mu
#'   # [1] 0.2223
#'   #
#'   # $ahp.sd
#'   # [1] 0.2229
#'   #
#'   # $hbonds
#'   # [1] 4
#'   }
#'
#' @export
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Bound Water Environment"
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
#'   Leslie A Kuhn, Craig A Swanson, Michael E Pique, John A Tainer,
#'   and Elizabeth D Getzof. Atomic and Residue Hydrophilicity in the Context of
#'   Folded Protein Structures. _PROTEINS: Structure, Function, and
#'   Genetics_, 1995, **2** (_4_), pp 536-547.
#'   [DOI: 10.1002/prot.340230408](http://doi.org/10.1002/prot.340230408)
#'   [PMID: 8749849](http://www.ncbi.nlm.nih.gov/pubmed/8749849)
#'
BoundWaterEnvironment.interact <- function(distances,
                                           set.oi.idc,
                                           names.atoms,
                                           names.res.atoms,
                                           radius = 3.6) {

  ##----- atomic density
  ##      protein atoms within user-defined distance (Angstroms)
  nearby.atoms.idc <- Nearby(distances, set.idc = set.oi.idc, radius = radius)
  num.nearby.atoms <- length(nearby.atoms.idc)
  names.res.nearby.atoms <- names.res.atoms[nearby.atoms.idc]

  ##----- local atomic hydrophilicity
  localHydratFract <- calcNearbyHydrationFraction(names.res.nearby.atoms)

  ##----- number of hydrogen bonds
  num.hbonds <- calcNumHydrogenBonds(distances,
                                     nearby.atoms.idc,
                                     names.atoms,
                                     set.oi.idc)

  ##----- return the results
  list(adn = num.nearby.atoms,
       ahp.sum = localHydratFract$ahp.sum,
       ahp.mu = localHydratFract$ahp.mu,
       ahp.sd = localHydratFract$ahp.sd,
       hbonds = num.hbonds)

}


## BoundWaterEnvironment.quality docs ------------------------------------------
#' @title Bound Water Environment (atomic quality)
#' @description Various enviroment counts for bound waters.
#' @details For the heavy atoms near each water molecule (oxygen atom) the bound
#'   water environment is calculated. These values are defined in the **Return**
#'   section. The default radius distance is 3.6 Angstroms. While it is possible
#'   to define the radius to a value other than 3.6 this value is hardcoded into
#'   the [ConservedWaters()] function. This might change in future versions.
#'
#'   _**NOTE**_: This function is designed to work with [ConservedWaters()] via
#'   the [base::apply()] function processing rows (the `MARGIN = 1` option). For
#'   this reason it is **NOT** a public function. The [Nearby()] is specifically
#'   designed to work with this function.
#'
#' @param distances Matrix of atomic pairwise distances
#' @param set.oi.idc Indices of atoms of interest; can be protein, water, or
#'   HETATMs if those are of interest
#' @param structure The protein structure of interest with its residue and atom
#'   names; X, Y, and Z coordinates; residue and atom numbers; and B-value,
#'   Normalized B-value, Occupancy, and Mobility values.
#' @param radius Distance in Angstroms between the atoms of interest; default:
#'   3.6 Angstroms
#'
#' @return
#'   A list of the bound water environment values for nearby heavy atoms.
#'   * **o.sum**: sum of occupancy values
#'   * **o.mu**: mean of occupancy values
#'   * **o.sd**: standard deviation of occupancy values
#'   * **b.exp.sum**: sum of experimental B-values
#'   * **b.exp.mu**: mean of experimental B-values
#'   * **b.exp.sd**: standard deviation of experimental B-values
#'   * **mobility.sum**: sum of mobility values
#'   * **mobility.mu**: mean of mobility values
#'   * **mobility.sd**: standard deviation of mobility values
#'   * **nBvalue.sum**: sum of normalized B-values
#'   * **nBvalue.mu**: mean of normalized B-values
#'   * **nBvalue.sd**: standard deviation of normalized B-values
#'
#' @examples
#'   \dontrun{
#'   distances <- PDB.1hai.h2o.prot.dists[3, ]
#'   set.oi.idc <- prot.idc
#'   structure <- PDB.1hai.aoi.clean
#'   BoundWaterEnvironment.quality(distances,
#'                                 set.oi.idc,
#'                                 structure,
#'                                 radius = 3.6)
#'   }
#'
#' @importFrom stats sd
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Bound Water Environment"
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
#'   Leslie A Kuhn, Craig A Swanson, Michael E Pique, John A Tainer,
#'   and Elizabeth D Getzof. Atomic and Residue Hydrophilicity in the Context of
#'   Folded Protein Structures. _PROTEINS: Structure, Function, and
#'   Genetics_, 1995, **2** (_4_), pp 536-547.
#'   [DOI: 10.1002/prot.340230408](http://doi.org/10.1002/prot.340230408)
#'   [PMID: 8749849](http://www.ncbi.nlm.nih.gov/pubmed/8749849)
#'
BoundWaterEnvironment.quality <- function(distances,
                                          set.oi.idc,
                                          structure,
                                          radius = 3.6) {

  ##_ atomic density -----
  ##  protein atoms within user-defined distance (Angstroms)
  nearby.atoms.idc <- Nearby(distances, set.idc = set.oi.idc, radius = radius)
  num.nearby.atoms <- length(nearby.atoms.idc)
  # names.res.nearby.atoms <- names.res.atoms[nearby.atoms.idc]

  ##_ calculate the quality -----
  if ( num.nearby.atoms > 0 ) {
    ##__ get the values -----
    b.values        <- structure$b[nearby.atoms.idc]
    o.values        <- structure$o[nearby.atoms.idc]
    mobility.values <- structure$mobility[nearby.atoms.idc]
    nBvalue.values  <- structure$nBvalue[nearby.atoms.idc]
    ##__ calculate the sum -----
    b.exp.sum    <- sum(b.values, na.rm = TRUE)
    o.sum        <- sum(o.values, na.rm = TRUE)
    mobility.sum <- sum(mobility.values, na.rm = TRUE)
    nBvalue.sum  <- sum(nBvalue.values, na.rm = TRUE)
    ##__ calculate the mean -----
    b.exp.mu    <- suppressWarnings( mean(b.values, na.rm = TRUE) )
    o.mu        <- suppressWarnings( mean(o.values, na.rm = TRUE) )
    mobility.mu <- suppressWarnings( mean(mobility.values, na.rm = TRUE) )
    nBvalue.mu  <- suppressWarnings( mean(nBvalue.values, na.rm = TRUE) )
    ##__ calculate the standard deviation -----
    b.exp.sd    <- suppressWarnings( sd(b.values, na.rm = TRUE) )
    o.sd        <- suppressWarnings( sd(o.values, na.rm = TRUE) )
    mobility.sd <- suppressWarnings( sd(mobility.values, na.rm = TRUE) )
    nBvalue.sd  <- suppressWarnings( sd(nBvalue.values, na.rm = TRUE) )
  } else {
    b.exp.sum <- o.sum <- mobility.sum <- nBvalue.sum <- 0
    b.exp.mu  <- o.mu  <- mobility.mu  <- nBvalue.mu  <- 0
    b.exp.sd  <- o.sd  <- mobility.sd  <- nBvalue.sd  <- 0
  }

  ##_ return the results -----
  list(o.sum = o.sum,
       o.mu = o.mu,
       o.sd = o.sd,
       b.exp.sum = b.exp.sum,
       b.exp.mu = b.exp.mu,
       b.exp.sd = b.exp.sd,
       mobility.sum = mobility.sum,
       mobility.mu = mobility.mu,
       mobility.sd = mobility.sd,
       nBvalue.sum = nBvalue.sum,
       nBvalue.mu = nBvalue.mu,
       nBvalue.sd = nBvalue.sd)

}


## Bound Water Environment docs ------------------------------------------------
#' @title Bound Water Environment
#' @description Various enviroment counts for bound waters.
#' @details For the heavy atoms near each water molecule (oxygen atom) the bound
#'   water environment is calculated. These values are defined in the **Return**
#'   section. The default radius distance is 3.6 Angstroms. While it is possible
#'   to define the radius to a value other than 3.6 this value is hardcoded into
#'   the [ConservedWaters()] function. This might change in future versions.
#'
#'   _**NOTE**_: This function is designed to work with [ConservedWaters()] via
#'   the [base::apply()] function processing rows (the `MARGIN = 1` option). For
#'   this reason it is **NOT** a public function. The [Nearby()] is specifically
#'   designed to work with this function.
#'
#' @param distances Matrix of atomic pairwise distances
#' @param set.oi.idc Indices of atoms of interest; can be protein, water, or
#'   HETATMs if those are of interest
#' @param names.atoms Atom names for the atoms of interest. Valid atom names are
#'   provided in the [names.backbone.atoms()] and [names.sidechain.atoms()]
#'   functions; _e.g._; "C", "O", "CB", "OG1", "CG2", "N"
#' @param names.res.atoms Residue and atom names of interest. Valid residue-atom
#'   names are provided in the [names.res.AtomTypes()] function; _e.g._; "THR
#'   C", "THR O", "THR CB", "THR OG1"
#' @param structure The protein structure of interest with its residue and atom
#'   names; X, Y, and Z coordinates; residue and atom numbers; and B-value,
#'   Normalized B-value, Occupancy, and Mobility values.
#' @param radius Distance in Angstroms between the atoms of interest; default:
#'   3.6 Angstroms
#'
#' @return
#'   A list of the bound water environment values for nearby heavy atoms.
#'   * **adn**: num of nearby heavy atoms
#'   * **ahp.sum**: sum of hydrodrophilicy values
#'   * **ahp.mu**: mean of hydrodrophilicy values
#'   * **ahp.sd**: standard deviation of hydrodrophilicy values
#'   * **hbonds**: number of possible hydrogen bonds
#'   * **o.sum**: sum of occupancy values
#'   * **o.mu**: mean of occupancy values
#'   * **o.sd**: standard deviation of occupancy values
#'   * **b.exp.sum**: sum of experimental B-values
#'   * **b.exp.mu**: mean of experimental B-values
#'   * **b.exp.sd**: standard deviation of experimental B-values
#'   * **mobility.sum**: sum of mobility values
#'   * **mobility.mu**: mean of mobility values
#'   * **mobility.sd**: standard deviation of mobility values
#'   * **nBvalue.sum**: sum of normalized Bvalues
#'   * **nBvalue.mu**: mean of normalized Bvalues
#'   * **nBvalue.sd**: standard deviation of normalized Bvalues
#'
#' @examples
#'   \dontrun{
#'   distances <- PDB.1hai.h2o.prot.dists[3, ]
#'   set.oi.idc <- prot.idc
#'   names.atoms <- PDB.1hai.aoi.clean$elety[prot.idc]
#'   names.res.atoms <- paste(PDB.1hai.aoi.clean$resid[prot.idc], names.atoms, sep =" ")
#'   structure <- PDB.1hai.aoi.clean
#'   BoundWaterEnvironment(distances,
#'                         set.oi.idc,
#'                         names.atoms,
#'                         names.res.atoms,
#'                         structure,
#'                         radius = 3.6)
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Bound Water Environment"
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
#'   Leslie A Kuhn, Craig A Swanson, Michael E Pique, John A Tainer,
#'   and Elizabeth D Getzof. Atomic and Residue Hydrophilicity in the Context of
#'   Folded Protein Structures. _PROTEINS: Structure, Function, and
#'   Genetics_, 1995, **2** (_4_), pp 536-547.
#'   [DOI: 10.1002/prot.340230408](http://doi.org/10.1002/prot.340230408)
#'   [PMID: 8749849](http://www.ncbi.nlm.nih.gov/pubmed/8749849)
#'
BoundWaterEnvironment <- function(distances,
                                  set.oi.idc,
                                  names.atoms,
                                  names.res.atoms,
                                  structure,
                                  radius = 3.6) {

  ##_ calculate the interactions -----
  interaction.values <- BoundWaterEnvironment.interact(distances,
                                                       set.oi.idc,
                                                       names.atoms,
                                                       names.res.atoms,
                                                       radius = radius)

  ##_ calculate the quality -----
  quality.values <- BoundWaterEnvironment.quality(distances,
                                                  set.oi.idc,
                                                  structure,
                                                  radius = radius)

  ##_ return the results -----
  list(adn = interaction.values$adn,
       ahp.sum = interaction.values$ahp.sum,
       ahp.mu = interaction.values$ahp.mu,
       ahp.sd = interaction.values$ahp.sd,
       hbonds = interaction.values$hbonds,
       o.sum = quality.values$o.sum,
       o.mu = quality.values$o.mu,
       o.sd = quality.values$o.sd,
       b.exp.sum = quality.values$b.exp.sum,
       b.exp.mu = quality.values$b.exp.mu,
       b.exp.sd = quality.values$b.exp.sd,
       mobility.sum = quality.values$mobility.sum,
       mobility.mu = quality.values$mobility.mu,
       mobility.sd = quality.values$mobility.sd,
       nBvalue.sum = quality.values$nBvalue.sum,
       nBvalue.mu = quality.values$nBvalue.mu,
       nBvalue.sd = quality.values$nBvalue.sd)

}
