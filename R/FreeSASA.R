## FreeSASA.R
##
## oct-13-2016 (exe) created
## nov-26-2016 (exe) updated documentation
## feb-06-2017 (exe) added examples
## feb-06-2017 (exe) updated documentation
## feb-09-2017 (exe) updated documentation (again...)
## mar-06-2017 (exe) updated for freeSASA 2.0
## jul-25-2017 (exe) updated documentation
## jul-31-2017 (exe) updated FreeSASA.diff() and FreeSASAcheck() documentation
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


## FreeSASA.diff docs ----------------------------------------------------------
#' @title Atomic SASA difference of hydrated PDB via FreeSASA
#' @description Calculates the atomic solvent accessible surface area (SASA) of
#'   the provided PDB (protein structure) using the FreeSASA application
#'   ([website](http://freesasa.github.io)).
#' @details The purpose of this function is to calculate and return the
#'   calculated atomic SASA for the provided PDB (protein structure) and the
#'   SASA of the protein when including the hydrating waters.
#'
#'   Several of the [FreeSASA](http://freesasa.github.io) options are set
#'   and _**NOT**_ user changeable. Specifically, no log information
#'   is returned; the `-L` ; the number of slices per atom is set to the
#'   [FreeSASA](http://freesasa.github.io) default of 20 (Lee & Richards
#'   algorithm); each [FreeSASA](http://freesasa.github.io) calculation
#'   uses four (4) threads; and the ProtOr atomic radii are used.
#'
#'   It might be too late if you are reading this, but it is
#'   _**strongly**_ encouraged to run [FreeSASAcheck()] to
#'   check if the [FreeSASA](http://freesasa.github.io) application is
#'   correctly installed.
#'
#' @param atoms.oi PDB structure read into R by [bio3d::read.pdb()];
#'   the [base::data.frame()] of `pdb$atom`
#' @param probeRadius Numerical values indicating the probe radius in Angstroms
#'   for the [FreeSASA](http://freesasa.github.io) application; default:
#'   1.4
#'
#' @return A PDB list with [FreeSASA](http://freesasa.github.io) (ProtOr)
#'   atomic radii placed in the _occupancy (o)_ column and SASA values
#'   calculated using the Lee & Richards method in the _b-value (b)_
#'   column.
#'
#' @export
#'
#' @import bio3d
#'
#' @examples
#' \dontrun{
#'   SASA.diff <- FreeSASA.diff(atoms.oi = thrombin.1hai$atom,
#'                              probeRadius = 1.4)
#'   head(SASA.diff)
#'   #   uniq.atom.ids SASA.prot SASA.hetatm SASA.lost
#'   # 1   THR_1_L_N_1      0.00        0.00         0
#'   # 2  THR_1_L_CA_2      0.00        0.00         0
#'   # 3   THR_1_L_C_3      0.00        0.00         0
#'   # 4   THR_1_L_O_4      0.00        0.00         0
#'   # 5  THR_1_L_CB_5      1.92        1.92         0
#'   # 6 THR_1_L_OG1_6     11.25       11.25         0
#'   #
#'   stem(SASA.diff$SASA.lost)
#'   #
#'   # The decimal point is at the |
#'   #
#'   #  0 | 00000000000000000000000000000000000000000000000000000000000000000000+1721
#'   #  2 | 00000000000001111111111122222222223333333333444444445555555555566666+88
#'   #  4 | 00001111111122222333333333334444444455555566677777777788899999000000+23
#'   #  6 | 00000111111222222222223333444455556666677778889900001112222333334455
#'   #  8 | 001122222333344566777888999001111222223444556667899
#'   # 10 | 000000001233445567888990001333344567899999
#'   # 12 | 00001233334446678800111233344788889
#'   # 14 | 00134448881223589
#'   # 16 | 014466389
#'   # 18 | 0945578888
#'   # 20 | 22347702
#'   # 22 | 246999
#'   # 24 | 09457
#'   # 26 | 44
#'   # 28 | 6
#'   # 30 |
#'   # 32 |
#'   # 34 | 9
#'   # 36 |
#'   # 38 |
#'   # 40 |
#'   # 42 |
#'   # 44 |
#'   # 46 | 9
#'   #
#' }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "Hydrophilicity Evaluation"
#'
#' @references
#'   **ProtOr (Protein-Organic) atomic radii:** \cr
#'   Jerry Tsai, Robin Taylor, Cyrus Chothia, and Mark Gerstein. The packing
#'   density in proteins: standard radii and volumes. _J Mol Biol_, 1999,
#'   **290** (_1_), pp 253-266.
#'   [DOI: 10.1006/jmbi.1999.2829](http://doi.org/10.1006/jmbi.1999.2829)
#'   [PMID: 10388571](https://www.ncbi.nlm.nih.gov/pubmed/10388571)
#'
#'   **SASA calculation method:** \cr
#'   B Lee, FM Richards. The interpretation of protein structures: estimation
#'   of static accessibility. _J Mol Biol_, 1971, **55** (_3_),
#'   pp 379-400.
#'   [DOI: 10.1016/0022-2836(71)90324-X](http://doi.org/10.1016/0022-2836(71)90324-X)
#'
#'   **FreeSASA application:** \cr
#'   Simon Mitternacht. FreeSASA: An open source C library for solvent
#'   accessible surface area calculations \[version 1; referees: 2 approved\].
#'   _F1000Research_, 2016, **5**:189
#'   [DOI: 10.12688/f1000research.7931.1](http://doi.org/10.12688/f1000research.7931.1)
#'   [PMCID: PMC4776673](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4776673/)
#'   [FreeSASA](http://freesasa.github.io)
#'
#'
FreeSASA.diff <- function(atoms.oi, probeRadius = 1.4) {

  ## FILES ---------------------------------------------------------------------
  ##----- temporary files
  infile.pdb <- tempfile()
  outfile.sasa.prot <- tempfile()
  outfile.sasa.hetatm <- tempfile()

  ##----- input file (PDB)
  write.basic.pdb(file = infile.pdb, atoms.oi)


  ## FreeSASA OUTPUT / RESULTS OPTIONS -----------------------------------------
  ##----- SASA for each atom
  # sasa.prot <- paste("--B-value-file=", outfile.sasa.prot, sep = "")
  # sasa.hetatm <- paste("--B-value-file=", outfile.sasa.hetatm, sep = "")
  sasa.prot   <- paste("--output=", outfile.sasa.prot, sep = "")
  sasa.hetatm <- paste("--output=", outfile.sasa.hetatm, sep = "")


  ## THE COMMAND AND EXECUTION -------------------------------------------------
  ##----- build the commands
  ##--- command for protein only
  cmd.prot <- paste("freesasa --lee-richards -p ", probeRadius,
                    " -n 20 -t 4 --format=pdb --radii=protor -w",
                    sasa.prot, infile.pdb, sep = " ")
  ##--- command for protein + het-atoms
  cmd.hetatm <- paste("freesasa --lee-richards -p ", probeRadius,
                      " -n 20 -t 4 --format=pdb --radii=protor -H -w",
                      sasa.hetatm, infile.pdb, sep = " ")

  ##----- execute the FreeSASA command
  os.type <- .Platform$OS.type
  if ( os.type == "windows" ) {
    FreeSASA.status.prot <- shell(cmd.prot,
                                  ignore.stderr = FALSE,
                                  ignore.stdout = FALSE)
    FreeSASA.status.hetatm <- shell(cmd.hetatm,
                                    ignore.stderr = FALSE,
                                    ignore.stdout = FALSE)
    FreeSASA.status <- FreeSASA.status.prot + FreeSASA.status.hetatm
  } else {
    FreeSASA.status.prot <- system(cmd.prot,
                                   ignore.stderr = FALSE,
                                   ignore.stdout = FALSE)
    FreeSASA.status.hetatm <- system(cmd.hetatm,
                                     ignore.stderr = FALSE,
                                     ignore.stdout = FALSE)
    FreeSASA.status <- FreeSASA.status.prot + FreeSASA.status.hetatm
  }

  if ( FreeSASA.status != 0 ) {
    mess <- paste("\nWhoa! Not cool. The FreeSASA application is not ",
                  "compiled/installed correctly. Please make sure FreeSASA ",
                  "compiled cleanly and the freesasa executable is in",
                  "your path. Please visit http://freesasa.github.io ",
                  "and following the posted instructions. Once completed ",
                  "run the FreeSASAcheck function.", sep = "")
    stop(mess)
  }


  ## READ IN THE SASA RESULTS --------------------------------------------------
  ##----- read the FreeSASA results
  SASA.results.prot <- bio3d::read.pdb(file = outfile.sasa.prot)$atom
  SASA.results.hetatm <- bio3d::read.pdb(file = outfile.sasa.hetatm)$atom


  ## ORGANIZE THE SASA RESULTS -------------------------------------------------
  ##----- create the SASA data.frame
  cols.oi <- c("resid", "resno", "chain", "elety", "eleno")
  SASA.prot.ids <- UniqueAtomHashes(atoms.oi = SASA.results.prot,
                                    cols.oi = cols.oi, separator = "_")
  SASA.hetatm.ids <- UniqueAtomHashes(atoms.oi = SASA.results.hetatm,
                                      cols.oi = cols.oi, separator = "_")

  SASA.results.prot <- cbind(SASA.results.prot,
                             uniq.atom.ids = SASA.prot.ids,
                             stringsAsFactors = FALSE)

  SASA.results.hetatm <- cbind(SASA.results.hetatm,
                               uniq.atom.ids = SASA.hetatm.ids,
                               stringsAsFactors = FALSE)

  ##----- merge the results
  SASA.merge <- merge(x = SASA.results.prot[, c("o", "b", "uniq.atom.ids")],
                      y = SASA.results.hetatm[, c("o", "b", "uniq.atom.ids")],
                      by.x = "uniq.atom.ids", by.y = "uniq.atom.ids",
                      suffixes = c(".prot", ".hetatm"),
                      all = FALSE, sort = FALSE)

  SASA.lost <- SASA.merge$b.prot - SASA.merge$b.hetatm

  SASA.results <- data.frame(uniq.atom.ids = SASA.merge$uniq.atom.ids,
                             SASA.prot = SASA.merge$b.prot,
                             SASA.hetatm = SASA.merge$b.hetatm,
                             SASA.lost = SASA.lost,
                             stringsAsFactors = FALSE)


  ## FINISH --------------------------------------------------------------------
  ##----- cleanup the temporary files
  unlink(c(infile.pdb, outfile.sasa.prot, outfile.sasa.hetatm))

  ##----- return the SASA results
  return(SASA.results)
}





## FreeSASAcheck docs ----------------------------------------------------------
#' @title FreeSASA Check
#' @description Determines if FreeSASA is (correctly) installed.
#' @details Because [FreeSASA](http://freesasa.github.io) is
#'   _**NOT**_ included with [vanddraabe] it
#'   is important to ensure the application has been installed and was correctly
#'   compiled.
#'
#' @return
#'   When [FreeSASA](http://freesasa.github.io) is correctly installed
#'   the current version and citation are returned to the user: \cr
#'   ```
#'   FreeSASA 2.0
#'   License: MIT <http://opensource.org/licenses/MIT>
#'   If you use this program for research, please cite:
#'     Simon Mitternacht (2016) FreeSASA: An open source C
#'     library for solvent accessible surface area calculations.
#'   F1000Research 5:189.
#'   ```
#'
#'   When [FreeSASA](http://freesasa.github.io) is _**NOT**_
#'   correctly installed the following are returned to the user:
#'   ```
#'   Error in FreeSASAcheck() :
#'   Uh-oh!!
#'   Please make sure FreeSASA is correctly installed! Please visit
#'   (http://freesasa.github.io) for
#'   instructions specific to your operating system.
#'   ```
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Result for correct installation
#'   FreeSASAcheck()
#'   # FreeSASA 2.0
#'   # License: MIT <http://opensource.org/licenses/MIT>
#'   # If you use this program for research, please cite:
#'   #   Simon Mitternacht (2016) FreeSASA: An open source C
#'   #   library for solvent accessible surface area calculations.
#'   #   F1000Research 5:189.
#'   #
#'   # Report bugs to <https://github.com/mittinatten/freesasa/issues>
#'   # Home page: <http://freesasa.github.io>
#'   #
#'   # Congratulations! FreeSASA is correctly installed!
#'   #
#'   # Result for incorrect installation
#'   FreeSASAcheck()
#'   # Error:
#'   # Uh-oh!!
#'   # Please make sure FreeSASA is correctly installed. Please visit
#'   # http://freesasa.github.io for instructions specific to your operating
#'   # system.
#' }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "vanddraabe utilities"
#'
#' @references
#'   Simon Mitternacht. FreeSASA: An open source C library for solvent
#'   accessible surface area calculations \[version 1; referees: 2 approved\].
#'   _F1000Research_, 2016, **5**:189
#'   [DOI: 10.12688/f1000research.7931.1](http://doi.org/10.12688/f1000research.7931.1)
#'   [PMCID: PMC4776673](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4776673/)
#'   [FreeSASA](http://freesasa.github.io)
#'
FreeSASAcheck <- function() {

  ##----- check if the program is available
  os.type <- .Platform$OS.type
  if ( os.type == "windows" ) {
    FreeSASA.status <- shell("freesasa --version",
                             ignore.stderr = FALSE,
                             ignore.stdout = FALSE)
  } else {
    FreeSASA.status <- system("freesasa --version",
                              ignore.stderr = FALSE,
                              ignore.stdout = FALSE)
  }


  ##----- inform user of FreeSASA's status
  if ( FreeSASA.status == 0 ) {
    message("\nCongratulations! FreeSASA is correctly installed!")
  } else {
    mess <- paste("\nUh-oh!!\n",
                  "Please make sure FreeSASA is correctly installed. Please ",
                  "visit http://freesasa.github.io for instructions specific ",
                  "to your operating system.", sep = "")
    stop(mess)
  }

}
