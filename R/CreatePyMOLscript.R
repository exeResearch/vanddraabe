## CreatePyMOLscript.R
##
## aug-01-2016 (exe) created
## aug-05-2016 (exe) working
## sep-13-2016 (exe) set pocket residues to a single color
## oct-07-2016 (exe) added PyMOL file location message
## jan-20-2017 (exe) corrected formatting based on lintr
## feb-06-2017 (exe) split CreatePyMOLfile into script creation and file writing
##                   functions to aid unit testing
## feb-28-2017 (exe) updated documentation
## apr-11-2017 (exe) stand-alone function
## jul-25-2017 (exe) updated documentation
## aug-02-2017 (exe) streamlined input parameters
## aug-09-2017 (exe) added logic for presence/absence of ligands (primarily for MDS analysis)
## aug-09-2017 (exe) added vanddraabe and R version to PyMOL script preamble
## aug-09-2017 (exe) added @importFrom utils ...
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



## create PyMOL script file docs -----------------------------------------------
#' @title Create PyMOL Script File
#' @description Create PyMOL script file to visualize conserved waters
#' @details The ability to visualize the conserved waters is important and their
#'   surroundings is when exploring conserved water results.
#'
#'   Conserved waters within 6 Angstroms of the PyMOL identified ligands are
#'   displayed. The conserved waters are colored based on their percent
#'   conservation range using the same color scheme as the Percent Conservation
#'   plot. Waters conserved less than 50\% are colored light grey, 50-69\% are
#'   red, 70-79\% are dark red, 80-89\% are light blue, 90-99\% are medium blue,
#'   and 100\% are dark blue. The conserved waters are labeled using their
#'   ranking based on percent conseration.
#'
#'   This function creates _**two**_ PyMOL script files; one with a
#'   black background and another with a white background. The color of the
#'   pocket residues is changed based on the background. The pocket residues are
#'   colored light-grey for the black background and dark-grey for the white
#'   background. The ligand is assigned the user-defined color for both
#'   representations. Pocket residues -- and associated molecular surface -- are
#'   defined as those within 5 Angstroms of the conserved waters. The depicted
#'   cartoon representation is for residues within 15 Angstroms of the
#'   ligand(s).
#'
#'   The potential hydrogen bonds are depicted between:
#'   - conserved waters and ligand: orange dashed line
#'   - conserved waters and protein: green dashed line
#'   - conserved waters: blue dashed line
#'
#' @param conservedWaters.data The `h2o.clusters.summary` data.frame from the
#'   [ConservedWaters()] function containing the `nBvalue.mu`
#'   information. This data.frame is found within the `h2o.cluster.passed`
#'   and `h2o.cluster.all`
#' @param passed.waters Logical indicator to plot results for waters **passing**
#'   [Mobility()] and [NormalizedBvalue()] _**OR**_ using **all** waters within
#'   the `PDB` files.
#' @param PDBid.ref name for reference structure in PyMOL; _e.g._,
#'   `"1hai"`
#' @param LigResname.ref PDB residue code for reference ligand; _e.g._, `"0g6"`
#' @param hbond The minimum distance between hydrogen bond acceptor and donor;
#'   default: `3.75`
#' @param lig.carbon.color One of the ten pre-defined carbon-color options using
#'   PyMOL's `util.cbaX` command. The `X` represents the user defined
#'   color of carbon atoms. `X` can be `g`: green; `c`: cyan;
#'   `m`: magenta; `y`: yellow; `s`: salmon; `w`: grey;
#'   `b`: slate; `o`: orange; `p`: purple; and `k`: pink;
#'   default: `cyan`
#' @param filename Prefix for the PyMOL script files. It is probably best to use
#'   the initial portion of the conserved waters PDB filename; _e.g._,
#'   `"thrombin10"`
#'
#' @export
#'
#' @import bio3d
#' @importFrom utils packageVersion
#'
#' @examples
#'   \dontrun{
#'   current.time <- Sys.time()
#'   CreatePyMOLfile(PDBid.ref = "Thrombin_initial10_alignedGood/1hai_aligned_pruned.pdb",
#'                   PDBid.ref = "1hai",
#'                   LigResname.ref = "0g6",
#'                   conserved.waters = "Thrombin_initial10_ConservedWaters_PASSED_mar292017_1535.pdb",
#'                   hbond = 3.75,
#'                   lig.carbon.color = "cyan",
#'                   filename = "thrombin10_ConservedWaters_PASSED")
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family Visualization
#'
CreatePyMOLscript <- function(conservedWaters.data,
                              passed.waters = TRUE,
                              PDBid.ref = "1hai",
                              LigResname.ref = NULL,
                              hbond = 3.75,
                              lig.carbon.color = "cyan",
                              filename = "thrombin10") {

  ## get the PASSED or ALL waters plot data -----------------------------------
  ##- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  is.MDS <- conservedWaters.data$MDS
  if (is.MDS == TRUE) {
    passed.waters <- FALSE
  }
  if (passed.waters == TRUE) {
    PDB.conserved.waters <- conservedWaters.data$filenames.used$conserved.passed
    PyMOL.filename <- paste0(filename, "_ConservedWaters_PASSED")
  } else {
    PDB.conserved.waters <- conservedWaters.data$filenames.used$conserved.all
    PyMOL.filename <- paste0(filename, "_ConservedWaters_ALL")
  }

  ## GET FILENAMES -------------------------------------------------------------
  ##- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  prefix <- conservedWaters.data$filenames.used$prefix
  prefix.PDB.files <- list.files(prefix)
  if ( conservedWaters.data$MDS == TRUE ){
    PDBid.ref.filename <- prefix.PDB.files[grepl(pattern=PDBid.ref, x=prefix.PDB.files)]
    PDBid.ref.filename <- file.path(prefix, PDBid.ref.filename)
  } else {
    PDBid.ref <- tolower(PDBid.ref)
    PDBid.ref.filename <- prefix.PDB.files[grepl(pattern=PDBid.ref, x=prefix.PDB.files)]
    PDBid.ref.filename <- file.path(prefix, PDBid.ref.filename)
  }

  ## CHECK USER PROVIDED PARAMETERS --------------------------------------------
  ##- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  stop.tf <- FALSE
  ##_ does PDBid.ref.filename exist? -----
  if( file.exists(PDBid.ref.filename) == FALSE ) {
    stop.tf <- TRUE
    mess <- paste("The provided PDBid.ref (PDB reference structure) file (",
                  PDBid.ref, ") does not exist. Please make sure the reference ",
                  "PDBid is present in the ", prefix, " directory. This function ",
                  "is trying to find the file: ", PDBid.ref.filename, sep = "")
    message(mess)
  }
  ##_ stop writing PyMOL script if files do NOT exist -----
  if (stop.tf == TRUE) {
    mess <- paste("The PyMOL script cannot be created until the above user",
                  "provided parameters are corrected", sep = " ")

    stop(mess)
  }

  ##_ does LigResname.ref exist? -----
  if ( is.null(LigResname.ref) ) {
    LigResname.provided.tf <- FALSE
    LigResname.ref.tf <- FALSE
    LIG.present <- FALSE
  } else {
    LigResname.provided.tf <- TRUE
    LIG.present <- TRUE
    LigResname.ref <- toupper(LigResname.ref)
    ##- load reference PDB
    PDB.ref <- bio3d::read.pdb2(PDBid.ref.filename)
    LigResname.ref.tf <- any(unique(PDB.ref$atom$resid) %in% LigResname.ref)
  }
  ##__ no LigResname.ref -----
  if( LigResname.provided.tf != TRUE ) {
    mess <- paste("<<<|||>>>The provided LigResname.ref (ligand residue name within ",
                  "the reference PDB; ", LigResname.ref, ") does not exist. ",
                  "All ligand-like (organic) structures will be identified as ",
                  "ligand.<<<|||>>>\n", sep = "")
    message(mess)
    LigResname.ref <- NULL
    LIG.present <- FALSE
  }


  ## DETERMINE THE "Color By Atom" COLOR SCHEME --------------------------------
  ##- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## http://www.pymolwiki.org/index.php/Advanced_Coloring
  ##__ the available carbon color schemes -----
  colors.available <- c("green",  "cyan",   "magenta",
                        "yellow", "salmon", "grey",
                        "slate",  "orange", "purple", "pink")
  cba.opts <- c("cmd.util.cbag", "cmd.util.cbac", "cmd.util.cbam",
                "cmd.util.cbay", "cmd.util.cbas", "cmd.util.cbaw",
                "cmd.util.cbab", "cmd.util.cbao", "cmd.util.cbap",
                "cmd.util.cbak")
  ##__ unify spelling of grey -----
  if (lig.carbon.color == "gray") {
    lig.carbon.color <- "grey"
  }
  # if (poc.carbon.color == "gray") poc.carbon.color <- "grey"

  ##__ if provided LIGAND color is not an option, set lig.carbon.color to cyan -----
  if (any(lig.carbon.color == colors.available) == FALSE) {
    color.mess <- paste("The provided LIGAND color ", lig.carbon.color,
                        " is not a valid option. In the future please select ",
                        "one of the following colors: ",
                        toString(colors.available), ". The carbon color has ",
                        "been set to 'cyan'; the default.", sep = "")
    message(color.mess)
    lig.carbon.color <- "cyan"
  }
  ##__ set the color by atom scheme -----
  lig.cba.scheme <- cba.opts[which(lig.carbon.color == colors.available)]


  ## BASIC SESSION INFORMATION FOR PyMOL FILE HEADER ---------------------------
  ##- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  system.info <- Sys.info()
  pymol.pre <- c("##-----------------")
  pymol.pre <- c(pymol.pre, paste("## created:", format(Sys.Date(), "%B %d, %Y"),
                                  "at", format(Sys.time(), "%H:%M"), sep = " ") )
  pymol.pre <- c(pymol.pre, paste("## user:", Sys.getenv("LOGNAME"), sep = " ") )
  pymol.pre <- c(pymol.pre, paste("## hostname:",
                                  as.vector(system.info[names(system.info) == "nodename"])
                                  , sep=" ") )
  pymol.pre <- c(pymol.pre, "## this PyMOL script file was created by the CreatePyMOLscript function")
  pymol.pre <- c(pymol.pre, "## of the vanddraabe R (r-project.org) package")
  pymol.pre <- c(pymol.pre, paste("## vanddraabe version:", packageVersion("vanddraabe"), sep =" ") )
  pymol.pre <- c(pymol.pre, paste("## R version:", getRversion(), sep =" ") )
  pymol.pre <- c(pymol.pre, "##-----------------")
  pymol.pre <- c(pymol.pre, "")

  ##__ combine the PyMOL script commands -----
  pymol.script.black <- c(pymol.pre)
  pymol.script.white <- c(pymol.pre)


  ## SETUP THE PyMOL SESSION ---------------------------------------------------
  ##- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pymol.setup <- c("from pymol import cmd")
  pymol.setup <- c(pymol.setup, "from pymol.cgo import *")
  pymol.setup <- c(pymol.setup, "")
  pymol.setup <- c(pymol.setup, "")
  pymol.setup <- c(pymol.setup, "##_ the setup")
  pymol.setup <- c(pymol.setup, "cmd.set('valence', 'on')")
  pymol.setup <- c(pymol.setup, "cmd.set('two_sided_lighting', 'on')")
  pymol.setup <- c(pymol.setup, "cmd.set('surface_quality', '1')")
  pymol.setup <- c(pymol.setup, "##__ set size of the PyMOL graphics window")
  pymol.setup <- c(pymol.setup, "cmd.viewport(1100, 850)" )
  pymol.setup <- c(pymol.setup, "##__ set background")
  ##-- combine the PyMOL script commands -----
  pymol.script.black <- c(pymol.script.black, pymol.setup)
  pymol.script.white <- c(pymol.script.white, pymol.setup)
  ##__ black background -----
  pymol.script.black <- c(pymol.script.black, paste("cmd.bg_color('black')", sep = "") )
  pymol.script.black <- c(pymol.script.black, "")
  pymol.script.black <- c(pymol.script.black, "")
  ##__ white background -----
  pymol.script.white <- c(pymol.script.white, "cmd.set('opaque_background', 'off')")
  pymol.script.white <- c(pymol.script.white, paste("cmd.bg_color('white')", sep = "") )
  pymol.script.white <- c(pymol.script.white, "")
  pymol.script.white <- c(pymol.script.white, "")


  ## CREATE THE COLORS FOR PERCENT CONSERVATION --------------------------------
  ##- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## these color values are based on the Water Conservation plot colors
  ## (cons.color6) from the constants.R file. the hex value was converted to
  ## RGB values using http://color-hex.com. the RGB values were then divided by
  ## 255 to obtain the PyMOL color values.
  pymol.colors <- c("##_ create colors for percent conservation to match the plots")
  ##__ less than 50% [grey; #d9d9d9; (217,217,217)] -----
  pymol.colors <- c(pymol.colors, "cmd.set_color('consLT50_grey', [0.851, 0.851, 0.851] )")
  ##__ greater than or equal to 50% and less than 60% [red; #a50f15; (165,15,21)] -----
  pymol.colors <- c(pymol.colors, "cmd.set_color('cons50_69_red', [0.6471, 0.05882, 0.08235] )")
  ##__ greater than or equal to 60% and less than 70% [medium red; #ef3b2c; (239,59,44)] -----
  pymol.colors <- c(pymol.colors, "cmd.set_color('cons70_79_mRed', [0.9373, 0.2314, 0.1725] )")
  ##__ greater than or equal to 80% and less than 90% [light blue; #9ecae1; (158,202,225)] -----
  pymol.colors <- c(pymol.colors, "cmd.set_color('cons80_89_lBlue', [0.6196, 0.7922, 0.8824] )")
  ##__ greater than or equal to 90% and less than 100% [medium blue; #4292c6; (66,146,198)] -----
  pymol.colors <- c(pymol.colors, "cmd.set_color('cons90_99_mBlue', [0.2588, 0.5725, 0.7765] )")
  ##__ equal to 100% [dark blue; #08519c; (8,81,156)] -----
  pymol.colors <- c(pymol.colors, "cmd.set_color('cons100_dBlue', [0.03137, 0.31765, 0.61176] )")
  pymol.colors <- c(pymol.colors, "")
  pymol.colors <- c(pymol.colors, "")

  ##__ combine the PyMOL script commands -----
  pymol.script.black <- c(pymol.script.black, pymol.colors)
  pymol.script.white <- c(pymol.script.white, pymol.colors)


  ## LOAD THE REFERENCE & CONSERVED WATER PDB FILE -----------------------------
  ##- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pymol.load <- c("##_ load the reference structure and the waters")
  pymol.load <- c(pymol.load, paste("cmd.load('", PDBid.ref.filename, "', '",
                                    PDBid.ref, "')", sep = ""))
  ##__ remove waters from reference structure -----
  pymol.load <- c(pymol.load, "cmd.remove('resn HOH+DOD+WAT')")
  pymol.load <- c(pymol.load, paste("cmd.load('", PDB.conserved.waters,
                                    "', 'ConservedWaters')", sep = ""))
  pymol.load <- c(pymol.load, "cmd.h_add()")
  pymol.load <- c(pymol.load, "")
  pymol.load <- c(pymol.load, "")

  ##__ combine the PyMOL script commands -----
  pymol.script.black <- c(pymol.script.black, pymol.load)
  pymol.script.white <- c(pymol.script.white, pymol.load)


  ## "SELECT" THE PROTEIN, LIGAND, & COUNTER IONS ------------------------------
  ##- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##__ the protein -----
  pymol.selections <- c("##_ create objects ")
  pymol.poly <- paste("cmd.select('ref_prot', 'polymer and ", PDBid.ref, "')", sep = "")

  ##__ the ligand(s) -----
  if ( (length(LigResname.ref) > 1) & ( LIG.present == TRUE ) ) {
    LigResname.ref <- paste0(LigResname.ref, collapse="+")
  }
  ##- if the ligand is present use it, if not all organic-like are the ligand
  if ( LigResname.ref.tf == TRUE ) {
    pymol.orgs <- paste("cmd.select('ligs', 'resn ", LigResname.ref, " and ", PDBid.ref, "')", sep = "")
  } else {
    pymol.orgs <- paste("cmd.select('ligs', 'organic and ", PDBid.ref, "')", sep = "")
  }

  ##__ the counter ions -----
  pymol.inorg <- paste("cmd.select('ions', 'inorganic and ", PDBid.ref, "')", sep = "")

  ##__ the pocket conserved waters -----
  ##__ the pocket residues -----
  if ( LIG.present == TRUE ) {
    pymol.pocWaters <- c("cmd.select('pocWaters', 'byres ligs around 6 and ConservedWaters')")
    pymol.pocket <- c("cmd.select('pocResidues', 'byres pocWaters around 5 and ref_prot')")
    pymol.pocket <- c(pymol.pocket, "cmd.select('pocCartoon', 'byres ligs around 15 and ref_prot')")
    pymol.pocket <- c(pymol.pocket, "cmd.select('pocSurface', 'byres ligs around 5 and ref_prot')")
  } else {
    pymol.pocWaters <- c("cmd.select('pocWaters', 'byres ref_prot and (q>0.399999) around 4 and ConservedWaters')")
    pymol.pocket <- c("cmd.select('pocResidues', 'byres pocWaters around 5 and ref_prot')")
    pymol.pocket <- c(pymol.pocket, "cmd.select('pocCartoon', 'byres pocWaters around 15 and ref_prot')")
    pymol.pocket <- c(pymol.pocket, "cmd.select('pocSurface', 'byres pocWaters around 5 and ref_prot')")
  }
  ##__ combine the different selections -----
  pymol.selections <- c(pymol.selections, pymol.poly, pymol.orgs, pymol.inorg,
                        pymol.pocWaters, pymol.pocket)
  ##__ deselect everything -----
  pymol.selections <- c(pymol.selections, "cmd.deselect()")
  ##__ hide everything -----
  pymol.selections <- c(pymol.selections, "##__ hide everything")
  pymol.selections <- c(pymol.selections, "cmd.hide()")

  ##__ combine the PyMOL script commands -----
  pymol.script.black <- c(pymol.script.black, pymol.selections)
  pymol.script.white <- c(pymol.script.white, pymol.selections)


  ## RENDER THE REFERENCE PROTEIN ----------------------------------------------
  ##- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pymol.refProt <- c("##_ render the reference protein")
  ##__ render the protein -----
  pymol.refProt <- c(pymol.refProt, "##__ the protein")
  pymol.refProt <- c(pymol.refProt, "cmd.show('cartoon', 'pocCartoon')")
  pymol.refProt <- c(pymol.refProt, "cmd.show('sticks', 'pocResidues')")
  pymol.refProt <- c(pymol.refProt, "cmd.hide('sticks', 'pocResidues and hydro')")
  pymol.refProt <- c(pymol.refProt, "##-- the pocket residues")

  ##-- combine the PyMOL script commands -----
  pymol.script.black <- c(pymol.script.black, pymol.refProt)
  pymol.script.white <- c(pymol.script.white, pymol.refProt)

  ##__ render the pocket residues -----
  ##    set color of pocket sticks based on background
  ##-- black background -----
  pymol.script.black <- c(pymol.script.black, "cmd.color('grey80', 'pocResidues')")
  ##-- white background -----
  pymol.script.white <- c(pymol.script.white, "cmd.color('grey40', 'pocResidues')")
  ##__ set the color of the non-carbon atoms -----
  # pymol.refProt <- c("cmd.util.cnc('pocResidues')")
  ##__ set the transparency value for the pocket residues -----
  pymol.refProt <- c(pymol.refProt, "cmd.set_bond('stick_transparency', '0.5', 'pocResidues')")

  ##__ render the secondary structure representation -----
  pymol.refProt <- c(pymol.refProt, "##__ render the secondary structure representation")
  pymol.refProt <- c(pymol.refProt, paste("cmd.set('cartoon_color', 'red', '", PDBid.ref, " and ss h')", sep = "") )
  pymol.refProt <- c(pymol.refProt, paste("cmd.set('cartoon_color', 'yellow', '", PDBid.ref, " and ss s')", sep = "") )
  pymol.refProt <- c(pymol.refProt, paste("cmd.set('cartoon_color', 'green', '", PDBid.ref, " and ss l+')", sep = "") )
  pymol.refProt <- c(pymol.refProt, paste("cmd.set('cartoon_transparency', '0.75', '", PDBid.ref, "')", sep = "") )

  ##__ render the ligand(s) -----
  if ( LIG.present == TRUE ) {
    pymol.refProt <- c(pymol.refProt, "##__ the ligand(s) with only polar-hydrogens")
    pymol.refProt <- c(pymol.refProt, "cmd.show('sticks', 'ligs')")
    pymol.refProt <- c(pymol.refProt, "cmd.hide('sticks', 'hydro')")
    pymol.refProt <- c(pymol.refProt, "cmd.show('sticks', 'ligs and ele h and neighbor (ele n+o)')")
    pymol.refProt <- c(pymol.refProt,
                       paste(lig.cba.scheme, "('ligs')", sep = "") )
  }

  ##__ render the ions -----
  pymol.refProt <- c(pymol.refProt, "##__ the ions")
  pymol.refProt <- c(pymol.refProt, "cmd.set('sphere_scale', '0.6', 'ions')")
  pymol.refProt <- c(pymol.refProt, "cmd.show('spheres', 'ions')")
  pymol.refProt <- c(pymol.refProt, "cmd.util.cnc('ions')")
  ##__ create the pocket and surface -----
  pymol.refProt <- c(pymol.refProt, "##__ the pocket and surface")
  pymol.refProt <- c(pymol.refProt, "cmd.show('surface', 'pocSurface')")
  pymol.refProt <- c(pymol.refProt, "cmd.set('surface_color', 'grey90')")
  pymol.refProt <- c(pymol.refProt, paste("cmd.set('transparency', '0.4', '", PDBid.ref, "')", sep = "") )
  pymol.refProt <- c(pymol.refProt, "")
  pymol.refProt <- c(pymol.refProt, "")

  ##-- combine the PyMOL script commands -----
  pymol.script.black <- c(pymol.script.black, pymol.refProt)
  pymol.script.white <- c(pymol.script.white, pymol.refProt)


  ## RENDER THE CONSERVED WATERS -----------------------------------------------
  ##- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pymol.consWat <- c("##_ render the conserved waters")
  ##__ render ths size of the conserved waters based on mean b-value -----
  pymol.consWat <- c(pymol.consWat, "##__ display as spheres with radius based on b-value")
  pymol.consWat <- c(pymol.consWat, "cmd.alter('ConservedWaters', 'vdw=((b/10) / (8 * 3.14**2))**0.5 * 3')")
  pymol.consWat <- c(pymol.consWat, "##__ only display conserved waters with 40% or greater conservation")
  if ( LIG.present == TRUE ) {
    pymol.consWat <- c(pymol.consWat, "cmd.show('spheres', 'pocWaters and (q>0.399999) and not hydro')")
  } else {
    pymol.consWat <- c(pymol.consWat, "cmd.show('spheres', 'ConservedWaters and (q>0.399999) and not hydro')")
  }
  # pymol.consWat <- c(pymol.consWat, "cmd.set('sphere_transparency', '0.20', 'pocWaters')")
  ##__ color the conserved waters based on percent conservation -----
  pymol.consWat <- c(pymol.consWat, "##__ color conserved waters base on percent conservation")
  pymol.consWat <- c(pymol.consWat, "cmd.color('consLT50_grey', 'ConservedWaters and (q<0.50)')")
  pymol.consWat <- c(pymol.consWat, "cmd.color('cons50_69_red', 'ConservedWaters and (q>0.49999 and q<0.70)')")
  pymol.consWat <- c(pymol.consWat, "cmd.color('cons70_79_mRed', 'ConservedWaters and (q>0.69999 and q<0.80)')")
  pymol.consWat <- c(pymol.consWat, "cmd.color('cons80_89_lBlue', 'ConservedWaters and (q>0.79999 and q<0.90)')")
  pymol.consWat <- c(pymol.consWat, "cmd.color('cons90_99_mBlue', 'ConservedWaters and (q>0.89999 and q<1.0)')")
  pymol.consWat <- c(pymol.consWat, "cmd.color('cons100_dBlue', 'ConservedWaters and (q=1.00)')")
  pymol.consWat <- c(pymol.consWat, "cmd.label('pocWaters and (q>0.399999) and not hydro', 'resi')")
  pymol.consWat <- c(pymol.consWat, "")
  pymol.consWat <- c(pymol.consWat, "")

  ##-- combine the PyMOL script commands -----
  pymol.script.black <- c(pymol.script.black, pymol.consWat)
  pymol.script.white <- c(pymol.script.white, pymol.consWat)


  ## RENDER THE HYDROGEN BONDS -------------------------------------------------
  ##- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pymol.hbonds <- c("##_render the hydrogen bonds")
  ##__ hydrogen bond display parameters -----
  pymol.hbonds <- c(pymol.hbonds, "##__ hydrogen bond display parameters")
  pymol.hbonds <- c(pymol.hbonds, "cmd.set('dash_color', 'grey80')")
  pymol.hbonds <- c(pymol.hbonds, "cmd.set('dash_gap', '0.3')")
  pymol.hbonds <- c(pymol.hbonds, "cmd.set('dash_length', '0.15')")
  pymol.hbonds <- c(pymol.hbonds, "cmd.set('dash_radius', '0.15')")
  pymol.hbonds <- c(pymol.hbonds, "cmd.set('dash_round_ends', 'on')")
  pymol.hbonds <- c(pymol.hbonds, "cmd.set('label_size', '-0.7')")
  pymol.hbonds <- c(pymol.hbonds, "cmd.set('label_color', 'white')")
  pymol.hbonds <- c(pymol.hbonds, "cmd.set('label_font_id', '7')")
  ##__ identify hydrogen bond acceptors and donors -----
  pymol.hbonds <- c(pymol.hbonds, "##__ identify hydrogen bond acceptors and donors")
  pymol.hbonds <- c(pymol.hbonds, "cmd.select('hbd', '(elem n,o and (neighbor hydro))')")
  pymol.hbonds <- c(pymol.hbonds, "cmd.select('hba', '(elem o or (elem n and not (neighbor hydro)))')")
  ##__ hydrogen bonds between conserved waters and ligand -----
  if ( LIG.present == TRUE ){
    pymol.hbonds <- c(pymol.hbonds, "##__ hydrogen bonds between conserved waters and ligand")
    pymol.hbonds <- c(pymol.hbonds, paste("cmd.distance('wl_hba', '(ligs and hba)', '(pocWaters and hbd and (q>0.399999))', '", hbond, "')", sep = "") )
    pymol.hbonds <- c(pymol.hbonds, paste("cmd.distance('wl_hbd', '(ligs and hbd)', '(pocWaters and hba and (q>0.399999))', '", hbond, "')", sep = "") )
    pymol.hbonds <- c(pymol.hbonds, "cmd.set('dash_color', 'brightorange', 'wl_hba')")
    pymol.hbonds <- c(pymol.hbonds, "cmd.set('dash_color', 'brightorange', 'wl_hbd')")
  }
  ##__ hydrogen bonds between conserved waters and protein -----
  pymol.hbonds <- c(pymol.hbonds, "##__ hydrogen bonds between conserved waters and protein")
  pymol.hbonds <- c(pymol.hbonds, paste("cmd.distance('wp_hba', '(ref_prot and hba)', '(pocWaters and hbd and (q>0.399999))', '", hbond, "')", sep = "") )
  pymol.hbonds <- c(pymol.hbonds, paste("cmd.distance('wp_hbd', '(ref_prot and hbd)', '(pocWaters and hba and (q>0.399999))', '", hbond, "')", sep = "") )
  pymol.hbonds <- c(pymol.hbonds, "cmd.set('dash_color', 'forest', 'wp_hba')")
  pymol.hbonds <- c(pymol.hbonds, "cmd.set('dash_color', 'forest', 'wp_hbd')")
  ##__ hydrogen bonds between conserved waters -----
  pymol.hbonds <- c(pymol.hbonds, "##__ hydrogen bonds between conserved waters")
  pymol.hbonds <- c(pymol.hbonds, paste("cmd.distance('ww_hba', '(pocWaters and hba and (q>0.399999))', '(pocWaters and hbd and (q>0.399999))', '", hbond, "')", sep = "") )
  pymol.hbonds <- c(pymol.hbonds, paste("cmd.distance('ww_hbd', '(pocWaters and hbd and (q>0.399999))', '(pocWaters and hba and (q>0.399999))', '", hbond, "')", sep = "") )
  pymol.hbonds <- c(pymol.hbonds, "cmd.set('dash_color', 'marine', 'ww_hba')")
  pymol.hbonds <- c(pymol.hbonds, "cmd.set('dash_color', 'marine', 'ww_hbd')")
  pymol.hbonds <- c(pymol.hbonds, "")
  pymol.hbonds <- c(pymol.hbonds, "")

  ##-- combine the PyMOL script commands -----
  pymol.script.black <- c(pymol.script.black, pymol.hbonds)
  pymol.script.white <- c(pymol.script.white, pymol.hbonds)


  ## CLEANUP REPRESENTATION ----------------------------------------------------
  ##- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pymol.cleanup <- c("##_ cleanup the representation")

  pymol.cleanup <- c(pymol.cleanup, "cmd.deselect()")

  pymol.cleanup <- c(pymol.cleanup, "cmd.hide('labels', 'w*hb*')")
  # pymol.cleanup <- c(pymol.cleanup, "cmd.delete('pocCartoon')")
  # pymol.cleanup <- c(pymol.cleanup, "cmd.delete('pocSurface')")
  # pymol.cleanup <- c(pymol.cleanup, "cmd.delete('hba')")
  # pymol.cleanup <- c(pymol.cleanup, "cmd.delete('hbd')")
  pymol.cleanup <- c(pymol.cleanup, "")
  pymol.cleanup <- c(pymol.cleanup, "")

  ##-- combine the PyMOL script commands -----
  pymol.script.black <- c(pymol.script.black, pymol.cleanup)
  pymol.script.white <- c(pymol.script.white, pymol.cleanup)


  ## ORIENT THE SYSTEM --------------------------------------------------------_
  ##- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pymol.view <- c("##_ orient the system")
  if ( LIG.present == TRUE ) {
    pymol.view <- c(pymol.view, "cmd.orient('ligs')")
    pymol.view <- c(pymol.view, "cmd.zoom('ligs', '10')")
    pymol.view <- c(pymol.view, "")
    pymol.view <- c(pymol.view, "")
  } else {
    pymol.view <- c(pymol.view, "cmd.orient('pocWaters')")
    pymol.view <- c(pymol.view, "cmd.zoom('pocWaters', '10')")
    pymol.view <- c(pymol.view, "")
    pymol.view <- c(pymol.view, "")
  }

  ##-- combine the PyMOL script commands -----
  pymol.script.black <- c(pymol.script.black, pymol.view)
  pymol.script.white <- c(pymol.script.white, pymol.view)


  ## WRITE PyMOL SCRIPT --------------------------------------------------------
  ##- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  date.time <- conservedWaters.data$filenames.used$when

  ##_ create the filenames -----
  filename.black.pml <- paste(PyMOL.filename,
                              "_PyMOL_black_background_",
                              date.time,
                              ".pml", sep = "")

  filename.white.pml <- paste(PyMOL.filename,
                              "_PyMOL_white_background_",
                              date.time,
                              ".pml", sep = "")

  ##_ write the scripts to the .py files -----
  write(pymol.script.black, file=filename.black.pml)
  write(pymol.script.white, file=filename.white.pml)


  ## WRITE SUMMARY BLOCK -------------------------------------------------------
  ##- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##_ input files provided -----

  ##_ ligand residue name used -----

  ##_ PyMOL session files written -----
  files.location <- paste("The PyMOL script files (", filename.black.pml,
                          " and ", filename.white.pml,
                          ") to easily view the conserved",
                          "waters are in the ", getwd(),
                          " directory (folder).", sep = "")
  message(files.location)

}
