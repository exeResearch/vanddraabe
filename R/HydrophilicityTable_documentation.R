## HydrophilicityTable_documentation.R :: hydrophilicity table documentation
##
## apr-03-2016 (exe) created
## jan-20-2017 (exe) constants for openxlsx cell styles
## jan-20-2017 (exe) data.frame for hydrophilicity values
## jan-20-2017 (exe) constant for plot color scheme values
## jan-20-2017 (exe) constants for atom and residue names
## jan-20-2017 (exe) moved to data-raw to construct R/sysdata.rda
## jan-20-2017 (exe) documentation is in R/Constants_docs.R
## mar-27-2017 (exe) updated documentation to reflect Top8000 based values
## jul-25-2017 (exe) updated documentation
## jul-25-2017 (exe) changed name (HydrophilicityTable.R -> HydrophilicityTable_documentation.R) for clarity
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


## Residue Atom Type Hydrophilicity Values docs --------------------------------
#' @title Residue Atom Type Hydrophilicity Values
#' @description Atomic hydrophilicity values for the 20 naturally occurring
#'   amino acids and water.
#' @details The Hydrophilicity Table is based on the work of Esposito (see
#'   reference below) in [vanddraabe] package. The hydrophilicity values are
#'   based on information from a 1995 analysis of published PDB structures and
#'   indicate how likely the individual atoms of the amino acid residues are to
#'   have a water molecule within 4.0 angstroms.
#'
#'   The data contained within the Hydrophilicity Table is based on ~7900
#'   experimentally determined crystallographic protein structures with
#'   resolution values less than or equal to x.x Angstroms, a R-factor less than
#'   or equal to 0.26, and 20 or more bound waters each. The protein structures
#'   are from the Top8000 ("a database of about 8000 high-resolution,
#'   quality-filtered protein chains"; reference below) high quality protein
#'   dataset from the [Kinemage laboratory](http://kinemage.biochem.duke.edu) at
#'   Duke University. The included structures had a range of B-values and
#'   occupancy values.
#'
#'   These values are based on the methods and protocols of Kuhn _et al_.
#'
#'   The Hydrophilicity Table contains:
#'   - **residueAtomName**: Contracted residue type and atom name to aid
#'     looking up hydrophilicity values.
#'   - **residue**: Three-letter residue name.
#'   - **atomName**: Atom name indicating the atom type and its position in
#'     the amino acid residue.
#'   - **surfaceOccurrences**: Number of times each
#'     atom has a defined solvent exposed surface area.
#'   - **hydratOccurrences**: Proportion of the solvent exposed residue-specific
#'     atom type with a water molecule closely bound (within 4.0 Angstroms).
#'
#' @name HydrophilicityTable
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @references
#'   Leslie A Kuhn, Craig A Swanson, Michael E Pique, John A Tainer,
#'   and Elizabeth D Getzof. Atomic and Residue Hydrophilicity in the Context of
#'   Folded Protein Structures. _PROTEINS: Structure, Function, and
#'   Genetics_, 1995, **23** (_4_), pp 536-547.
#'   [DOI: 10.1002/prot.340230408](http://doi.org/10.1002/prot.340230408)
#'   [PMID: 8749849](http://www.ncbi.nlm.nih.gov/pubmed/8749849)
#'
#'   Bradley J Hintze, Steven M Lewis, Jane S Richardson, and David C
#'   Richardson. Molprobity's ultimate rotamer-library distributions for model
#'   validation. _Proteins: Structure, Function, and Bioinformatics_, 2016,
#'   **84** (_9_), pp 1177-1189.
#'   [DOI: 10.1002/prot.25039](https://dx.doi.org/10.1002/prot.25039)
#'   [Top8000 webpage](http://kinemage.biochem.duke.edu/databases/top8000.php)
#'
NULL
