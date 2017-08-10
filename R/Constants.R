## Constants.R :: constant values for ConservedWaters
##
## apr-03-2016 (exe) created
## jan-20-2017 (exe) constants for openxlsx cell styles
## jan-20-2017 (exe) data.frame for hydrophilicity values
## jan-20-2017 (exe) constant for plot color scheme values
## jan-20-2017 (exe) constants for atom and residue names
## jan-20-2017 (exe) moved to data-raw to construct R/sysdata.rda
## jan-20-2017 (exe) documentation is in R/Constants_docs.R
## jan-23-2017 (exe) exported names.*
## mar-27-2017 (exe) update documentation
## jul-07-2017 (exe) update documentation
## jul-25-2017 (exe) updated documentation
## jul-25-2017 (exe) constants are now accessible to the user, no longer in sysdata.rda
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


## backbone atom names docs ----------------------------------------------------
#' @title Backbone Atom Names
#' @description Backbone atom names based on PDB atom naming conventions.
#' @details Protein backbone atom names based on the PDB atom naming
#'   conventions.
#'
#'   - **N**: Nitrogen backbone atom; amide, "leading" functional group
#'   - **CA**: alpha-Carbon backbone atom; bonds/connects the side chain to
#'     the backbone
#'   - **C**: Carbon backbone atom; carboxyl, "tail" functional group
#'   - **O**: Oxygen backbone atom double bonded to the carbon backbone (C)
#'     atom; part of the carboxyl, "tail" functional groups
#'
#' @name names.backbone.atoms
#'
#' @examples
#'   names.backbone.atoms
#'   # [1] "N"  "CA" "C"  "O"
#'
#' @export
#'
#' @family constants
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
names.backbone.atoms <- c("N", "CA", "C", "O")


## sidechain atom names docs ---------------------------------------------------
#' @title Sidechain Atom Names
#' @description Sidechain atom names based on PDB atom naming conventions.
#' @details The 32 unique sidechain atom names. The first character is the
#'   element and the second character is the Greek letter (B=beta, D=delta,
#'   E=epsilon, G=gamma, Z=zeta) defining the specific position within the
#'   sidechain. The exception to the use of Greek letters is `OH`
#'   indicating a hydroxyl group at the _para_ position of the six-member
#'   ring of tyrosine. Some sidechain atom names have a number in the third
#'   character position when there are mirrored/symetrical atoms; _e.g._,
#'   `CG1` and `CG2` of valine.
#'
#' @name names.sidechain.atoms
#'
#' @examples
#'   names.sidechain.atoms
#'   # [1] "CB"  "CG"  "CD"  "NE"  "CZ"  "NH1" "NH2" "OD1" "ND2" "OD2" "SG"
#'   #     "OE1" "NE2" "OE2" "CD2" "ND1" "CE1" "CG1" "CG2" "CD1" "CE"  "NZ"
#'   #     "SD"  "CE2" "OG"  "OG1" "NE1" "CE3" "CZ2" "CZ3" "CH2" "OH"
#'
#' @export
#'
#' @family constants
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
names.sidechain.atoms <- c("CB", "CG", "CD", "NE", "CZ", "NH1", "NH2", "OD1",
                           "ND2", "OD2", "SG", "OE1", "NE2", "OE2", "CD2",
                           "ND1", "CE1", "CG1", "CG2", "CD1", "CE", "NZ", "SD",
                           "CE2", "OG", "OG1", "NE1", "CE3", "CZ2", "CZ3",
                           "CH2", "OH")


## polar atom names docs -------------------------------------------------------
#' @title Polar Atom Names
#' @description Polar atom names based on PDB atom naming conventions.
#' @details Polar atoms are those possessing a lone pair(s) of elections able to
#'   participate in hydrogen bonds with hydrogen atoms within 3.5 Angstroms and
#'   XX degrees of the lone pair containing atom. Traditionally, nitrogen,
#'   oxygen, and sulfur atoms possess lone pair(s) of electrons and participate
#'   in hydrogen bonds in biological sytems. Water molecules are able to
#'   hydrogen bond with and participate in hydrogen bonds.
#'
#' @name names.polar.atoms
#'
#' @examples
#'   names.polar.atoms
#'   # [1] "N"   "NE"  "NH1" "NH2" "ND2" "NE2" "ND1" "NZ"  "NE1" "O"   "OD1"
#'   #     "OD2" "OE1" "OE2" "OG"  "OG1" "OH"  "S"   "SD"  "SG"
#'
#' @export
#'
#' @family constants
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
names.polar.atoms <- c("N", "NE", "NH1", "NH2", "ND2", "NE2", "ND1", "NZ",
                       "NE1", "O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1",
                       "OH", "S", "SD", "SG")


## neutral oxygen residue-atomtype names docs ----------------------------------
#' @title Neutral Oxygen Residue-AtomType Names
#' @description Neutral oxygen residue-atomtype names based on PDB atom naming
#'   conventions.
#' @details These residue-atomtype names indicate oxygen atoms with a neutral
#'   charge.
#'
#' @examples
#'  names.resATs.oxy.neut
#'  # [1] "ALA O"   "ARG O"   "ASN O"   "ASN OD1" "ASP O"   "CYS O"   "GLN O"
#'  # "GLN OE1" "GLU O"   "GLY O"   "HIS O"   "ILE O"   "LEU O"   "LYS O"   "MET O"
#'  # [16] "PHE O"   "PRO O"   "SER O"   "SER OG"  "THR O"   "THR OG1" "TRP O"
#'  # "TYR O"   "TYR OH"  "VAL O"
#'
#' @export
#'
#' @family constants
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
names.resATs.oxy.neut <- c("ALA O",   "ARG O", "ASN O",   "ASN OD1", "ASP O",
                           "CYS O",   "GLN O", "GLN OE1", "GLU O",   "GLY O",
                           "HIS O",   "ILE O", "LEU O",   "LYS O",   "MET O",
                           "PHE O",   "PRO O", "SER O",   "SER OG",  "THR O",
                           "THR OG1", "TRP O", "TYR O",   "TYR OH",  "VAL O")


## negative oxygen residue-atomtype names docs ---------------------------------
#' @title Negative Oxygen Residue-AtomType Names
#' @description Negatvie oxygen residue-atomtype names based on PDB atom naming
#'   conventions.
#' @details These residue-atomtype names indicate oxygen atoms with a negative
#'   charge.
#'
#' @examples
#'  names.resATs.oxy.neg
#'  # [1] "ASP OD1" "ASP OD2" "GLU OE1" "GLU OE2"
#'
#' @export
#'
#' @family constants
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
names.resATs.oxy.neg <- c("ASP OD1", "ASP OD2", "GLU OE1", "GLU OE2")


## neutral nitrogen residue-atomtype names docs --------------------------------
#' @title Neutral Nitrogen Residue-AtomType Names
#' @description Neutral nitrogen residue-atomtype names based on PDB atom naming
#'   conventions.
#' @details These residue-atomtype names indicate nitrogen atoms with a neutral
#'   charge.
#'
#' @examples
#'  names.resATs.nitro.neut
#'  # [1] "ALA N"   "ARG N"   "ARG NE"  "ASN N"   "ASN ND2" "ASP N"   "CYS N"
#'  # "GLN N"   "GLN NE2" "GLU N"   "GLY N"   "HIS N"   "ILE N"   "LEU N"   "LYS N"
#'  # [16] "MET N"   "PHE N"   "PRO N"   "SER N"   "THR N"   "TRP N"   "TRP NE1"
#'  # "TYR N"   "VAL N"
#'
#' @export
#'
#' @family constants
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
names.resATs.nitro.neut <- c("ALA N", "ARG N",   "ARG NE", "ASN N", "ASN ND2",
                             "ASP N", "CYS N",   "GLN N",  "GLN NE2", "GLU N",
                             "GLY N", "HIS N",   "ILE N",  "LEU N", "LYS N",
                             "MET N", "PHE N",   "PRO N",  "SER N", "THR N",
                             "TRP N", "TRP NE1", "TYR N",  "VAL N"  )


## positive nitrogen residue-atomtype names docs -------------------------------
#' @title Positive Nitrogen Residue-AtomType Names
#' @description Positive nitrogen residue-atomtype names based on PDB atom
#'   naming conventions.
#' @details These residue-atomtype names indicate nitrogen atoms with a positive
#'   charge.
#'
#' @examples
#'  names.resATs.nitro.pos
#'  # [1] "ARG NH1" "ARG NH2" "HIS ND1" "HIS NE2" "LYS NZ"
#'
#' @export
#'
#' @family constants
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
names.resATs.nitro.pos <- c("ARG NH1", "ARG NH2",
                            "HIS ND1", "HIS NE2",
                            "LYS NZ")


## carbon and sulfur residue-atomtype names docs -------------------------------
#' @title Carbon and Sulfur Residue-AtomType Names
#' @description Carbon and sulfur residue-atomtype names based on PDB atom
#'   naming conventions.
#' @details These residue-atomtype names indicate carbon and sulfur atoms with a
#'   neutral charge.
#'
#' @examples
#'  names.resATs.carb.sulf
#'  # [1] "ALA CA"  "ALA C"   "ALA CB"  "ARG CA"  "ARG C"   "ARG CB"  "ARG CG"
#'  # "ARG CD"  "ARG CZ"  "ASN CA"  "ASN C"   "ASN CB"  "ASN CG"  "ASP CA"  "ASP C"
#'  # [16] "ASP CB"  "ASP CG"  "CYS CA"  "CYS C"   "CYS CB"  "CYS SG"  "GLN CA"
#'  # "GLN C"   "GLN CB"  "GLN CG"  "GLN CD"  "GLU CA"  "GLU C"   "GLU CB"  "GLU CG"
#'  # [31] "GLU CD"  "GLY CA"  "GLY C"   "HIS CA"  "HIS C"   "HIS CB"  "HIS CG"
#'  # "HIS CD2" "HIS CE1" "ILE CA"  "ILE C"   "ILE CB"  "ILE CG1" "ILE CG2" "ILE CD1"
#'  # [46] "LEU CA"  "LEU C"   "LEU CB"  "LEU CG"  "LEU CD1" "LEU CD2" "LYS CA"
#'  # "LYS C"   "LYS CB"  "LYS CG"  "LYS CD"  "LYS CE"  "MET CA"  "MET C"   "MET CB"
#'  # [61] "MET CG"  "MET SD"  "MET CE"  "PHE CA"  "PHE C"   "PHE CB"  "PHE CG"
#'  # "PHE CD1" "PHE CD2" "PHE CE1" "PHE CE2" "PHE CZ"  "PRO CA"  "PRO C"   "PRO CB"
#'  # [76] "PRO CG"  "PRO CD"  "SER CA"  "SER C"   "SER CB"  "THR CA"  "THR C"
#'  # "THR CB"  "THR CG2" "TRP CA"  "TRP C"   "TRP CB"  "TRP CG"  "TRP CD1" "TRP CD2"
#'  # [91] "TRP CE2" "TRP CE3" "TRP CZ2" "TRP CZ3" "TRP CH2" "TYR CA"  "TYR C"
#'  # "TYR CB"  "TYR CG"  "TYR CD1" "TYR CD2" "TYR CE1" "TYR CE2" "TYR CZ"  "VAL CA"
#'  # [106] "VAL C"   "VAL CB"  "VAL CG1" "VAL CG2"
#'
#' @export
#'
#' @family constants
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
names.resATs.carb.sulf <- c("ALA CA", "ALA C", "ALA CB",
                            "ARG CA", "ARG C", "ARG CB", "ARG CG", "ARG CD", "ARG CZ",
                            "ASN CA", "ASN C", "ASN CB", "ASN CG",
                            "ASP CA", "ASP C", "ASP CB", "ASP CG",
                            "CYS CA", "CYS C", "CYS CB", "CYS SG",
                            "GLN CA", "GLN C", "GLN CB", "GLN CG", "GLN CD",
                            "GLU CA", "GLU C", "GLU CB", "GLU CG", "GLU CD",
                            "GLY CA", "GLY C",
                            "HIS CA", "HIS C", "HIS CB", "HIS CG", "HIS CD2", "HIS CE1",
                            "ILE CA", "ILE C", "ILE CB", "ILE CG1", "ILE CG2", "ILE CD1",
                            "LEU CA", "LEU C", "LEU CB", "LEU CG", "LEU CD1", "LEU CD2",
                            "LYS CA", "LYS C", "LYS CB", "LYS CG", "LYS CD", "LYS CE",
                            "MET CA", "MET C", "MET CB", "MET CG", "MET SD", "MET CE",
                            "PHE CA", "PHE C", "PHE CB", "PHE CG", "PHE CD1", "PHE CD2", "PHE CE1", "PHE CE2", "PHE CZ",
                            "PRO CA", "PRO C", "PRO CB", "PRO CG", "PRO CD",
                            "SER CA", "SER C", "SER CB",
                            "THR CA", "THR C", "THR CB", "THR CG2",
                            "TRP CA", "TRP C", "TRP CB", "TRP CG", "TRP CD1", "TRP CD2", "TRP CE2", "TRP CE3", "TRP CZ2", "TRP CZ3", "TRP CH2",
                            "TYR CA", "TYR C", "TYR CB", "TYR CG", "TYR CD1", "TYR CD2", "TYR CE1", "TYR CE2", "TYR CZ",
                            "VAL CA", "VAL C", "VAL CB", "VAL CG1", "VAL CG2")


## residue names docs ----------------------------------------------------------
#' @title Residue Names
#' @description Residue names based on PDB atom naming conventions.
#' @details The three (3) letter abbreviation for the twenty (20) naturally
#'   occurring amino acid residues.
#'
#' @name names.residues
#'
#' @examples
#'   names.residues
#'   # [1] "ALA" "ARG" "ASN" "ASP" "CYS" "GLN" "GLU" "GLY" "HIS" "ILE"
#'   #     "LEU" "LYS" "MET" "PHE" "PRO" "SER" "THR" "TRP" "TYR" "VAL"
#'
#' @export
#'
#' @family constants
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
names.residues <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
                    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                    "THR", "TRP", "TYR", "VAL")


## residue and atomtype names docs ---------------------------------------------
#' @title Residue and AtomType Names
#' @description Residue and AtomType names based on PDB atom naming conventions.
#' @details The 167 residue-atomtype names based on the 20 naturally occurring
#'   amino acids.
#'
#' @name names.res.AtomTypes
#'
#' @examples
#'   names.res.AtomTypes
#'   # [1] "ALA C"   "ALA CA"  "ALA CB"  "ALA N"   "ALA O"   "ARG C"   "ARG CA"
#'   #     "ARG CB"  "ARG CD"  "ARG CG"  "ARG CZ"  "ARG N"   "ARG NE"  "ARG NH1"
#'   #     "ARG NH2" "ARG O"   "ASN C"   "ASN CA"  "ASN CB"  "ASN CG"  "ASN N"
#'   #     "ASN ND2" "ASN O"   "ASN OD1" "ASP C"   "ASP CA"  "ASP CB"  "ASP CG"
#'   #     "ASP N"   "ASP O"   "ASP OD1" "ASP OD2" "CYS C"   "CYS CA"  "CYS CB"
#'   #     "CYS N"   "CYS O"   "CYS SG"  "GLN C"   "GLN CA"  "GLN CB"  "GLN CD"
#'   #     "GLN CG"  "GLN N"   "GLN NE2" "GLN O"   "GLN OE1" "GLU C"   "GLU CA"
#'   #     "GLU CB"  "GLU CD"  "GLU CG"  "GLU N"   "GLU O"   "GLU OE1" "GLU OE2"
#'   #     "GLY C"   "GLY CA"  "GLY N"   "GLY O"   "HIS C"   "HIS CA"  "HIS CB"
#'   #     "HIS CD2" "HIS CE1" "HIS CG"  "HIS N"   "HIS ND1" "HIS NE2" "HIS O"
#'   #     "ILE C"   "ILE CA"  "ILE CB"  "ILE CD1" "ILE CG1" "ILE CG2" "ILE N"
#'   #     "ILE O"   "LEU C"   "LEU CA"  "LEU CB"  "LEU CD1" "LEU CD2" "LEU CG"
#'   #     "LEU N"   "LEU O"   "LYS C"   "LYS CA"  "LYS CB"  "LYS CD"  "LYS CE"
#'   #     "LYS CG"  "LYS N"   "LYS NZ"  "LYS O"   "MET C"   "MET CA"  "MET CB"
#'   #     "MET CE"  "MET CG"  "MET N"   "MET O"   "MET SD"  "PHE C"   "PHE CA"
#'   #     "PHE CB"  "PHE CD1" "PHE CD2" "PHE CE1" "PHE CE2" "PHE CG"  "PHE CZ"
#'   #     "PHE N"   "PHE O"   "PRO C"   "PRO CA"  "PRO CB"  "PRO CD"  "PRO CG"
#'   #     "PRO N"   "PRO O"   "SER C"   "SER CA"  "SER CB"  "SER N"   "SER O"
#'   #     "SER OG"  "THR C"   "THR CA"  "THR CB"  "THR CG2" "THR N"   "THR O"
#'   #     "THR OG1" "TRP C"   "TRP CA"  "TRP CB"  "TRP CD1" "TRP CD2" "TRP CE2"
#'   #     "TRP CE3" "TRP CG"  "TRP CH2" "TRP CZ2" "TRP CZ3" "TRP N"   "TRP NE1"
#'   #     "TRP O"   "TYR C"   "TYR CA"  "TYR CB"  "TYR CD1" "TYR CD2" "TYR CE1"
#'   #     "TYR CE2" "TYR CG"  "TYR CZ"  "TYR N"   "TYR O"   "TYR OH"  "VAL C"
#'   #     "VAL CA"  "VAL CB"  "VAL CG1" "VAL CG2" "VAL N"   "VAL O"
#'
#' @export
#'
#' @family constants
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
names.res.AtomTypes <- sort(c("ALA N", "ALA CA", "ALA C", "ALA O", "ALA CB",

                              "ARG N", "ARG CA", "ARG C", "ARG O", "ARG CB",
                              "ARG CG", "ARG CD", "ARG NE", "ARG CZ", "ARG NH1",
                              "ARG NH2",

                              "ASN N", "ASN CA", "ASN C", "ASN O", "ASN CB",
                              "ASN CG", "ASN OD1", "ASN ND2",

                              "ASP N", "ASP CA", "ASP C", "ASP O", "ASP CB",
                              "ASP CG", "ASP OD1", "ASP OD2",

                              "CYS N", "CYS CA", "CYS C", "CYS O", "CYS CB",
                              "CYS SG",

                              "GLN N", "GLN CA", "GLN C", "GLN O", "GLN CB",
                              "GLN CG", "GLN CD", "GLN OE1", "GLN NE2",

                              "GLU N", "GLU CA", "GLU C", "GLU O", "GLU CB",
                              "GLU CG", "GLU CD", "GLU OE1", "GLU OE2",

                              "GLY N", "GLY CA", "GLY C", "GLY O",

                              "HIS N", "HIS CA", "HIS C", "HIS O", "HIS CB",
                              "HIS CG", "HIS CD2", "HIS ND1", "HIS CE1",
                              "HIS NE2",

                              "ILE N", "ILE CA", "ILE C", "ILE O", "ILE CB",
                              "ILE CG1", "ILE CG2", "ILE CD1",

                              "LEU N", "LEU CA", "LEU C", "LEU O", "LEU CB",
                              "LEU CG", "LEU CD1", "LEU CD2",

                              "LYS N", "LYS CA", "LYS C", "LYS O", "LYS CB",
                              "LYS CG", "LYS CD", "LYS CE", "LYS NZ",

                              "MET N", "MET CA", "MET C", "MET O", "MET CB",
                              "MET CG", "MET SD", "MET CE",

                              "PHE N", "PHE CA", "PHE C", "PHE O", "PHE CB",
                              "PHE CG", "PHE CD1", "PHE CD2", "PHE CE1",
                              "PHE CE2", "PHE CZ",

                              "PRO N", "PRO CA", "PRO C", "PRO O", "PRO CB",
                              "PRO CG", "PRO CD",

                              "SER N", "SER CA", "SER C", "SER O", "SER CB",
                              "SER OG",

                              "THR N", "THR CA", "THR C", "THR O", "THR CB",
                              "THR OG1", "THR CG2",

                              "TRP N", "TRP CA", "TRP C", "TRP O", "TRP CB",
                              "TRP CG", "TRP CD1", "TRP CD2", "TRP NE1",
                              "TRP CE2", "TRP CE3", "TRP CZ2", "TRP CZ3",
                              "TRP CH2",

                              "TYR N", "TYR CA", "TYR C", "TYR O", "TYR CB",
                              "TYR CG", "TYR CD1", "TYR CD2",
                              "TYR CE1", "TYR CE2", "TYR CZ", "TYR OH",

                              "VAL N", "VAL CA", "VAL C", "VAL O", "VAL CB",
                              "VAL CG1", "VAL CG2"))


## water residue names docs ----------------------------------------------------
#' @title Water Residue Names
#' @description Water residue names based on PDB naming conventions.
#' @details The three (3) letter abbreviation for the three (3) commonly used
#'   abbreviations for water residues.
#'
#' @name names.waters
#'
#' @examples
#'   names.waters
#'   # [1] "HOH", "DOD", "WAT"
#'
#' @export
#'
#' @family constants
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
names.waters <- c("HOH", "DOD", "WAT")


## openxlsx Cell Style docs ----------------------------------------------------
#' @title openxlsx Cell Style
#' @description A collection of cell style formats for the [openxlsx]
#'   package.
#' @details A centralized location defining the cell styles removes the need to
#'   change the formatting in several functions and provides a way to
#'   standardize cell formatting throughout the results.
#'
#'   The cell styles for the [openxlsx] package are defined within the
#'   `openxlsxCellStyle.R` file.
#'
#'   The defined cell styles are:
#'   * **cs.green**: background: lime, font: green and bold
#'   * **cs.pink**: background: pink, font: red and bold
#'   * **cs.amber**: background: amber, font orange and bold
#'   * **cs.0digits**: integer?
#'   * **cs.comma**: comma delineated values; _e.g._, `1,234`
#'   * **cs.date**: date formatted
#'   * **cs.1digits**: one digit after the decimal point
#'   * **cs.2digits**: two digits after the decimal point
#'   * **cs.3digits**: three digits after the decimal point
#'   * **cs.4digits**: four digits after the decimal point
#'   * **cs.header**: top row of table; font: black, bold, centered, with a
#'     line along the bottom of the cell
#'   * **cs.titles.tables**: top row of table; font: black, bold, and centered
#'
#' @name openxlsxCellStyles
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "openxlsx functions"
#'
NULL
##_ cell colors -----
cs.green <- openxlsx::createStyle(textDecoration = "bold",
                                  fontColour = "#006100",
                                  bgFill = "#c6efce",
                                  fgFill = "#c6efce")


cs.pink <- openxlsx::createStyle(textDecoration = "bold",
                                 fontColour = "#9c0006",
                                 bgFill = "#ffc7ce",
                                 fgFill = "#ffc7ce")


cs.amber <- openxlsx::createStyle(textDecoration = "bold",
                                  fontColour = "darkgoldenrod",
                                  bgFill = "#ffed98",
                                  fgFill = "#ffed98")

##_ number of digits -----
cs.comma <- openxlsx::createStyle(numFmt = "COMMA")

cs.date <- openxlsx::createStyle(numFmt = "DATE")

cs.0digits <- openxlsx::createStyle(numFmt = "GENERAL")

cs.1digits <- openxlsx::createStyle(numFmt = paste0("0.", paste0(rep(0, 1),
                                                                 collapse = ""))
                                    )

cs.2digits <- openxlsx::createStyle(numFmt = paste0("0.", paste0(rep(0, 2),
                                                                 collapse = ""))
                                    )

cs.3digits <- openxlsx::createStyle(numFmt = paste0("0.", paste0(rep(0, 3),
                                                                 collapse = ""))
                                    )

cs.4digits <- openxlsx::createStyle(numFmt = paste0("0.", paste0(rep(0, 4),
                                                                 collapse = ""))
                                    )

##_ headers and titles -----
cs.header <- openxlsx::createStyle(textDecoration = "bold",
                                   halign = "center",
                                   border = "Bottom")

cs.titles.tables <- openxlsx::createStyle(textDecoration = "bold",
                                          halign = "center")


## color values for plots ------------------------------------------------------
#' @title Color Values for Plots
#' @description Color values for plots with percent waters conserved plots.
#' @details The five (5) and six (6) color palettes are to used to color-code
#'   the plots illustrating percent water conserved (water conservation). The
#'   five color palette is for conservation values between 50\% to 100\% and the
#'   six color palette includes a color for less than 50\% conservation.
#'
#'   The colors are based on "percent conservation" with light grey dots
#'   indicating clusters with less than 50\% conservation, dark red dots
#'   representing clusters with 50\% to 69\% conservations, red dots are
#'   clusters with 70\% to 79\% conservation, light blue dots have 80\% to 89\%
#'   conservation, blue dots are clusters with 90\% to 99\% conservation, and
#'   dark blue dots are 100\% conserved water clusters (all structures
#'   contribute to the water cluster).
#'
#'   The defined colors are:
#'   - **cons.color5**: red, medium red, light blue, medium blue, and dark blue
#'   - **cons.color6**: light grey, red, medium red, light blue, medium blue,
#'     and dark blue
#'
#'   The defined legend titles are:
#'   - **cons.color5.legend**: Water Conservation
#'   - **cons.color6.legend**: Water Conservation
#'
#'   The defined break titles are:
#'   - **cons.color5.breaks**: set1, set2, set3, set4, and set5
#'   - **cons.color6.breaks**: set0, set1, set2, set3, set4, and set5
#'
#'   The defined labels are:
#'   - **cons.color5.labels**: 50-69\%, 70-79\%, 80-89\%, 90-99\%, 100\%
#'   - **cons.color6.labels**: < 50\%, 50-69\%, 70-79\%, 80-89\%, 90-99\%, 100\%
#'
#' @name colorPalettes
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "vanddraabe plots"
#'
#' @references
#'   Cynthia A Brewer. 200x.
#'   [http://www.ColorBrewer.org](http://www.ColorBrewer.org), accessed
#'   March 27, 2017.
#'
NULL
##_ five (5) colors (red to dark blue) -----
cons.color5 <- c("#a50f15", ## red
                 "#ef3b2c", ## medium red
                 "#9ecae1", ## light blue
                 "#4292c6", ## medium blue
                 "#08519c"  ## dark blue
)
##__ legend information -----
cons.color5.legend  <- "Water Conservation"
cons.color5.breaks <- c("set1", "set2", "set3", "set4", "set5")
cons.color5.labels <- c("50-69%", "70-79%", "80-89%", "90-99%", "100%")

##_ six (6) colors (red to dark blue) -----
cons.color6 <- c("#d9d9d9", ## light grey
                 "#a50f15", ## red
                 "#ef3b2c", ## medium red
                 "#9ecae1", ## light blue
                 "#4292c6", ## medium blue
                 "#08519c"  ## dark blue
)
##__ legend information -----
cons.color6.legend  <- "Water Conservation"
cons.color6.breaks <- c("set0", "set1", "set2", "set3", "set4", "set5")
cons.color6.labels <- c("< 50%", "50-69%", "70-79%", "80-89%", "90-99%", "100%")
