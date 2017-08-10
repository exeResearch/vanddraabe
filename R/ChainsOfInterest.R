## ChainsOfInterest.R
##
## oct-07-2016 (exe) created
## oct-07-2016 (exe) RetainChainsOfInterest accepts pdb$atoms instead of pdb
## jul-25-2017 (exe) updated documentation
## jul-31-2017 (exe) updated RetainChainsOfInterest() documentation
## jul-31-2017 (exe) changed file name: ChainsOfInterest_functions.R -> ChainsOfInterest.R
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




## DetermineChainsOfInterest docs ----------------------------------------------
#' @title Determine Chains Of Interest
#' @description Determine the chains identification
#' @details Standardizes user provided chain(s) of interest. This function
#'   simply standardizes the user provided chains of interest. Acceptable values
#'   are: - **first**: alphabetically the first chain - **all**: all chains
#'   within a structure file - **user defined**: a single letter or a set of
#'   letters; _e.g._; `"A"` or `c("H", "L")`
#'
#'   _**NOTE**_: This is a _**non-public**_ function and is _**NOT**_ available
#'   for general use. Please contact the author if you believe this function
#'   should be available for general use.
#'
#' @param chains.to.explore _**NOTE**_: `"first"` is alphabetically first. Thus
#'   if the order within the original PDB file is `L` and then `H`, this
#'   function will return `H` because it is alphebetically first.
#'
#' @return string indicating which chain designation (_e.g._, `"first"` chain,
#'   `"all"` chains, or `"user"` defined) to include in the conserved water
#'   analysis
#'
#' @examples
#'   \dontrun{
#'   DetermineChainsOfInterest("first")
#'   # [1] "first"
#'   DetermineChainsOfInterest("ALL")
#'   # [1] "all"
#'   DetermineChainsOfInterest("D")
#'   # [1] "user"
#'   DetermineChainsOfInterest(c("H", "L"))
#'   # [1] "user"
#'   DetermineChainsOfInterest("vanddraabe")
#'   # The provided chain ID VANDDRAABE is not valid and the first chain will
#'   # be used; likely chain A.
#'   # [1] "first"
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family utilities
#'
DetermineChainsOfInterest <- function(chains.to.explore) {

  ##----- convert provided chains to upper case
  chains.explore <- toupper(chains.to.explore)

  ##----- check if provided is possible choice
  possible.chains <- c("FIRST", "ALL", LETTERS)
  if ( any(possible.chains == chains.explore) == FALSE ) {
    mess <- paste("The provided chain ID ", chains.explore, " is not valid ",
                  "and the first chain will be used; likely chain A.", sep = "")
    chains.explore <- "FIRST"
    message(mess)
  }


  ##------ determine the chains of interest classification
  if ( length(chains.explore) == 1 ) {
    if ( chains.explore == "FIRST" ) { chains.oi <- "first" }
    if ( chains.explore == "ALL" )   { chains.oi <- "all" }
    if ( (chains.explore != "FIRST") & (chains.explore != "ALL") ) {
      chains.oi <- "user"
      }
  } else {
    chains.oi <- "user"
  }

  ##------ return chains of interest
  return(chains.oi)

}


## RetainChainsOfInterest docs -------------------------------------------------
#' @title Retain Chains Of Interest
#' @description Retain chains of interest based on user input parameters
#' @details Using the user provided chains of interest, indicate the PDB chains
#'   to retain.
#'
#'   _**NOTE**_: This is a _**non-public**_ function and is _**NOT**_ available
#'   for general use. Please contact the author if you believe this function
#'   should be available for general use.
#'
#' @param atoms.oi `data.frame` containing the PDB of the protein
#' @param chains.explore String of chains to explore
#' @param chains.oi String from the [DetermineChainsOfInterest()] function
#'   indicating if `"first"`, `"all"`, or a `"user"` defined set of chains
#'   should be used. _**NOTE**_: `"first"` is alphabetically first. Thus if the
#'   order within the original PDB file is `L` and then `H`, this function will
#'   return `H` because it is alphebetically first.
#'
#' @return `data.frame` of the protein atoms retained based on the indicated
#'   chains of interest
#'
#' @examples
#'   \dontrun{
#'   RetainChainsOfInterest(PDB.4dfr$atom, "B", "user")
#'   RetainChainsOfInterest(PDB.1hai$atom, c("H", "L"), "user")
#'   RetainChainsOfInterest(PDB.4dfr$atom, "A", "first")
#'   }
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family utilities
#'
RetainChainsOfInterest <- function(atoms.oi, chains.explore, chains.oi) {

  ##----- determine the chains within the structure
  atoms.chains <- atoms.oi$chain
  atoms.chains.uniq <- unique(atoms.chains)

  ##----- first chain (alphabetically)
  if ( chains.oi == "first") {
    chains.oi.letters <- sort(atoms.chains.uniq)[1]
  }

  ##----- all chains
  if ( chains.oi == "all") {
    chains.oi.letters <- atoms.chains.uniq
  }

  ##----- user defined chains
  if ( chains.oi == "user" ) {
    ##--- check if the user defined chains are available
    chains.oi.letters <- chains.explore[chains.explore %in% atoms.chains.uniq]
    ##--- indicated if the structure does NOT have all the chains of interest
    if ( sum(chains.explore %in% atoms.chains.uniq) < length(chains.explore) ) {
      message("Some of the requested chains are NOT present in this structure.")
    }
  }

  ##----- get the atoms for the chains of interest
  atoms.chains.oi <- atoms.oi[atoms.chains %in% chains.oi.letters, ]

  ##----- return the atoms for the chain of interest
  return(atoms.chains.oi)

}





