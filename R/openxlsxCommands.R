## openxlsxCommands.R :: the molecular similarity within a dataset
##
## NOTE!!!! THESE FUNCTIONS ARE HIGHLY SPECIFIC AND ARE MEANT TO BE USED WITH
##          THE INDICATED FUNCTION(S) AND data.frames.
##
## apr-15-2016 (exe) created
## apr-25-2016 (exe) added the AlignOverlap worksheet
## apr-25-2016 (exe) updated documentation
## feb-09-2017 (exe) updated documentation (rMarkdown)
## mar-15-2017 (exe) updated documentation (rMarkdown)
## apr-14-2017 (exe) added the Cleaned PDB Summary worksheet for
##                   CleanProteinStructures()
## apr-14-2017 (exe) added oxPlainDataSheet() to construct a plain worksheet
## jul-25-2017 (exe) updated documentation
## jul-31-2017 (exe) updated oxPDBcleanedSummarySheet() and oxPlainDataSheet()
##                   documentation
## jul-31-2017 (exe) updated documentation to indicate cell formatting
## aug-02-2017 (exe) updated documentation
## aug-07-2017 (exe) added cell formatting for "*.calc" to oxClusterSummarySheet()
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


## oxClusterStatsSheet docs ---------------------------------------------------------
#' @title openxlsx Water Cluster Statistics
#' @description Constructs the [openxlsx] worksheet for the Water
#'   Cluster statistics.
#' @details **This function is to _ONLY_ be used with the results of
#'   [ConservedWaterStats()]**. Specific aspects of how the
#'   returned `data.frame` will be formatted are **hard-coded** into this
#'   function.
#'
#'   This [openxlsx] function is _**NOT**_ exported.
#'
#' @param wb.name Name of the workbook for the results; _e.g._, results.wb
#' @param sheet.name Name of the worksheet being formatted; default:
#'   "ClusterStatistics"
#' @param df data.frame containing the results of the `GetSimilarityPairs`
#'   function; _e.g._, `h2o.cluster.stats`
#'
#' @return The workbook containing the indicated and newly formatted worksheet.
#'
#' @import openxlsx
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "openxlsx functions"
#'
oxClusterStatsSheet <- function(wb.name,
                                sheet.name = "ClusterStatistics",
                                df) {

  ##----- add sheet
  openxlsx::addWorksheet(wb = wb.name, sheetName = sheet.name, gridLines = TRUE)

  ##----- add data
  openxlsx::writeData(wb = wb.name,
                      sheet = sheet.name,
                      x = df,
                      startCol = 1,
                      startRow = 3,
                      colNames = FALSE,
                      rowNames = TRUE,
                      keepNA = FALSE,
                      withFilter = FALSE,
                      headerStyle = NULL)

  ##----- merge cells for titles
  openxlsx::mergeCells(wb = wb.name,
                       sheet = sheet.name,
                       cols = 2:3,
                       rows = 1)
  openxlsx::mergeCells(wb = wb.name,
                       sheet = sheet.name,
                       cols = 4:5,
                       rows = 1)

  ##----- add table titles
  openxlsx::writeData(wb = wb.name,
                      sheet = sheet.name,
                      startCol = 2,
                      startRow = 1,
                      x = "All Waters")
  openxlsx::writeData(wb = wb.name,
                      sheet = sheet.name,
                      startCol = 4,
                      startRow = 1,
                      x = "Passed Waters")

  openxlsx::writeData(wb = wb.name,
                      sheet = sheet.name,
                      startCol = 2,
                      startRow = 2,
                      x = "Number")
  openxlsx::writeData(wb = wb.name,
                      sheet = sheet.name,
                      startCol = 3,
                      startRow = 2,
                      x = "Percentage (%)")
  openxlsx::writeData(wb = wb.name,
                      sheet = sheet.name,
                      startCol = 4,
                      startRow = 2,
                      x = "Number")
  openxlsx::writeData(wb = wb.name,
                      sheet = sheet.name,
                      startCol = 5,
                      startRow = 2,
                      x = "Percentage (%)")

  ##----- set column widths
  openxlsx::setColWidths(wb = wb.name,
                         sheet = sheet.name,
                         cols = 1:5,
                         widths = "auto")

  ##----- number of digits for "Average conservation"
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.2digits,
                     cols = c(2, 4),
                     rows = 7,
                     gridExpand = TRUE,
                     stack = TRUE)

  ##----- number of digits for percentage of clusters with various proportions
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.titles.tables,
                     cols = c(2:5),
                     rows = c(1:2),
                     gridExpand = TRUE,
                     stack = TRUE)

  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.1digits,
                     cols = c(3, 5),
                     rows = c(8:13),
                     gridExpand = TRUE,
                     stack = TRUE)


  ##----- number of digits for calculation times
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.comma,
                     cols = c(2, 4),
                     rows = c(3:6,14),
                     gridExpand = TRUE,
                     stack = TRUE)

  ##----- number of digits for calculation times
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.4digits,
                     cols = c(2, 4),
                     rows = c(16:17),
                     gridExpand = TRUE,
                     stack = TRUE)

  ##----- return the workbook
  return(wb.name)

}


## oxWaterOccurrenceSheet docs -------------------------------------------------------
#' @title openxlsx Water Occurrence Summary
#' @description Constructs the [openxlsx] worksheet for the Water Occurrence
#'   summary.
#' @details **This function is to _ONLY_ be used with the results of
#'   [ConservedWaters()]**. Specific aspects of how the
#'   returned `data.frame` will be formatted are **hard-coded** into this
#'   function.
#'
#'   This [openxlsx] function is _**NOT**_ exported.
#'
#' @param wb.name Name of the workbook for the results; _e.g._, results.wb
#' @param sheet.name Name of the worksheet being formatted; default:
#'   `"WaterOccurrenceSummary"`
#' @param df data.frame containing the water occurrence results of the
#'   [ConservedWaters()] function; _e.g._, `h2o.occurrence`
#'
#' @return The workbook containing the indicated and newly formatted worksheet.
#'
#' @import openxlsx
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "openxlsx functions"
#'
oxWaterOccurrenceSheet <- function(wb.name,
                                   sheet.name = "WaterOccurrenceSummary",
                                   df) {

  ##----- data.frame dimensions and other information
  df.size <- dim(df)
  df.row.idc <- 2:(df.size[1] + 1)
  df.col.idc <- 1:(df.size[2] + 1)
  df.colnames <- colnames(df)

  ##----- add sheet
  openxlsx::addWorksheet(wb = wb.name,
                         sheetName = sheet.name,
                         gridLines = TRUE)

  ##----- add data
  openxlsx::writeDataTable(wb = wb.name,
                           sheet = sheet.name,
                           x = as.data.frame(df, stringsAsFactors = FALSE),
                           colNames = TRUE,
                           rowNames = TRUE,
                           tableStyle = "none",
                           withFilter = FALSE,
                           headerStyle = cs.header)

  ##----- freeze the panes
  openxlsx::freezePane(wb = wb.name,
                       sheet = sheet.name,
                       firstRow = TRUE, firstCol = TRUE)

  ##----- set column widths
  openxlsx::setColWidths(wb = wb.name,
                         sheet = sheet.name,
                         cols = df.col.idc,
                         widths = "auto")

  ##----- number of digits
  ##--- the resolution column
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.2digits,
                     cols = 2,
                     rows = df.row.idc,
                     gridExpand = TRUE,
                     stack = TRUE)
  ##--- the rObserved and rFree columns
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.3digits,
                     cols = c(3, 4),
                     rows = df.row.idc,
                     gridExpand = TRUE,
                     stack = TRUE)
  ##--- the means columns
  mu.cols.idc <- grep(pattern = ".mu", df.colnames) + 1
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.2digits,
                     cols = mu.cols.idc,
                     rows = df.row.idc,
                     gridExpand = TRUE,
                     stack = TRUE)
  ##--- the standard deviation columns
  sd.cols.idc <- grep(pattern = ".sd", df.colnames) + 1
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.3digits,
                     cols = sd.cols.idc,
                     rows = df.row.idc,
                     gridExpand = TRUE,
                     stack = TRUE)
  ##--- the percent
  pct.cols.idc <- grep(pattern = "pct.", df.colnames) + 1
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.1digits,
                     cols = pct.cols.idc,
                     rows = df.row.idc,
                     gridExpand = TRUE,
                     stack = TRUE)

  ##----- color-code the cluster true-false table
  cluster.cols.idc <- grep(pattern = "clust_", df.colnames) + 1
  openxlsx::conditionalFormatting(wb = wb.name,
                                  sheet = sheet.name,
                                  cols = cluster.cols.idc,
                                  rows = df.row.idc,
                                  type = "contains",
                                  rule = "TRUE",
                                  style = cs.green)
  openxlsx::conditionalFormatting(wb = wb.name,
                                  sheet = sheet.name,
                                  cols = cluster.cols.idc,
                                  rows = df.row.idc,
                                  type = "contains",
                                  rule = "FALSE",
                                  style = cs.pink)

  ##----- return the workbook
  return(wb.name)

}




## oxClusterSummarySheet docs --------------------------------------------------
#' @title openxlsx Cluster Summary Sheet
#' @description Constructs the [openxlsx] worksheet for the Cluster Summary
#'   analysis.
#' @details **This function is to _ONLY_ be used with the results of
#'   [ConservedWaters()]**. Specific aspects of how the
#'   returned `data.frame` will be formatted are **hard-coded** into this
#'   function.
#'
#'   This [openxlsx] function is _**NOT**_ exported.
#'
#' @param wb.name Name of the workbook for the results; _e.g._, results.wb
#' @param sheet.name Name of the worksheet being formatted; default:
#'   `"ClusterSummary"`
#' @param df data.frame containing the cluster summary from the
#'   [ConservedWaters()] function; _e.g._, `h2o.clusters.summary`
#'
#' @return The workbook containing the indicated and newly formatted worksheet.
#'
#' @import openxlsx
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "openxlsx functions"
#'
oxClusterSummarySheet <- function(wb.name,
                                  sheet.name = "ClusterSummary",
                                  df) {

  ##----- data.frame dimensions and other information
  df.size <- dim(df)
  df.row.idc <- 2:(df.size[1] + 1)
  df.col.idc <- 1:df.size[2]
  df.colnames <- colnames(df)


  ##----- add sheet
  openxlsx::addWorksheet(wb = wb.name,
                         sheetName = sheet.name,
                         gridLines = TRUE)


  ##----- add data
  openxlsx::writeDataTable(wb = wb.name,
                           sheet = sheet.name,
                           x = as.data.frame(df, stringsAsFactors = FALSE),
                           colNames = TRUE,
                           rowNames = FALSE,
                           tableStyle = "none",
                           withFilter = FALSE,
                           headerStyle = cs.header)


  ##----- freeze the panes
  openxlsx::freezePane(wb = wb.name,
                       sheet = sheet.name,
                       firstRow = FALSE, firstActiveRow = 2,
                       firstCol = FALSE, firstActiveCol = 4)


  ##----- set column widths
  openxlsx::setColWidths(wb = wb.name,
                         sheet = sheet.name,
                         cols = df.col.idc,
                         widths = "auto")


  ##----- number of digits
  ##--- the x, y, and z columns
  xyz.cols.idc <- which(df.colnames %in% c("x", "y", "z"))
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.3digits,
                     cols = xyz.cols.idc,
                     rows = df.row.idc,
                     gridExpand = TRUE,
                     stack = TRUE)
  ##--- the calc columns
  mu.cols.idc <- grep(pattern = ".calc", df.colnames)
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.2digits,
                     cols = mu.cols.idc,
                     rows = df.row.idc,
                     gridExpand = TRUE,
                     stack = TRUE)
  ##--- the means columns
  mu.cols.idc <- grep(pattern = ".mu", df.colnames)
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.2digits,
                     cols = mu.cols.idc,
                     rows = df.row.idc,
                     gridExpand = TRUE,
                     stack = TRUE)
  ##--- the standard deviation columns
  sd.cols.idc <- grep(pattern = ".sd", df.colnames)
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.3digits,
                     cols = sd.cols.idc,
                     rows = df.row.idc,
                     gridExpand = TRUE,
                     stack = TRUE)
  ##--- the percent
  pct.cols.idc <- grep(pattern = "pct.", df.colnames)
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.1digits,
                     cols = pct.cols.idc,
                     rows = df.row.idc,
                     gridExpand = TRUE,
                     stack = TRUE)
  ##--- the rmsf column
  rmsf.cols.idc <- grep(pattern = "rmsf", df.colnames)
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.4digits,
                     cols = rmsf.cols.idc,
                     rows = df.row.idc,
                     gridExpand = TRUE,
                     stack = TRUE)
  ##--- the means columns
  mobility.cols.idc <- grep(pattern = "mobility", df.colnames)
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.2digits,
                     cols = mobility.cols.idc,
                     rows = df.row.idc,
                     gridExpand = TRUE,
                     stack = TRUE)


  ##----- return the workbook
  return(wb.name)

  }


## oxRCSBinfoSheet docs --------------------------------------------------
#' @title openxlsx PDB/RCSB Summary Sheet
#' @description Constructs the [openxlsx] worksheet for the Similarity
#'   Summary analysis.
#' @details **This function is to _ONLY_ be used with the results of
#'   [ConservedWaters()]**. Specific aspects of how the
#'   returned `data.frame` will be formatted are **hard-coded** into this
#'   function.
#'
#'   This [openxlsx] function is _**NOT**_ exported.
#'
#' @param wb.name Name of the workbook for the results; _e.g._, results.wb
#' @param sheet.name Name of the worksheet being formatted; default:
#'   "PDB_information"
#' @param df data.frame containing the PDB/RCSB information obtatined within the
#'   [ConservedWaters()] function; _e.g._, `pdbs.information`
#'
#' @return The workbook containing the indicated and newly formatted worksheet.
#'
#' @import openxlsx
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "openxlsx functions"
#'
oxRCSBinfoSheet <- function(wb.name,
                            sheet.name = "RCSB_information",
                            df) {

  ##----- data.frame dimensions and other information
  df.size <- dim(df)
  df.row.idc <- 2:(df.size[1] + 1)
  df.col.idc <- 1:df.size[2]
  df.colnames <- colnames(df)

  ##----- add sheet
  openxlsx::addWorksheet(wb = wb.name,
                         sheetName = sheet.name,
                         gridLines = TRUE)

  ##----- add data
  openxlsx::writeDataTable(wb = wb.name,
                           sheet = sheet.name,
                           x = as.data.frame(df, stringsAsFactors = FALSE),
                           colNames = TRUE,
                           rowNames = FALSE,
                           tableStyle = "none",
                           withFilter = FALSE,
                           headerStyle = cs.header)

  ##----- freeze the panes
  openxlsx::freezePane(wb = wb.name,
                       sheet = sheet.name,
                       firstRow = FALSE, firstActiveRow = 2,
                       firstCol = FALSE, firstActiveCol = 5)

  ##----- set column widths
  openxlsx::setColWidths(wb = wb.name,
                         sheet = sheet.name,
                         cols = df.col.idc,
                         widths = "auto")

  ##----- set the date format
  depositionDate.idx <- which(df.colnames == "depositionDate")
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.date,
                     cols = depositionDate.idx,
                     rows = df.row.idc,
                     gridExpand = TRUE,
                     stack = TRUE)

  ##----- return the workbook
  return(wb.name)

}


## oxInitWaterDataSheet docs ---------------------------------------------------
#' @title Initial Water Data Sheet
#' @description Constructs the [openxlsx] worksheet for the initial water data.
#' @details **This function is to _ONLY_ be used with the results of
#'   [ConservedWaters()]**. Specific aspects of how the
#'   returned `data.frame` will be formatted are **hard-coded** into this
#'   function.
#'
#'   This [openxlsx] function is _**NOT**_ exported.
#'
#' @param wb.name Name of the workbook for the results; _e.g._, results.wb
#' @param sheet.name Name of the worksheet being formatted; default:
#'   `"InitialWaterData"`
#' @param df data.frame containing the concatenate initial waters with
#'   experimental and experimentally derived values obtatined within the
#'   [ConservedWaters()] function; _e.g._, `h2o.df`
#'
#' @return The workbook containing the indicated and newly formatted worksheet.
#'
#' @import openxlsx
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "openxlsx functions"
#'
oxInitWaterDataSheet <- function(wb.name,
                                 sheet.name = "InitialWaterData",
                                 df) {

  ##_ data.frame dimensions and other information -----
  df.size <- dim(df)
  df.row.idc <- 2:(df.size[1] + 1)
  df.col.idc <- 1:(df.size[2] + 1)
  df.colnames <- colnames(df)

  ##_ add sheet -----
  openxlsx::addWorksheet(wb = wb.name,
                         sheetName = sheet.name,
                         gridLines = TRUE)

  ##_ add data -----
  openxlsx::writeDataTable(wb = wb.name,
                           sheet = sheet.name,
                           x = as.data.frame(df, stringsAsFactors = FALSE),
                           colNames = TRUE,
                           rowNames = TRUE,
                           tableStyle = "none",
                           withFilter = FALSE,
                           headerStyle = cs.header)

  ##_ freeze the panes -----
  openxlsx::freezePane(wb = wb.name, sheet = sheet.name,
                       firstRow = TRUE, firstCol = TRUE)

  ##_ set column widths -----
  openxlsx::setColWidths(wb = wb.name, sheet = sheet.name,
                         cols = df.col.idc, widths = "auto")

  ##_ number of digits -----
  ##--- the x, y, and z columns
  xyz.cols.idc <- which(df.colnames %in% c("x", "y", "z")) + 1
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.3digits,
                     cols = xyz.cols.idc,
                     rows = df.row.idc,
                     gridExpand = TRUE,
                     stack = TRUE)
  ##--- atom number column
  eleno.cols.idc <- which(df.colnames %in% c("eleno")) + 1
  openxlsx::addStyle(wb = wb.name,
                     sheet = sheet.name,
                     style = cs.0digits,
                     cols = eleno.cols.idc,
                     rows = df.row.idc,
                     gridExpand = TRUE,
                     stack = TRUE)
  ##--- occupation, B-value, mobility, and normalized B-values column
  ob.cols.idc <- which(df.colnames %in% c("o","b","mobility","nBvalue")) + 1
  openxlsx::addStyle(wb = wb.name, sheet = sheet.name, style = cs.2digits,
                     cols = ob.cols.idc, rows = df.row.idc,
                     gridExpand = TRUE, stack = TRUE)
  ##--- the resolution column
  res.col.idc <- grep(pattern = "resolution", df.colnames) + 1
  openxlsx::addStyle(wb = wb.name, sheet = sheet.name, style = cs.2digits,
                     cols = res.col.idc, rows = df.row.idc,
                     gridExpand = TRUE, stack = TRUE)
  ##--- the rObserved and rFree columns
  res.col.idc <- c(grep(pattern = "rObserved", df.colnames),
                   grep(pattern = "rFree", df.colnames) ) + 1
  openxlsx::addStyle(wb = wb.name, sheet = sheet.name, style = cs.3digits,
                     cols = c(3, 4), rows = df.row.idc,
                     gridExpand = TRUE, stack = TRUE)
  ##--- the summation columns
  sum.cols.idc <- grep(pattern = ".sum", df.colnames) + 1
  openxlsx::addStyle(wb = wb.name, sheet = sheet.name, style = cs.3digits,
                     cols = sum.cols.idc, rows = df.row.idc,
                     gridExpand = TRUE, stack = TRUE)
  ##--- the means columns
  mu.cols.idc <- grep(pattern = ".mu", df.colnames) + 1
  openxlsx::addStyle(wb = wb.name, sheet = sheet.name, style = cs.2digits,
                     cols = mu.cols.idc, rows = df.row.idc,
                     gridExpand = TRUE, stack = TRUE)
  ##--- the standard deviation columns
  sd.cols.idc <- grep(pattern = ".sd", df.colnames) + 1
  openxlsx::addStyle(wb = wb.name, sheet = sheet.name, style = cs.3digits,
                     cols = sd.cols.idc, rows = df.row.idc,
                     gridExpand = TRUE, stack = TRUE)


  ##_ color-code the cluster true-false table -----
  tf.cols.idc <- c(grep(pattern = ".keep", df.colnames),
                   grep(pattern = ".cutoffs", df.colnames) ) + 1
  openxlsx::conditionalFormatting(wb = wb.name, sheet = sheet.name,
                                  cols = tf.cols.idc, rows = df.row.idc,
                                  type = "contains", rule = "TRUE",
                                  style = cs.green)
  openxlsx::conditionalFormatting(wb = wb.name, sheet = sheet.name,
                                  cols = tf.cols.idc, rows = df.row.idc,
                                  type = "contains", rule = "FALSE",
                                  style = cs.pink)


  ##_ return the workbook -----
  return(wb.name)

}




## oxAlignOverlapSheet docs ----------------------------------------------------
#' @title Align Overlap Data Sheet
#' @description Constructs the [openxlsx] worksheet for the
#'   [AlignOverlap()] results.
#' @details **This function is to _ONLY_ be used with the results of
#'   [AlignOverlap()]**. Specific aspects of how the
#'   returned `data.frame` will be formatted are **hard-coded** into this
#'   function.
#'
#'   Notable formatting:
#'   - Top row frozen
#'   - Column widths are set based on column content
#'   - Structures passing the [AlignOverlap()] evaluation are hightlighted
#'     lime green
#'   - Structures failing the [AlignOverlap()] evaluation are hightlighted
#'     pink
#'
#'   This [openxlsx] function is _**NOT**_ exported.
#'
#' @param wb.name Name of the workbook for the results; _e.g._, results.wb
#' @param sheet.name Name of the worksheet being formatted; default:
#'   `"AlignOverlap"`
#' @param df data.frame containing the summary of [AlignOverlap()]; _e.g._,
#'   `df.results`
#'
#' @return The workbook containing the indicated and newly formatted worksheet.
#'
#' @import openxlsx
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "openxlsx functions"
#'
oxAlignOverlapSheet <- function(wb.name,
                                sheet.name = "AlignOverlap",
                                df) {

  ##----- data.frame dimensions and other information
  df.size <- dim(df)
  df.row.idc <- 2:(df.size[1] + 1)
  df.col.idc <- 1:(df.size[2] + 1)
  df.colnames <- colnames(df)

  ##----- add sheet
  openxlsx::addWorksheet(wb = wb.name,
                         sheetName = sheet.name,
                         gridLines = TRUE)

  ##----- add data
  openxlsx::writeDataTable(wb = wb.name,
                           sheet = sheet.name,
                           x = as.data.frame(df, stringsAsFactors = FALSE),
                           colNames = TRUE,
                           rowNames = FALSE,
                           tableStyle = "none",
                           withFilter = FALSE,
                           headerStyle = cs.header)

  ##----- freeze the panes
  openxlsx::freezePane(wb = wb.name,
                       sheet = sheet.name,
                       firstRow = TRUE,
                       firstCol = FALSE)

  ##----- set column widths
  openxlsx::setColWidths(wb = wb.name,
                         sheet = sheet.name,
                         cols = df.col.idc,
                         widths = "auto")

  ##----- color-code the cluster true-false table
  openxlsx::conditionalFormatting(wb = wb.name,
                                  sheet = sheet.name,
                                  cols = c(2, 5),
                                  rows = df.row.idc,
                                  type = "contains",
                                  rule = "TRUE",
                                  style = cs.green)
  openxlsx::conditionalFormatting(wb = wb.name,
                                  sheet = sheet.name,
                                  cols = c(2, 5),
                                  rows = df.row.idc,
                                  type = "contains",
                                  rule = "FALSE",
                                  style = cs.pink)

  ##----- return the workbook
  return(wb.name)

}


## oxPDBcleanedSummarySheet docs -----------------------------------------------
#' @title Cleaned PDB Structures Data Sheet
#' @description Constructs the [openxlsx] worksheet for the
#'   [CleanProteinStructures()] results.
#' @details **This function is to _ONLY_ be used with the results of
#'   [CleanProteinStructures()]**. Specific aspects of how the
#'   returned `data.frame` will be formatted are **hard-coded** into this
#'   function.
#'
#'   Notable formatting:
#'   - Top row frozen
#'   - Column widths are set based on column content
#'   - Structures with hydrogen atoms removed are highlighted with amber cell color
#'   - Structures with OoR values, modeled atoms, and removed waters are highlighted
#'     with amber cell color
#'
#'   This [openxlsx] function is _**NOT**_ exported.
#'
#' @param wb.name Name of the workbook for the results; _e.g._, results.wb
#' @param sheet.name Name of the worksheet being formatted; default:
#'   `"PDBcleanedSummary"`
#' @param df data.frame containing the summary of
#'   [CleanProteinStructures()]; _e.g._, `df.results`
#'
#' @return The workbook containing the indicated and newly formatted worksheet.
#'
#' @import openxlsx
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "openxlsx functions"
#'
oxPDBcleanedSummarySheet <- function(wb.name,
                                sheet.name = "PDBcleanedSummary",
                                df) {

  ##----- data.frame dimensions and other information
  df.size <- dim(df)
  df.row.idc <- 2:(df.size[1] + 1)
  df.col.idc <- 1:(df.size[2] + 1)
  df.colnames <- colnames(df)

  ##----- add sheet
  openxlsx::addWorksheet(wb = wb.name,
                         sheetName = sheet.name,
                         gridLines = TRUE)

  ##----- add data
  openxlsx::writeDataTable(wb = wb.name,
                           sheet = sheet.name,
                           x = as.data.frame(df, stringsAsFactors = FALSE),
                           colNames = TRUE,
                           rowNames = FALSE,
                           tableStyle = "none",
                           withFilter = FALSE,
                           headerStyle = cs.header)

  ##----- freeze the panes
  openxlsx::freezePane(wb = wb.name,
                       sheet = sheet.name,
                       firstRow = TRUE,
                       firstCol = FALSE)

  ##----- set column widths
  openxlsx::setColWidths(wb = wb.name,
                         sheet = sheet.name,
                         cols = df.col.idc,
                         widths = "auto")

  ##----- color structures with removed hydrogens amber
  openxlsx::conditionalFormatting(wb = wb.name,
                                  sheet = sheet.name,
                                  cols = c(2),
                                  rows = df.row.idc,
                                  type = "contains",
                                  rule = "TRUE",
                                  style = cs.amber)

  ##----- color values for structures with OoR values, modeled atoms, and
  ##      removed waters
  openxlsx::conditionalFormatting(wb = wb.name,
                                  sheet = sheet.name,
                                  cols = c(3, 4, 5),
                                  rows = df.row.idc,
                                  type = "expression",
                                  rule = "!=0",
                                  style = cs.amber)
  openxlsx::conditionalFormatting(wb = wb.name,
                                  sheet = sheet.name,
                                  cols = c(7),
                                  rows = df.row.idc,
                                  type = "expression",
                                  rule = "!=0",
                                  style = cs.amber)

  ##----- return the workbook
  return(wb.name)

}


## oxPlainDataSheet docs -------------------------------------------------------
#' @title Plain Data Sheet
#' @description Constructs a plain Excel worksheet via the [openxlsx] package.
#' @details **This function creates a basic Excel worksheet with minimal
#'   formatting.**
#'
#'   Notable formatting:
#'   - Top row frozen
#'   - Column widths are set based on column content
#'
#'   This [openxlsx] function is _**NOT**_ exported.
#'
#' @param wb.name Name of the workbook for the results; _e.g._, results.wb
#' @param sheet.name Name of the worksheet being formatted; default:
#'   `"basic"`
#' @param df data.frame containing the data to be written; _e.g._, `df.results`
#'
#' @return The workbook containing the indicated and newly formatted worksheet.
#'
#' @import openxlsx
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @family "openxlsx functions"
#'
oxPlainDataSheet <- function(wb.name,
                             sheet.name = "basic",
                             df) {

  ##----- data.frame dimensions and other information
  df.size <- dim(df)
  df.row.idc <- 2:(df.size[1] + 1)
  df.col.idc <- 1:(df.size[2] + 1)
  df.colnames <- colnames(df)

  ##----- add sheet
  openxlsx::addWorksheet(wb = wb.name,
                         sheetName = sheet.name,
                         gridLines = TRUE)

  ##----- add data
  openxlsx::writeDataTable(wb = wb.name,
                           sheet = sheet.name,
                           x = as.data.frame(df, stringsAsFactors = FALSE),
                           colNames = TRUE,
                           rowNames = FALSE,
                           tableStyle = "none",
                           withFilter = FALSE,
                           headerStyle = cs.header)

  ##----- freeze the panes
  openxlsx::freezePane(wb = wb.name,
                       sheet = sheet.name,
                       firstRow = TRUE,
                       firstCol = TRUE)

  ##----- set column widths
  openxlsx::setColWidths(wb = wb.name,
                         sheet = sheet.name,
                         cols = df.col.idc,
                         widths = "auto")

  ##----- return the workbook
  return(wb.name)

}


