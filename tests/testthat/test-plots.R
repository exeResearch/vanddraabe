## test-plots.R
##
## jul-07-2017 (exe) created
## jul-25-2017 (exe) updated documentation
## jul-26-2017 (exe) finish the tests for plots
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


context("vanddraabe plots")

## load data file --------------------------------------------------------------
load(file="Thrombin10_data.rda")


## OccupancyBarplot.summ -------------------------------------------------------
test_that("the Occupancy Summary Barplots (initial structures) are the same", {
  expect_equal_to_reference(OccupancyBarplot.summ(data=thrombin10.CLEANED),
                            file="plot-test_OccupancyBarplot_SUMMARY.rds")
})


## BvalueBarplot.summ ----------------------------------------------------------
test_that("the b-value Summary Barplots (initial structures) are the same", {
  expect_equal_to_reference(BvalueBarplot.summ(data=thrombin10.CLEANED),
                            file="plot-test_BvalueBarplot_SUMMARY.rds")
})


## normBvalueBarplot.summ ------------------------------------------------------
test_that("the Normalized b-value Summary Barplots (initial structures) are the same", {
  expect_equal_to_reference(normBvalueBarplot.summ(data=thrombin10.CLEANED),
                            file="plot-test_normBvalueBarplot_SUMMARY.rds")
})


## MobilityBarplot.summ --------------------------------------------------------
test_that("the Mobility Summary Barplots (initial structures) are the same", {
  expect_equal_to_reference(MobilityBarplot.summ(data=thrombin10.CLEANED),
                            file="plot-test_MobilityBarplot_SUMMARY.rds")
})


## ConservationPlot ------------------------------------------------------------
# test_that("the Conservation Plot is the same", {
#   expect_equal_to_reference(ConservationPlot(data=thrombin10.conservedWaters,
#                                              passed.waters=TRUE),
#                             file="plot-test_ConservationPlot.rds")
# })


## OccupancyBarplot ----------------------------------------------------------
# test_that("the Occupancy Barplot is the same", {
#   expect_equal_to_reference(OccupancyBarplot(data=thrombin10.conservedWaters,
#                                              passed.waters=TRUE),
#                             file="plot-test_OccupancyBarplot.rds")
# })


## MobilityBarplot -------------------------------------------------------------
# test_that("the Mobility Barplot is the same", {
#   expect_equal_to_reference(MobilityBarplot(data=thrombin10.conservedWaters,
#                                             passed.waters=TRUE),
#                             file="plot-test_MobilityBarplot.rds")
# })


## BvalueBarplot ---------------------------------------------------------------
# test_that("the b-value Barplot is the same", {
#   expect_equal_to_reference(BvalueBarplot(data=thrombin10.conservedWaters,
#                                           passed.waters=TRUE, calc.values=TRUE),
#                             file="plot-test_BvalueBarplot_calcValuesTRUE.rds")
#   expect_equal_to_reference(BvalueBarplot(data=thrombin10.conservedWaters,
#                                           passed.waters=TRUE, calc.values=FALSE),
#                             file="plot-test_BvalueBarplot_calcValuesFALSE.rds")
# })


## nBvalueBarplot --------------------------------------------------------------
# test_that("the normalized b-value Barplot is the same", {
#   expect_equal_to_reference(nBvalueBarplot(data=thrombin10.conservedWaters,
#                                            passed.waters=TRUE),
#                             file="plot-test_nBvalueBarplot.rds")
# })


## ClusterSummaryPlots ----------------------------------------------------------
test_that("the Cluster Summary Plots is the same", {
  expect_equal_to_reference(ClusterSummaryPlots(data=thrombin10.conservedWaters,
                                               passed.waters=TRUE),
                            file="plot-test_ClusterSummaryPlots.rds")
})


## MobNormBvalEvalPlots --------------------------------------------------------
test_that("the Mobility and Normalized b-value Plots are the same", {
  expect_equal_to_reference(MobNormBvalEvalPlots(data=thrombin10.conservedWaters,
                                                 passed.waters=TRUE),
                            file="plot-test_MobNormBvalEvalPlots.rds")
})


## BoundWaterEnvPlots ----------------------------------------------------------
test_that("the Bound Water Environment Plots are the same", {
  expect_equal_to_reference(BoundWaterEnvPlots(data=thrombin10.conservedWaters,
                                               passed.waters=TRUE),
                            file="plot-test_BoundWaterEnvPlots.rds")
})


## BoundWaterEnvSummaryPlot ----------------------------------------------------
test_that("the Bound Water Environment Summary Plot is the same", {
  expect_equal_to_reference(BoundWaterEnvSummaryPlot(data=thrombin10.conservedWaters,
                                                     passed.waters=TRUE),
                            file="plot-test_BoundWaterEnvSummaryPlot.rds")
})
