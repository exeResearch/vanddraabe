## ClusterWaters.R
##
## dec-28-2015 (exe) created
## may-05-2016 (exe) ConservedWaters: copy PDB files to used and rejected
##                   directories
## may-05-2016 (exe) ConservedWaters: rewrote PDB quality evaluations
## jul-30-2016 (exe) ConservedWaters: remove modeled heavy atoms
## sep-16-2016 (exe) split original file to create ClusterWaters.R
## apr-26-2017 (exe) added ClusterWaters.MDS() fxn
## apr-26-2017 (exe) updated documentation
## jul-25-2017 (exe) updated documentation
## aug-03-2017 (exe) added teh check.cluster.method() function
## aug-03-2017 (exe) calculate B-value from RMSF of waters in a cluster for ClusterWaters()
## aug-08-2017 (exe) replaced stats::hclust with fastcluster::hclust (see fxn documentation below)
## aug-09-2017 (exe) added @importFrom stats ... and utils ...
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



## check.cluster.method docs ---------------------------------------------------
#' @title Check Clustering Method
#' @description Ensures the user provided clustering method is a valid choice.
#' @details A simple check and reformatting of the clustering method indicated
#'   by the user in the [ConservedWaters()] and [ConservedWaters.MDS()]
#'   parameters.
#'
#' @param cluster.method The user defined clustering method for the
#'   [ConservedWaters()] and [ConservedWaters.MDS()].
#'
#' @return Correctly formatted clustering method or a stop error
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
check.cluster.method <- function(cluster.method) {

  cluster.method <- tolower(cluster.method)
  cluster.methods.avail <- c("ward.d", "ward.d2", "single", "complete")

  ##_ check of cluster.method is available -----
  if ( !any(cluster.method == cluster.methods.avail) ) {
    cluster.stop <- paste("Please select one of the following clustering",
                          "methods: \"complete\" (the default), \"ward.D\",",
                          "\"ward.D2\", or \"single\".", sep = " ")
    stop(cluster.stop)
  }

  ##_ correct case for ward methods -----
  if ( cluster.method == "ward.d" ) cluster.method <- "ward.D"
  if ( cluster.method == "ward.d2" ) cluster.method <- "ward.D2"

  ##_ return cluster.method -----
  return(cluster.method)

}


## cluster waters docs ---------------------------------------------------------
#' @title Cluster Conserved Waters
#' @description Cluster the conserved waters.
#' @details Calculate the conserved waters using a collection of
#'   crystallographic protein structures.
#'
#' @param data The water oxygens' X, Y, and Z coordinates, B-values, and
#'   occupancy values.
#' @param cutoff.cluster Numerical value provided by the user for the distance
#'   between water oxygen atoms to form a cluster; default: 2.4 Angstroms.
#' @param cluster.method Method of clustering the waters; default is
#'   "`complete`". Any other method accepted by the `hclust` function is
#'   appropriate. The original method used by Sanschagrin and Kuhn is the
#'   complete linkage clustering method and is the default. Other options
#'   include "`ward.D`" (equivilant to the only Ward option in `R` versions
#'   3.0.3 and earlier), "`ward.D2`" (implements Ward's 1963 criteria; see
#'   Murtagh and Legendre 2014), "`single`" (related to the minimal spanning
#'   tree method and adopts a "friend of friends" clustering method), along with
#'   "`average`" (= UPGMA), "`mcquitty`" (= WPGMA), "`median`" (= WPGMC) or
#'   "`centroid`" (= UPGMC). Due to size limitations with [stats::hclust()] --
#'   specifically the "`size cannot be NA nor exceed 65536`" --
#'   [fastcluster::hclust()] is being used because it is a complete replacement
#'   of [stats::hclust()], is fast (compared to [stats::hclust()]), and is able
#'   to accommodate dissimilarity matrices with more than 2^16 (65,536)
#'   observations.
#'
#' @return
#'   This function returns:
#'   - **h2o.clusters.raw**: Initial waters with assigned cluster ID
#'   - **h2o.clusters.summary**: Each cluster's:
#'       + cluster ID
#'       + number of waters
#'       + percent conservation
#'       + X, Y, and Z cooridinates
#'       + bound water environment measurements
#'       + mean distance between waters comprising the cluster
#'       + mean distance between waters comprising the cluster and the cluster's
#'         centroid
#'   - **h2o.occurrence**: A table indicating the structures (PDBs) contributing
#'       to each cluster. This summary table includes the PDB structure's:
#'       + resolution
#'       + R-free value
#'       + occupancy (mean and standard deviation)
#'       + mobility (mean and standard deviation)
#'       + B-value (mean and standard deviation)
#'       + number of waters in each cluster
#'       + number of waters passing the mobility cutoff
#'       + number of waters passing the normalized B-value
#'       + number of waters passing both cutoff values
#'       + percentage of waters passing both cutoffs
#'       + number of clusters the structure contributes to
#'       + True/False table indicating if the protein structure contributed to
#'         the water cluster
#'   - **clustering.info**: size and timing information
#'
#' @export
#'
#' @import fastcluster
#' @importFrom stats aggregate cutree dist sd
#' @importFrom utils object.size
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @references
#'   Paul C Sanschagrin and Leslie A Kuhn. Cluster analysis of
#'   consensus water sites in thrombin and trypsin shows conservation between
#'   serine proteases and contributions to ligand specificity. _Protein
#'   Science_, 1998, **7** (_10_), pp 2054-2064. \cr
#'   [DOI: 10.1002/pro.5560071002](http://doi.org/10.1002/pro.5560071002) \cr
#'   [PMID: 9792092](http://www.ncbi.nlm.nih.gov/pubmed/9792092) \cr
#'   [WatCH webpage](http://www.kuhnlab.bmb.msu.edu/software/watch/index.html) \cr
#'
#'   Fionn Murtagh and Pierre Legendre. Ward's Hierarchical Agglomerative
#'   Clustering Method: Which Algorithms Implement Ward's Criterion?
#'   _Journal of Classification_, 2014, **31**, (_3_), pp
#'   274-295. \cr
#'   [DOI: 10.1007/s00357-014-9161-z](http://doi.org/10.1007/s00357-014-9161-z)
#'
#'   Daniel Müllner. fastcluster: Fast Hierarchical, Agglomerative Clustering
#'   Routines for R and Python. _Journal of Statistical Software_, 2013, **53**
#'   (_9_) \cr
#'   [DOI: 10.18637/jss.v053.i09](http://dx.doi.org/10.18637/jss.v053.i09)
#'   [fastcluster webpage](http://danifold.net/fastcluster.html)
#'
ClusterWaters <- function(data, cutoff.cluster, cluster.method="complete") {

  ##_ basic information -----
  pdb.ids <- unique(data$PDBid)
  num.structures <- length(pdb.ids)
  colNames.data <- colnames(data)


  ##_ calculate the distance data -----
  ##--- message about number of pairwise distances being calculated
  num.h2o <- nrow(data)
  num.pairwise.distances.raw <- choose(num.h2o, 2)
  num.pairwise.distances <- prettyNum(num.pairwise.distances.raw,
                                      big.mark = ",")
  message.num.pairs <- paste("Calculating", num.pairwise.distances,
                             "pairwise distances for",
                             prettyNum(num.h2o, big.mark = ","),
                             "water molecules.", sep = " ")
  message(message.num.pairs)
  ##--- calulate the distances for ALL imported waters
  time.start <- Sys.time()
  h2o.distances <- dist(data[, c("x", "y", "z")],
                        method = "euclidean", diag = FALSE, upper = FALSE)
  time.pairwise.distances <- prettyNum(difftime(Sys.time(),
                                                time.start,
                                                units = "secs"))
  h2o.distances.size <- object.size(h2o.distances)
  h2o.distances.size.frmt <- format(h2o.distances.size,
                                    units = "auto")
  message.dist.time <- paste("The pairwise distance calculations for",
                             num.pairwise.distances, "water molecules took",
                             time.pairwise.distances, "seconds and is",
                             h2o.distances.size.frmt, "in size.", sep = " ")
  message(message.dist.time)


  ##_ cluster the waters -----
  message("Clustering the individual waters...")
  # h2o.hclust <- hclust(h2o.distances, method = cluster.method)
  h2o.hclust <- fastcluster::hclust(d = h2o.distances,
                                    method = cluster.method,
                                    members = NULL)

  ##_ determine clusters based on cluster.cutoff -----
  message("Constructing clusters...")
  h2o.clusters <- cutree(h2o.hclust, h = cutoff.cluster)
  cluster.idc <- unique(h2o.clusters)
  num.clusters <- max(h2o.clusters)
  ##--- add the cluster indices and percent conservation
  data <- cbind(data,
                cluster=h2o.clusters,
                stringsAsFactors = FALSE)
  ##--- reorder the waters based on cluster number
  data <- data[order(data$cluster, decreasing = FALSE), ]


  ##_ percent conservation -----
  h2o.per.cluster <- as.numeric(table(h2o.clusters))
  pct.conserved <- h2o.per.cluster / num.structures * 100
  pct.conserved <- rep(pct.conserved, h2o.per.cluster)

  ##--- add the cluster indices and percent conservation
  data <- cbind(data,
                pct.conserved = pct.conserved,
                stringsAsFactors = FALSE)

  ##--- reorder the waters based on percent conservation and cluster number
  data <- data[order(data$pct.conserved,
                     partial = data$cluster,
                     decreasing = TRUE), ]


  ##_ renumber the clusters so the lower the cluster number is correlated with -----
  ##  the more conserved clusters
  clusters.orig.idc <- data$cluster
  clusters.orig.unique <- unique(clusters.orig.idc)
  clusters.new.idc <- rep(NA, length(clusters.orig.idc))
  for ( cluster.idx in 1:num.clusters) {
    match.tf <- clusters.orig.idc %in% clusters.orig.unique[cluster.idx]

    clusters.new.idc[match.tf] <- cluster.idx
  }
  data$cluster <- clusters.new.idc


  ##_ table of structures with cluster present -----
  message("Constructing table of structures with cluster present...")
  cluster.occurrence <- as.data.frame(matrix(data=NA, nrow=num.structures, ncol=num.clusters))
  rownames(cluster.occurrence) <- pdb.ids
  colnames(cluster.occurrence) <- paste("clust", formatC(cluster.idc, flag="0", format="d", width=nchar(num.clusters)), sep="_")
  ##--- determine which structures have a water that is part of a cluster
  for (cluster.idx in cluster.idc) {
    cluster.strks <- data$PDBid[data$cluster == cluster.idx]
    cluster.occurrence[, cluster.idx] <- pdb.ids %in% cluster.strks
  }
  ##--- add experimental data about structures and passing of cutoffs
  list.PDBid <- list(data$PDBid)
  resolution <- as.numeric(as.vector(aggregate(data[, "resolution"], list.PDBid, unique)[, 2]))
  rObserved <- as.numeric(as.vector(aggregate(data[, "rObserved"], list.PDBid, unique)[, 2]))
  rFree <- as.numeric(as.vector(aggregate(data[, "rFree"], list.PDBid, unique)[, 2]))
  occupancy.mu <- aggregate(data[, "o"], list.PDBid, mean)[, 2]
  occupancy.sd <- aggregate(data[, "o"], list.PDBid, sd)[, 2]
  mobility.mu <- aggregate(data[, "mobility"], list.PDBid, mean)[, 2]
  mobility.sd <- aggregate(data[, "mobility"], list.PDBid, sd)[, 2]
  Bvalue.mu <- aggregate(data[, "b"], list.PDBid, mean)[, 2]
  Bvalue.sd <- aggregate(data[, "b"], list.PDBid, sd)[, 2]
  nBvalue.mu <- aggregate(data[, "nBvalue"], list.PDBid, mean, na.rm=TRUE)[, 2]
  nBvalue.sd <- aggregate(data[, "nBvalue"], list.PDBid, sd, na.rm=TRUE)[, 2]
  num.h2o <- aggregate(data[, "PDBid"], list.PDBid, length)[, 2]
  num.pass.mobility <- aggregate(data[, "mobility.keep"], list.PDBid, sum)[, 2]
  num.pass.nBvalue <- aggregate(data[, "nBvalue.keep"], list.PDBid, sum)[, 2]
  num.pass.h2o.evaluations <- aggregate(data[, "passed.cutoffs"], list.PDBid, sum)[, 2]
  pct.pass.cutoffs <- (num.pass.h2o.evaluations / num.h2o) * 100
  num.clusters.per.strk <- rowSums(cluster.occurrence)
  ##--- create the cluster occurrence results data.frame
  cluster.occurrence <- cbind(resolution=resolution,
                              rObserved=rObserved,
                              rFree=rFree,
                              occupancy.mu=occupancy.mu,
                              occupancy.sd=occupancy.sd,
                              mobility.mu=mobility.mu,
                              mobility.sd=mobility.sd,
                              Bvalue.mu=Bvalue.mu,
                              Bvalue.sd=Bvalue.sd,
                              nBvalue.mu=nBvalue.mu,
                              nBvalue.sd=nBvalue.sd,
                              num.waters=num.h2o,
                              num.pass.mobility=num.pass.mobility,
                              num.pass.nBvalue=num.pass.nBvalue,
                              num.pass.cutoffs=num.pass.h2o.evaluations,
                              pct.pass.cutoffs=pct.pass.cutoffs,
                              num.clusters=num.clusters.per.strk,
                              cluster.occurrence,
                              stringsAsFactors=FALSE)


  ##_ construct the centroid for each cluster along with calculate the -----
  ##  various parameters for each cluster
  ##- columns of interest for the mean and sd values
  message("Constructing data.frame with cluster information...")
  mean.coi <- c("x", "y", "z",
                "o", "b", "mobility", "nBvalue", "pct.conserved")
  sd.coi <- c("o", "b", "mobility", "nBvalue")
  ##--- mean and standard deviation
  cluster.mu <- aggregate(data[, mean.coi], list(data$cluster), mean)
  cluster.sd <- aggregate(data[, sd.coi], list(data$cluster), sd)
  ##--- number of waters in each cluster
  h2o.per.cluster <- as.numeric(table(data$cluster))

  ##_ calculate the RMSF, B-value, and mobility for each cluster -----
  ##--- calculate rmsf for each cluster
  num.conserved.clusters <- max(data$cluster)
  cluster.rmsfValues <- rep(NA, num.clusters)
  for ( cluster in 1:num.conserved.clusters ) {
    cluster.rmsfValues[cluster] <- bio3d::rmsf(data[data$cluster == cluster, c("x","y","z")])
  }

  ##--- calculate B-value for each cluster
  b.calc <- calcBvalue(cluster.rmsfValues)
  b.calc.gt100.tf <- b.calc > 100.00000
  b.calc[b.calc.gt100.tf] <- 100

  ##--- calculate Mobility values for each cluster
  cluster.oValues <- cluster.mu$o
  mobility.calc <- Mobility(b.calc, cluster.oValues)

  ##--- update the B-values and Mobility values
  cluster.mu <- cbind(cluster.mu, b.calc, mobility.calc)

  ##_ calculate and add the bound water environment data -----
  list.cluster <- list(data$cluster)
  adn.sums <- aggregate(data[, grep(pattern = ".adn", colNames.data)], list.cluster, sum)
  hbonds.sums <- aggregate(data[, grep(pattern = ".hbonds", colNames.data)], list.cluster, sum)
  cluster.sums <- aggregate(data[, grep(pattern = ".sum", colNames.data)], list.cluster, sum)
  cluster.prot.mu <- cluster.sums[, grep(pattern = "prot.", colnames(cluster.sums))] / adn.sums[, "prot.adn"]
  cluster.h2o.mu <- cluster.sums[, grep(pattern = "h2o.", colnames(cluster.sums))] / adn.sums[, "h2o.adn"]
  cluster.het.mu <- cluster.sums[, grep(pattern = "het.", colnames(cluster.sums))] / adn.sums[, "het.adn"]
  cluster.all.mu <- cluster.sums[, grep(pattern = "all.", colnames(cluster.sums))] / adn.sums[, "all.adn"]
  ##-- combine the results
  df.bwe.cluster <- cbind(prot.adn = adn.sums[, "prot.adn"], prot.hbonds = hbonds.sums[, "prot.hbonds"], cluster.prot.mu,
                          h2o.adn = adn.sums[, "h2o.adn"], h2o.hbonds = hbonds.sums[, "h2o.hbonds"], cluster.h2o.mu,
                          het.adn = adn.sums[, "het.adn"], het.hbonds = hbonds.sums[, "het.hbonds"], cluster.het.mu,
                          all.adn = adn.sums[, "all.adn"], all.hbonds = hbonds.sums[, "all.hbonds"], cluster.all.mu)
  colnames(df.bwe.cluster) <- gsub(pattern = ".sum", replacement = ".mu",
                                   x = colnames(df.bwe.cluster))

  ##----- mean distance (and standard deviation) between waters in the clusters
  ##      AND
  ##----- mean distance (and standard deviation) between cluster centroid and
  ##      waters in the cluster
  time.start <- Sys.time()
  message("Calculating the mean distance between waters AND between the centroid within each cluster...\n")
  dist.mu <- dist.sd <- dist.centroid.mu <- dist.centroid.sd <- rep(NA, num.clusters)
  for (h2o.cluster in cluster.idc) {

    ##--- when there is only ONE water in the cluster
    if (h2o.per.cluster[h2o.cluster] == 1) {
      dist.mu[h2o.cluster] <- 0
      dist.sd[h2o.cluster] <- 0
      dist.centroid.mu[h2o.cluster] <- 0
      dist.centroid.sd[h2o.cluster] <- 0
    } else {
    ##--- when there are multiple waters in the cluster
      ##--- the distance between waters within the cluster
      h2o.xyz <- data[data$cluster == h2o.cluster, c("x", "y", "z")]
      h2o.cluster.dist <- as.vector(dist(h2o.xyz,
                                         method = "euclidean",
                                         diag = FALSE, upper = FALSE))
      ##- insert the values
      dist.mu[h2o.cluster] <- mean(h2o.cluster.dist, na.rm=TRUE)
      dist.sd[h2o.cluster] <- sd(h2o.cluster.dist, na.rm=TRUE)

      ##--- the distance between waters within the cluster AND the centroid
      cluster.xyz <- rbind(cluster.mu[h2o.cluster, c("x", "y", "z")],
                           h2o.xyz)
      cluster.centroid.dist <- as.matrix(dist(cluster.xyz,
                                              method = "euclidean",
                                              diag = TRUE, upper = TRUE))
      cluster.centroid.dists <- as.vector(cluster.centroid.dist[-1, 1])
      ##- insert the values
      dist.centroid.mu[h2o.cluster] <- mean(cluster.centroid.dists)
      dist.centroid.sd[h2o.cluster] <- sd(cluster.centroid.dists)
    }

  }
  time.cluster.distances <- prettyNum(difftime(Sys.time(),
                                               time.start, units="secs"))

  ##_ the summary -----
  data.summary <- data.frame(cluster=cluster.idc,
                             num.waters=h2o.per.cluster,
                             pct.conserved=cluster.mu$pct.conserved,
                             x=cluster.mu$x,
                             y=cluster.mu$y,
                             z=cluster.mu$z,
                             o.mu=cluster.mu$o, o.sd=cluster.sd$o,
                             b.exp.mu=cluster.mu$b, b.exp.sd=cluster.sd$b,
                             b.calc=cluster.mu$b.calc,
                             mobility.mu=cluster.mu$mobility, mobility.sd=cluster.sd$mobility,
                             mobility.calc=cluster.mu$mobility.calc,
                             nBvalue.mu=cluster.mu$nBvalue, nBvalue.sd=cluster.sd$nBvalue,
                             df.bwe.cluster,
                             dist.mu=dist.mu, dist.sd=dist.sd,
                             dist.centroid.mu=dist.centroid.mu, dist.centroid.sd=dist.centroid.sd,
                             stringsAsFactors=FALSE)
  ##--- convert NaNs to NAs
  data.summary[is.na(data.summary)] <- NA


  ##_ timings and size -----
  clustering.info <- list(num.pairwise.dist=num.pairwise.distances.raw,
                          size.pairwise.dist=h2o.distances.size.frmt,
                          time.pairwise.dist=time.pairwise.distances,
                          time.cluster.dist=time.cluster.distances)


  ##_ return the data -----
  list(h2o.clusters.raw=data,
       h2o.clusters.summary=data.summary,
       h2o.occurrence=cluster.occurrence,
       clustering.info=clustering.info)

}





## cluster MDS waters docs -----------------------------------------------------
#' @title Cluster Conserved Waters (MDS)
#' @description Cluster the conserved waters from a molecular dynamics
#'   simulation trajectory.
#' @details Calculate the conserved waters using a molecular dynamics simulation
#'   trajectory.
#'
#' @param data The water oxygens' X, Y, and Z coordinates.
#' @param cutoff.cluster Numerical value provided by the user for the distance
#'   between water oxygen atoms to form a cluster; default: 2.4 Angstroms.
#' @param cluster.method Method of clustering the waters; default is
#'   "`complete`". Any other method accepted by the `hclust` function is
#'   appropriate. The original method used by Sanschagrin and Kuhn is the
#'   complete linkage clustering method and is the default. Other options
#'   include "`ward.D`" (equivilant to the only Ward option in `R` versions
#'   3.0.3 and earlier), "`ward.D2`" (implements Ward's 1963 criteria; see
#'   Murtagh and Legendre 2014), "`single`" (related to the minimal spanning
#'   tree method and adopts a "friend of friends" clustering method), along with
#'   "`average`" (= UPGMA), "`mcquitty`" (= WPGMA), "`median`" (= WPGMC) or
#'   "`centroid`" (= UPGMC). Due to size limitations with [stats::hclust()] --
#'   specifically the "`size cannot be NA nor exceed 65536`" --
#'   [fastcluster::hclust()] is being used because it is a complete replacement
#'   of [stats::hclust()], is fast (compared to [stats::hclust()]), and is able
#'   to accommodate dissimilarity matrices with more than 2^16 (65,536)
#'   observations.
#'
#' @return
#'   This function returns:
#'   - **h2o.clusters.raw**: Initial waters with assigned cluster ID
#'   - **h2o.clusters.summary**: Each cluster's:
#'       + cluster ID
#'       + number of waters
#'       + percent conservation
#'       + X, Y, and Z cooridinates
#'       + bound water environment measurements
#'       + mean distance between waters comprising the cluster
#'       + mean distance between waters comprising the cluster and the cluster's
#'         centroid
#'   - **h2o.occurrence**: A table indicating the structures (PDBs) contributing
#'       to each cluster. This summary table includes the PDB structure's:
#'       + number of waters in each cluster
#'       + number of clusters the structure contributes to
#'       + True/False table indicating if the protein structure contributed to
#'         the water cluster
#'   - **clustering.info**: size and timing information
#'
#' @export
#'
#' @import bio3d
#' @import fastcluster
#' @importFrom stats aggregate cutree dist sd
#' @importFrom utils object.size
#'
#' @author Emilio Xavier Esposito \email{emilio@@exeResearch.com}
#'
#' @references
#'   Paul C Sanschagrin and Leslie A Kuhn. Cluster analysis of
#'   consensus water sites in thrombin and trypsin shows conservation between
#'   serine proteases and contributions to ligand specificity. _Protein
#'   Science_, 1998, **7** (_10_), pp 2054-2064. \cr
#'   [DOI: 10.1002/pro.5560071002](http://doi.org/10.1002/pro.5560071002) \cr
#'   [PMID: 9792092](http://www.ncbi.nlm.nih.gov/pubmed/9792092) \cr
#'   [WatCH webpage](http://www.kuhnlab.bmb.msu.edu/software/watch/index.html) \cr
#'
#'   Fionn Murtagh and Pierre Legendre. Ward's Hierarchical Agglomerative
#'   Clustering Method: Which Algorithms Implement Ward's Criterion?
#'   _Journal of Classification_, 2014, **31**, (_3_), pp
#'   274-295. \cr
#'   [DOI: 10.1007/s00357-014-9161-z](http://doi.org/10.1007/s00357-014-9161-z)
#'
#'   Daniel Müllner. fastcluster: Fast Hierarchical, Agglomerative Clustering
#'   Routines for R and Python. _Journal of Statistical Software_, 2013, **53**
#'   (_9_) \cr
#'   [DOI: 10.18637/jss.v053.i09](http://dx.doi.org/10.18637/jss.v053.i09)
#'   [fastcluster webpage](http://danifold.net/fastcluster.html)
#'
ClusterWaters.MDS <- function(data,
                              cutoff.cluster,
                              cluster.method = "complete") {

  ##----- basic information
  pdb.ids <- unique(data$PDBid)
  num.structures <- length(pdb.ids)
  colNames.data <- colnames(data)


  ##----- calculate the distance data
  ##--- message about number of pairwise distances being calculated
  num.h2o <- nrow(data)
  num.pairwise.distances.raw <- choose(num.h2o, 2)
  num.pairwise.distances <- prettyNum(num.pairwise.distances.raw,
                                      big.mark = ",")
  message.num.pairs <- paste("Calculating", num.pairwise.distances,
                             "pairwise distances for",
                             prettyNum(num.h2o, big.mark = ","),
                             "water molecules.", sep = " ")
  message(message.num.pairs)
  ##--- calulate the distances for ALL imported waters
  time.start <- Sys.time()
  h2o.distances <- dist(data[, c("x", "y", "z")],
                        method = "euclidean", diag = FALSE, upper = FALSE)
  time.pairwise.distances <- prettyNum(difftime(Sys.time(),
                                                time.start,
                                                units = "secs"))
  h2o.distances.size <- object.size(h2o.distances)
  h2o.distances.size.frmt <- format(h2o.distances.size,
                                    units = "auto")
  message.dist.time <- paste("The pairwise distance calculations for",
                             num.pairwise.distances, "water molecules took",
                             time.pairwise.distances, "seconds and is",
                             h2o.distances.size.frmt, "in size.", sep = " ")
  message(message.dist.time)


  ##----- cluster the waters
  message("Clustering the individual waters...")
  # h2o.hclust <- hclust(h2o.distances, method = cluster.method)
  h2o.hclust <- fastcluster::hclust(d = h2o.distances,
                                    method = cluster.method,
                                    members = NULL)


  ##----- determine clusters based on cluster.cutoff
  message("Constructing clusters...")
  h2o.clusters <- cutree(h2o.hclust, h = cutoff.cluster)
  cluster.idc <- unique(h2o.clusters)
  num.clusters <- max(h2o.clusters)
  ##--- add the cluster indices and percent conservation
  data <- cbind(data,
                cluster=h2o.clusters,
                stringsAsFactors = FALSE)
  ##--- reorder the waters based on cluster number
  data <- data[order(data$cluster, decreasing = FALSE), ]


  ##----- percent conservation
  h2o.per.cluster <- as.numeric(table(h2o.clusters))
  pct.conserved <- h2o.per.cluster / num.structures * 100
  pct.conserved <- rep(pct.conserved, h2o.per.cluster)

  ##--- add the cluster indices and percent conservation
  data <- cbind(data,
                pct.conserved = pct.conserved,
                stringsAsFactors = FALSE)

  ##--- reorder the waters based on percent conservation and cluster number
  data <- data[order(data$pct.conserved,
                     partial = data$cluster,
                     decreasing = TRUE), ]


  ##----- renumber the clusters so the lower the cluster number is correlated with
  ##      the more conserved clusters
  clusters.orig.idc <- data$cluster
  clusters.orig.unique <- unique(clusters.orig.idc)
  clusters.new.idc <- rep(NA, length(clusters.orig.idc))
  for ( cluster.idx in 1:num.clusters) {
    match.tf <- clusters.orig.idc %in% clusters.orig.unique[cluster.idx]

    clusters.new.idc[match.tf] <- cluster.idx
  }
  data$cluster <- clusters.new.idc


  ##----- table of structures with cluster present
  message("Constructing table of structures with cluster present...")
  cluster.occurrence <- as.data.frame(matrix(data = NA,
                                             nrow = num.structures,
                                             ncol = num.clusters))
  rownames(cluster.occurrence) <- pdb.ids
  colnames(cluster.occurrence) <- paste("clust",
                                        formatC(cluster.idc,
                                                flag = "0",
                                                format = "d",
                                                width = nchar(num.clusters)),
                                        sep = "_")
  ##--- determine which structures have a water that is part of a cluster
  for (cluster.idx in cluster.idc) {
    cluster.strks <- data$PDBid[data$cluster == cluster.idx]
    cluster.occurrence[, cluster.idx] <- pdb.ids %in% cluster.strks
  }
  ##--- add experimental data about structures and passing of cutoffs
  list.PDBid <- list(data$PDBid)
  num.h2o <- aggregate(data[, "PDBid"], list.PDBid, length)[, 2]
  num.clusters.per.strk <- rowSums(cluster.occurrence)
  ##--- create the cluster occurrence results data.frame
  cluster.occurrence <- cbind(num.waters = num.h2o,
                              num.clusters = num.clusters.per.strk,
                              cluster.occurrence)


  ##----- construct the centroid for each cluster along with calculate the
  ##      various parameters for each cluster
  ##- columns of interest for the mean and sd values
  message("Constructing data.frame with cluster information...")
  mean.coi <- c("x", "y", "z",
                "pct.conserved")
  ##--- mean and standard deviation
  cluster.mu <- aggregate(data[, mean.coi], list(data$cluster), mean)
  ##--- number of waters in each cluster
  h2o.per.cluster <- as.numeric(table(data$cluster))

  ##_ calculate the RMSF, B-value, and mobility for each cluster -----
  ##--- calculate rmsf for each cluster
  num.conserved.clusters <- max(data$cluster)
  cluster.rmsfValues <- rep(NA, num.clusters)
  for ( cluster in 1:num.conserved.clusters ) {
    cluster.rmsfValues[cluster] <- bio3d::rmsf(data[data$cluster == cluster, c("x","y","z")])
  }

  ##--- calculate B-value for each cluster
  cluster.Bvalues <- calcBvalue(cluster.rmsfValues)
  cluster.Bvalues.gt100.tf <- cluster.Bvalues > 100.00000
  cluster.Bvalues[cluster.Bvalues.gt100.tf] <- 100

  ##--- calculate Mobility values for each cluster
  # cluster.oValues <- rep(1, num.conserved.clusters)
  cluster.oValues <- cluster.mu$pct.conserved
  cluster.MobilityValues <- Mobility(cluster.Bvalues, cluster.oValues)


  ##--- add the bound water environment data
  list.cluster <- list(data$cluster)
  adn.sums <- aggregate(data[, grep(pattern = ".adn", colNames.data)], list.cluster, sum)
  hbonds.sums <- aggregate(data[, grep(pattern = ".hbonds", colNames.data)], list.cluster, sum)
  cluster.sums <- aggregate(data[, grep(pattern = ".sum", colNames.data)], list.cluster, sum)
  cluster.prot.mu <- cluster.sums[, grep(pattern = "prot.", colnames(cluster.sums))] / adn.sums[, "prot.adn"]
  cluster.h2o.mu <- cluster.sums[, grep(pattern = "h2o.", colnames(cluster.sums))] / adn.sums[, "h2o.adn"]
  cluster.het.mu <- cluster.sums[, grep(pattern = "het.", colnames(cluster.sums))] / adn.sums[, "het.adn"]
  cluster.all.mu <- cluster.sums[, grep(pattern = "all.", colnames(cluster.sums))] / adn.sums[, "all.adn"]
  ##-- combine the results
  df.bwe.cluster <- cbind(prot.adn = adn.sums[, "prot.adn"],
                          prot.hbonds = hbonds.sums[, "prot.hbonds"],
                          cluster.prot.mu,
                          h2o.adn = adn.sums[, "h2o.adn"],
                          h2o.hbonds = hbonds.sums[, "h2o.hbonds"],
                          cluster.h2o.mu,
                          het.adn = adn.sums[, "het.adn"],
                          het.hbonds = hbonds.sums[, "het.hbonds"],
                          cluster.het.mu,
                          all.adn = adn.sums[, "all.adn"],
                          all.hbonds = hbonds.sums[, "all.hbonds"],
                          cluster.all.mu)
  colnames(df.bwe.cluster) <- gsub(pattern = ".sum", replacement = ".mu",
                                   x = colnames(df.bwe.cluster))

  ##----- mean distance (and standard deviation) between waters in the clusters
  ##      AND
  ##----- mean distance (and standard deviation) between cluster centroid and
  ##      waters in the cluster
  time.start <- Sys.time()
  message("Calculating the mean distance between waters AND between the centroid within each cluster...\n")
  dist.mu <- dist.sd <- dist.centroid.mu <- dist.centroid.sd <- rep(NA, num.clusters)
  for (h2o.cluster in cluster.idc) {

    ##--- when there is only ONE water in the cluster
    if (h2o.per.cluster[h2o.cluster] == 1) {
      dist.mu[h2o.cluster] <- 0
      dist.sd[h2o.cluster] <- 0
      dist.centroid.mu[h2o.cluster] <- 0
      dist.centroid.sd[h2o.cluster] <- 0
    } else {
      ##--- when there are multiple waters in the cluster
      ##--- the distance between waters within the cluster
      h2o.xyz <- data[data$cluster == h2o.cluster, c("x", "y", "z")]
      h2o.cluster.dist <- as.vector(dist(h2o.xyz,
                                         method = "euclidean",
                                         diag = FALSE, upper = FALSE))
      ##- insert the values
      dist.mu[h2o.cluster] <- mean(h2o.cluster.dist, na.rm=TRUE)
      dist.sd[h2o.cluster] <- sd(h2o.cluster.dist, na.rm=TRUE)

      ##--- the distance between waters within the cluster AND the centroid
      cluster.xyz <- rbind(cluster.mu[h2o.cluster, c("x", "y", "z")],
                           h2o.xyz)
      cluster.centroid.dist <- as.matrix(dist(cluster.xyz,
                                              method = "euclidean",
                                              diag = TRUE, upper = TRUE))
      cluster.centroid.dists <- as.vector(cluster.centroid.dist[-1, 1])
      ##- insert the values
      dist.centroid.mu[h2o.cluster] <- mean(cluster.centroid.dists)
      dist.centroid.sd[h2o.cluster] <- sd(cluster.centroid.dists)
    }

  }
  time.cluster.distances <- prettyNum(difftime(Sys.time(),
                                               time.start, units = "secs"))

  ##--- the summary
  data.summary <- data.frame(cluster = cluster.idc,
                             num.waters = h2o.per.cluster,
                             pct.conserved = cluster.mu$pct.conserved,
                             x = cluster.mu$x,
                             y = cluster.mu$y,
                             z = cluster.mu$z,
                             o = cluster.oValues,
                             b.calc = cluster.Bvalues,
                             rmsf = cluster.rmsfValues,
                             mobility = cluster.MobilityValues,
                             df.bwe.cluster,
                             dist.mu = dist.mu,
                             dist.sd = dist.sd,
                             dist.centroid.mu = dist.centroid.mu,
                             dist.centroid.sd = dist.centroid.sd,
                             stringsAsFactors = FALSE)
  ##--- convert NaNs to NAs
  data.summary[is.na(data.summary)] <- NA


  ##----- timings and size
  clustering.info <- list(num.pairwise.dist = num.pairwise.distances.raw,
                          size.pairwise.dist = h2o.distances.size.frmt,
                          time.pairwise.dist = time.pairwise.distances,
                          time.cluster.dist = time.cluster.distances)


  ##----- return the data
  list(h2o.clusters.raw = data,
       h2o.clusters.summary = data.summary,
       h2o.occurrence = cluster.occurrence,
       clustering.info = clustering.info)

}


