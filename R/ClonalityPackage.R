#' Cluster mutations within ditanceTolerance
#'
#' @param eventDataFrame Table of mutations in the dataframe format
#' @param identityColumns The columns that identify a event. All events with the same values in the
#'   identyColumns will be grouped together.
#' @param ditanceTolerance Distance tolerance
#' @export
clusterCloseMutations <- function(eventDataFrame, identityColumns, distanceTolerance)
{
  events <- unique(sapply(1:nrow(eventDataFrame),function(i){paste(eventDataFrame[i,identityColumns[identityColumns!="pos"]], collapse="_")}))
  coarseGrainedTable <- foreach (event = events, .combine = rbind, .init=data.frame(),.errorhandling = c("stop", "remove", "pass")[1], .inorder=FALSE) %dopar%
  {
    eventDataFrame_sub <- eventDataFrame[which(events==event),]
    if (nrow(eventDataFrame_sub)<2)
    {
      eventDataFrame_sub$cluster <- 1
      return(eventDataFrame_sub)
    }
    else
    {
      D <- dist(eventDataFrame_sub[,"pos"], method = "manhattan")
      clusters <- hclust(D, method="complete")
      eventDataFrame_sub$cluster <- cutree(clusters, h = distanceTolerance) 
      return(eventDataFrame_sub)
    }
  }
  return(coarseGrainedTable)
}

#' Adds mutation IDs to the table of mutations using only identityColumns
#'
#' @param eventDataFrame Table of mutations in the dataframe format
#' @param identityColumns The columns that identify a event. All events with the same values in the
#'   identyColumns will be grouped together.
#' @export
addMutIDtoTableAndGetMappingMutIDtoMut <- function(eventDataFrame, identityColumns, samples)
{
  mapMutID <- unique(eventDataFrame[,identityColumns])
  rownames(mapMutID) <- 1:nrow(mapMutID)
  mapMutID$mutID <- 1:nrow(mapMutID)
  eventDataFrame <- merge(eventDataFrame,mapMutID, by=identityColumns)
  eventDataFrame$freq <- sapply(eventDataFrame$mutID, function(i) {sum(eventDataFrame$mutID == i)})/length(samples)
  return(eventDataFrame)
}

#' Converts a dataframe of events to a matrix of events by samples.
#'
#' @param eventDataFrame Table of mutations in the dataframe format
#' @param identityColumns The columns that identify a event. All events with the same values in the
#'   identyColumns will be grouped together.
#' @param samples The column that contains the sample identifiers.
#'
#' @export
eventDataFrameToMatrix <- function(pairs, eventDataFrame, identityColumns, samples=NULL)
{
  if (is.null(samples))
  {
    samples <- unique(c(pairs[, 1], pair[, 2]))
  }
  eventMatrix <- as.matrix(matrix(0, nrow = length(samples), ncol=length(unique(eventDataFrame$mutID))))
  rownames(eventMatrix) <- samples
  colnames(eventMatrix) <- unique(eventDataFrame$mutID)
  
  for (i in 1:nrow(eventDataFrame)) 
  {
    eventMatrix[eventDataFrame$SampleID[i], eventDataFrame$mutID[i]] <- 1
  }
  return(eventMatrix)
}

makeEventExpectation <- function(eventMatrix)
{
  library(discover)
  return(t(discover.matrix(t(eventMatrix))$bg))
}

#' Score the clonality of a single pair of samples.
#' 
#' @param events1 A binary vector with the events of the first sample.
#' @param events2 A binary vector with the events of the second sample.
#' @param eventExpectation1 A binary vector with the expectation of all events of the first sample, row of the eventExpectation for the first sample.
#' @param eventExpectation2 A binary vector with the expectation of all events of the second sample, row of the eventExpectation for the second sampke.
#'
#' @export
scoreClonality <- function(events1, events2, eventExpectation1, eventExpectation2)
{
  mutNum1 <- length(events1)
  mutNum2 <- length(events2)
  indSharedMut <- intersect(events1,events2)
  sharedMutNum <- length(indSharedMut)
  p1 <- eventExpectation1[indSharedMut]
  p2 <- eventExpectation2[indSharedMut]
  score <- sum( - log(p1*p2/(1-(1-p1)*(1-p2))) )
  return(data.frame(score = score, mutNum1 = mutNum1, mutNum2 = mutNum2, sharedMutNum = sharedMutNum))
}

#'  Compute clonality scores for the pairs
#' 
#' @param pairs A dataframe of sample pairs with sample identifiers in columns sample1 and sample2.
#' @param eventMatrix A binary matrix indicating events (columns) in samples (rows).
#' @param eventExpection A matrix indicating the expectation of events (columns) in samples (rows).
#' @param scorefun The function used to score a pair of samples.
#' @return pairs with an added column scores containing the clonality scores.
#'
#' @export
computeClonalityScores <- function(pairs, eventMatrix, eventExpectation, scoreFun=scoreClonality)
{
  clonalityScores <- foreach(i = 1:nrow(pairs), .combine = rbind, .inorder=TRUE, .init=data.frame()) %dopar%
  {
    cbind(pairs[i, ], scoreFun(which(eventMatrix[pairs[i, 1], ] == 1), which(eventMatrix[pairs[i, 2], ] == 1), eventExpectation[pairs[i, 1], ], eventExpectation[pairs[i, 2], ]) )
  }
  return(clonalityScores)
}    


#'  Compute table of shared mutations
#' 
#' @param pairs A dataframe of sample pairs with sample identifiers in columns sample1 and sample2.
#' @param eventMatrix A binary matrix indicating events (columns) in samples (rows).
#' @param eventExpection A matrix indicating the expectation of events (columns) in samples (rows).
#' @param scorefun The function used to score a pair of samples.
#' @return pairs with an added column scores containing the clonality scores.
#'
#' @export
computeSharedMutations <- function(pairs, eventDataFrame, identityColumns)
{
  sharedMutations  <- foreach(i = 1:nrow(pairs), .combine = rbind, .inorder=TRUE, .init=data.frame()) %dopar%
  {
    data <- merge(eventDataFrame[eventDataFrame$SampleID== pairs[i, 1], ], eventDataFrame[eventDataFrame$SampleID== pairs[i, 2], ], by = c(identityColumns, "mutID", "freq"))
  }
  colnames(sharedMutations)[colnames(sharedMutations)=="SampleID.x"] <- "sample1"
  colnames(sharedMutations)[colnames(sharedMutations)=="SampleID.y"] <- "sample2"
  return(sharedMutations)
}    


#' Generate a reference distribution by computing the score for all random pairs from different
#' patients and calculate p-value.
#' 
#' @param pairs A dataframe of sample pairs with sample identifiers in columns sample1 and sample2.
#'   Column Patient is optional. If it is not available it is assumed all pairs are from different
#'   patients. Column score gives the score for the pairs.
#' @param events A binary matrix indicating events (columns) in samples (rows).
#' @param events A matrix indicating the expectation of events (columns) in samples (rows).
#' @param scorefun The function used to score a pair of samples.
#' @return pairs with an added column p containing p-values from the randomization test.
#'
#' @export
testClonalityScores  <- function(pairs, eventMatrix, eventExpectation, scoreFun=scoreClonality) 
{
  if (is.null(pairs$patient))
  {
    pairs$patient <- 1:nrow(pairs)
  }
  refPairs <- expand.grid(list(Sample1=unique(pairs[[1]]), Sample2=unique(pairs[[2]])), stringsAsFactors = FALSE)
  refPairs <- refPairs[which(patients[refPairs$Sample1,1] != patients[refPairs$Sample2,1]),]
  return(computeClonalityScores(refPairs, eventMatrix, eventExpectation, scoreFun=scoreClonality) )
}      


#' Compute clonality score and p-values of paired samples.
#'
#' @param eventDataFrame table of mutations in the dataframe format
#' @param pairs A dataframe of sample pairs with sample identifiers in columns sample1 and sample2.
#'   Column Patient is optional. If it is not available it is assumed all pairs are from different
#'   patients.
#' @param identityColumns relevant columns
#' @param samples samples to analyse
#' @param ditanceTolerance Distance tolerance to cluster mutations
#' @return pairs with extra columns scorei, number of shared breaks and p-value
#'
#' @export
clonality <- function(pairs, eventDataFrame, identityColumns, samples, distanceTolerance=NULL)
{
  print(paste0("Start analyzing the data: ", nrow(pairs), " pairs, ", nrow(eventDataFrame), " mutations, ", length(samples), " samples, identytiy columns are ", paste0(identityColumns, collapse=' ')))
  if (!is.null(distanceTolerance) )
  {
    if (distanceTolerance > 0)
    {
      eventDataFrame <- clusterCloseMutations(eventDataFrame, identityColumns, distanceTolerance)
      identityColumns[identityColumns=="pos"] <- "cluster"
      print("Clustered mutations and made cluster as identityColumns instead of pos")
    }
  }
  eventDataFrame <- addMutIDtoTableAndGetMappingMutIDtoMut(eventDataFrame, identityColumns, samples)
  print("Added mutations IDs to the mutation data.frame")
  eventMatrix <- eventDataFrameToMatrix(pairs, eventDataFrame, identityColumns, samples)
  print("Made mutation matrix")
  eventExpectation <- makeEventExpectation(eventMatrix)
  print("Made mutation expectation matrix (using DISCOVER)")
  referenceScores <- testClonalityScores(pairs, eventMatrix, eventExpectation, scoreFun=scoreClonality) 
  clonalityScores <- computeClonalityScores(pairs, eventMatrix, eventExpectation, scoreFun=scoreClonality)
  clonalityScores$p.value <- sapply(clonalityScores$score, function(s) {mean(s <= c(referenceScores$score,s))} )
  sharedMutations <- computeSharedMutations(pairs, eventDataFrame, identityColumns)
  sharedMutations <- merge(sharedMutations, clonalityScores[, c("sample1", "sample2", "p.value")], by=c("sample1", "sample2"))
  return(list(clonalityScores, sharedMutations))
}
