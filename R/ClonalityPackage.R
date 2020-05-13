clusterCloseMutations <- function(eventDataFrame, relevantColumns, gap)
{
  events <- unique(sapply(1:nrow(eventDataFrame),function(i){paste(eventDataFrame[i,relevantColumns[relevantColumns!="pos"]], collapse="_")}))
  CoarseGrainedTable <- foreach (event = events, .combine = rbind, .init=data.frame(),.errorhandling = c("stop", "remove", "pass")[1], .inorder=FALSE) %dopar%
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
      eventDataFrame_sub$cluster <- cutree(clusters, h = gap) 
      return(eventDataFrame_sub)
    }
  }
  return(CoarseGrainedTable)
}


addMutIDtoTableAndGetMappingMutIDtoMut <- function(eventDataFrame, identityColumns)
{
  mapMutID <- unique(eventDataFrame[,identityColumns])
  rownames(mapMutID) <- 1:nrow(mapMutID)
  eventDataFrame$mutID <- 0
  for (i in 1:nrow(mapMutID))
  {
    Ind <- sapply(1:nrow(eventDataFrame),function(j){all(mapMutID[i,] == eventDataFrame[j, identityColumns])})
    eventDataFrame$mutID[Ind] <- i
  }
  return(eventDataFrame)
}

#' Converts a dataframe of events to a matrix of events by samples.
#'
#' @param eventDataFrame Table of mutations in the dataframe format
#' @param identityColumns The columns that identify a event. All events with the same values in the
#'   identyColumns will be grouped together.
#' @param sampleColumn The column that contains the sample identifiers.
#'
#' @export
eventDataFrameToMatrix <- function(eventDataFrame, identityColumns, sampleColumn=NULL)
{
  eventDataFrame <- addMutIDtoTableAndGetMappingMutIDtoMut(eventDataFrame, identityColumns)
  if (is.null(sampleColumn))
  {
    sampleColumn <- unique(eventDataFrame$SampleID)
  }
  eventMatrix <- as.matrix(matrix(0, nrow = length(sampleColumn), ncol=length(unique(eventDataFrame$mutID))))
  rownames(eventMatrix) <- sampleColumn
  colnames(eventMatrix) <- unique(eventDataFrame$mutID)
  
  for (i in 1:nrow(eventDataFrame)) 
  {
    eventMatrix[eventDataFrame$SampleID[i],eventDataFrame$mutID[i]] <- 1
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
  clonalityScores <- foreach(i = 1:nrow(pairs), .combine = rbind, .inorder=TRUE, .init=data.frame()) %do%
  {
    output <- scoreFun(which(eventMatrix[pairs[i,1],]==1), which(eventMatrix[pairs[i,2],]==1), eventExpectation[pairs[i,1]], eventExpectation[pairs[i,2]]) 
    cbind(pairs[i,],output)
  }  
  referenceScores <- testClonalityScores(pairs, eventMatrix, eventExpectation, scoreFun=scoreClonalityn) 
  clonalityScores$p.value <- sapply(clonalityScores$score, function(s) {mean(s <= c(referenceScores$score,s))} )
  
  return(clonalityScores)
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
testClonalityScores  <- function(pairs, eventMatrix, eventExpectation, scoreFun=scoreClonalityn) 
{
  if (is.null(pairs$patient))
  {
    pairs$patient <- 1:nrow(pairs)
  }
  refPairs <- expand.grid(list(Sample1 = unique(pairs[[1]]), Sample2 = unique(pairs[[2]])), stringsAsFactors = FALSE)
  refPairs <- refPairs[which(patients[refPairs$Sample1,1] != patients[refPairs$Sample2,1]),]
  return(ClonalityPairs(refPairs, eventMatrix, eventExpectation) )
}      



