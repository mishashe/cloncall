% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClonalityPackage.R
\name{clonality}
\alias{clonality}
\title{Compute clonality score and p-values of paired samples.}
\usage{
clonality(
  pairs,
  eventDataFrame,
  identityColumns,
  samples,
  ditanceTolerance = NULL
)
}
\arguments{
\item{pairs}{A dataframe of sample pairs with sample identifiers in columns sample1 and sample2.
Column Patient is optional. If it is not available it is assumed all pairs are from different
patients.}

\item{eventDataFrame}{table of mutations in the dataframe format}

\item{identityColumns}{relevant columns}

\item{samples}{samples to analyse}

\item{ditanceTolerance}{Distance tolerance to cluster mutations}
}
\value{
pairs with extra columns scorei, number of shared breaks and p-value
}
\description{
Compute clonality score and p-values of paired samples.
}
