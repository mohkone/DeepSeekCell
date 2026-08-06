\name{extract_reliability_features}
\alias{extract_reliability_features}
\title{Extract reliability-model features from first-pass annotations}
\usage{
extract_reliability_features(
  annotations,
  markers = NULL,
  tissue = NULL,
  truth = NULL,
  full_refined_annotations = NULL,
  use_ontology_evidence = TRUE
)
}
\arguments{
\item{annotations}{Annotation data frame.}

\item{markers}{Optional named marker list. When supplied and evidence columns
are absent, deterministic evidence scores are computed first.}

\item{tissue}{Optional tissue context.}

\item{truth}{Optional named truth vector keyed by cluster id. When supplied,
\code{InitiallyCorrect} and \code{FirstPassIncorrect} are added.}

\item{full_refined_annotations}{Optional FullRefined annotations used to
derive \code{FullRefinedCorrect} and \code{RefinementBeneficial}.}

\item{use_ontology_evidence}{Whether ontology evidence contributes if
evidence scores must be computed.}
}
\value{
Data frame with numeric reliability features and optional targets.
}
\description{
Converts first-pass annotations and deterministic evidence scores into a
numeric feature matrix suitable for learning first-pass error risk or
expected refinement benefit.
}
