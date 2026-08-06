\name{run_reliability_feature_ablation}
\alias{run_reliability_feature_ablation}
\title{Run reliability feature-ablation analysis}
\usage{
run_reliability_feature_ablation(
  training_features,
  evaluation_features = NULL,
  target = c("first_pass_error", "refinement_benefit"),
  feature_groups = NULL,
  budget = NULL
)
}
\arguments{
\item{training_features}{Development feature table.}

\item{evaluation_features}{Optional evaluation feature table. Defaults to
\code{training_features}.}

\item{target}{One of \code{"first_pass_error"} or
\code{"refinement_benefit"}.}

\item{feature_groups}{Named list of feature groups to remove.}

\item{budget}{Optional fixed k per block. If \code{NULL}, each block uses the
number of deterministic v1.0 \code{RequiresRefinement} clusters.}
}
\value{
Data frame comparing full and ablated models.
}
\description{
Trains full and feature-group-ablated Risk-k models, then evaluates
calibration and fixed-budget selection behavior on a validation feature table.
}
