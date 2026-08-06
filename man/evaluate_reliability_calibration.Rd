\name{evaluate_reliability_calibration}
\alias{evaluate_reliability_calibration}
\title{Evaluate reliability-model calibration}
\usage{
evaluate_reliability_calibration(
  model,
  features,
  target = c("first_pass_error", "refinement_benefit"),
  target_column = NULL,
  n_bins = 10
)
}
\arguments{
\item{model}{A model returned by \code{\link[=train_reliability_model]{train_reliability_model()}}.}

\item{features}{Reliability feature data frame with target columns.}

\item{target}{One of \code{"first_pass_error"} or \code{"refinement_benefit"}.}

\item{target_column}{Optional explicit binary target column.}

\item{n_bins}{Number of equal-frequency calibration bins.}
}
\value{
List with \code{metrics}, \code{bins}, and per-row \code{predictions}.
}
\description{
Computes Brier score, binary negative log-likelihood, expected calibration
error, calibration slope/intercept, and equal-frequency reliability bins for a
trained Risk-k model.
}
