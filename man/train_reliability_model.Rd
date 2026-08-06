\name{train_reliability_model}
\alias{train_reliability_model}
\title{Train a learned reliability model}
\usage{
train_reliability_model(
  features,
  target = c("first_pass_error", "refinement_benefit"),
  target_column = NULL,
  feature_columns = RELIABILITY_MODEL_FEATURES
)
}
\arguments{
\item{features}{Data frame from \code{\link[=extract_reliability_features]{extract_reliability_features()}}
or another table containing the required numeric feature columns.}

\item{target}{One of \code{"first_pass_error"} or
\code{"refinement_benefit"}.}

\item{target_column}{Optional explicit binary target column.}

\item{feature_columns}{Numeric feature columns used to train the model.}
}
\value{
Object of class \code{deepseekcell_reliability_model}.
}
\description{
Fits a transparent logistic model for either first-pass error risk or
expected refinement benefit. The returned object can be passed to
\code{select_refinement_candidates(strategy = "risk")}.
}
