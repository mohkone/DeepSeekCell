\name{compute_reliability_contributions}
\alias{compute_reliability_contributions}
\title{Compute linear SHAP-style reliability contributions}
\usage{
compute_reliability_contributions(model, features)
}
\arguments{
\item{model}{A model returned by \code{\link[=train_reliability_model]{train_reliability_model()}}.}

\item{features}{Reliability feature data frame.}
}
\value{
Long data frame of cluster-feature log-odds contributions.
}
\description{
For logistic regression, the per-feature log-odds contribution is computed as
\code{beta_j * (x_ij - mean_j)}. This is a transparent linear analogue of SHAP
for the fitted Risk-k model.
}
