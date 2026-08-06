\name{summarise_reliability_contributions}
\alias{summarise_reliability_contributions}
\title{Summarise global reliability contribution importance}
\usage{
summarise_reliability_contributions(contributions)
}
\arguments{
\item{contributions}{Output from \code{\link[=compute_reliability_contributions]{compute_reliability_contributions()}}.}
}
\value{
Feature-level mean absolute and signed contributions.
}
\description{
Summarises linear SHAP-style reliability contributions into global feature
importance scores.
}
