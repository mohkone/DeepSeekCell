\name{get_reliability_spec}
\alias{get_reliability_spec}
\title{Get the frozen DeepSeekCell reliability specification}
\usage{
get_reliability_spec(include_marker_profiles = TRUE)
}
\arguments{
\item{include_marker_profiles}{Include marker profiles and aliases.}
}
\value{
Named list describing the frozen reliability specification.
}
\description{
Returns the named constants used by the frozen v1.0 reliability layer. The
companion JSON file in \code{inst/extdata/reliability_spec_v1.0.json} records
the same values together with external-validation provenance.
}
