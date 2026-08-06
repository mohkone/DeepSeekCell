\name{explain_reliability_model}
\alias{explain_reliability_model}
\title{Explain a learned reliability model}
\usage{
explain_reliability_model(model, conf_level = 0.95)
}
\arguments{
\item{model}{A model returned by \code{\link[=train_reliability_model]{train_reliability_model()}}.}

\item{conf_level}{Confidence level for Wald intervals.}
}
\value{
Data frame of feature-level model explanations.
}
\description{
Reports logistic coefficients, standardized coefficients, odds ratios and
Wald confidence intervals for a trained Risk-k model.
}
