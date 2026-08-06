\name{predict_reliability_risk}
\alias{predict_reliability_risk}
\title{Predict reliability risk or expected refinement benefit}
\usage{
predict_reliability_risk(model, features)
}
\arguments{
\item{model}{A model returned by
\code{\link[=train_reliability_model]{train_reliability_model()}}.}

\item{features}{Data frame containing reliability feature columns.}
}
\value{
Numeric vector of predicted risk/benefit probabilities in \code{[0, 1]}.
}
\description{
Predicts first-pass error risk or expected refinement benefit for candidate
annotations using a learned DeepSeekCell reliability model.
}
