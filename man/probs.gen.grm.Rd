% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robust.R
\name{probs.gen.grm}
\alias{probs.gen.grm}
\title{Response Probability Calculation (2PL) and response category probabilities (P)}
\usage{
probs.gen.grm(thetas, a, b)
}
\arguments{
\item{thetas}{An array of person abilities}

\item{a}{A matrix of slope parameters}

\item{b}{Array of intercept parameters}
}
\value{
list of threshold probabilities (p-star) and the corresponding response category probabilities (P)
}
\description{
This function generates threshold probabilities (p-star) and response category probabilities (P)
}