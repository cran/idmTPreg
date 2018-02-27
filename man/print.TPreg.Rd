\name{print.TPreg}
\alias{print.TPreg}
\title{
Printing the model fits
}
\description{
Provides the details about the estimated effect of covariate(s) on pre-specified
transition probabilities for given \code{s} and \code{t}.
}
\usage{
\method{print}{TPreg}(x, \dots)
}
\arguments{
  \item{x}{an object of class \code{TPreg}.}
  \item{\dots}{arguments to be passed to methods.}
}

\details{
\code{print.TPreg} formats the coefficients, standard errors, etc. 
The \code{coefficients} component of the result gives the estimated coefficients, their estimated standard errors,  their confidence bands, and the p-values. 
 The print method returns the aforementioned results for the specified \code{s} and \code{t}.
}
\value{
  \code{print.TPreg} returns a list with the following components:

    \item{Call}{ the matched call.}
    \item{Transition }{the specified transition(s).}
    \item{(s,t)}{ the vector of current time and future time for the transition probabilities. }
    \item{Coefficients}{ the matrix of coefficients, standard errors, lower confidence limits, upper confidence limits and p-values.}
 Also, \code{trans= "all"} will return \code{Coefficients 11}, \code{Coefficients 12}, \code{Coefficients 13}, and \code{Coefficients 23}.}
\seealso{
 \code{\link{TPreg}}, and  generic functions \code{\link{summary.TPreg}}, and \code{\link{plot.TPreg}}.
}
\examples{
## For examples see example TPreg.
}