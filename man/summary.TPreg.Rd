\name{summary.TPreg}
\alias{summary.TPreg}
\title{
Summarizing Details of the Estimated Effects on the Transition Probabilities in a Progressive Illness-Death Model
}
\description{
 Gives details about the estimated effect of covariate(s) on  pre-specified
transition probabilities  for a progressive illness-death model for a sequence of time from a given \code{s} to a given \code{t}.
}
\usage{
\method{summary}{TPreg}(object, \dots)
}
\arguments{
  \item{object}{ an object of class \code{TPreg}.}
  \item{\dots}{arguments to be passed to methods.}
}
\details{
\code{summary.TPreg} formats the coefficients, standard errors, etc. for  the jump times already fixed by 's', 't',  'by' arguments, corresponding to the possible transition. The jump times come out in the 
first column of these components. 
}
\value{
\code{summary.TPreg} returns a list with at least the following components:

 \item{Call }{the matched call.}
 \item{Transition(s) }{the specified transition(s).}
 \item{(s,t)}{the vector of current time and future time for the transition probabilities.}
 \item{Coefficients }{the matrix of selected jump times and  their corresponding coefficients.}
 \item{Standard Errors }{the matrix of selected jump times and their corresponding standard errors.}
 \item{Lower limits }{the matrix of selected jump times and their corresponding lower confidence bands.}
 \item{Upper limits }{the matrix of selected jump times and their corresponding upper confidence bands.}
 \item{p.values }{the matrix of selected jump times and their corresponding p-values.}
 
for \code{Transitions = "11", "12", "13", "23"},  \code{summary.TPreg} gives \code{Coefficients}, \code{Standard Errors}, \code{Lower limits}, \code{Upper limits}, and \code{p-values} for each of \code{Transition 11}, \code{Transition 12}, \code{Transition 13} and \code{Transition 23}. }

\seealso{
 \code{\link{TPreg}},  and  generic functions \code{\link{print.TPreg}} and \code{\link{plot.TPreg}}.
}
\examples{
## For examples see example TPreg.
}