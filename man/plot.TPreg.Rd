\name{plot.TPreg}
\alias{plot.TPreg}

\title{
Plotting a TPreg object
}
\description{
 Makes a plot for the estimated effect of covariate(s) on  pre-specified transition probabilities together with \eqn{95\%}
confidence bands along time, from time \code{s} to time \code{t}.
}
\usage{
\method{plot}{TPreg}(x, covar, rug = TRUE, main, ylab, xlab, Ylim, \ldots)
}
\arguments{
  \item{x}{an object of class \code{TPreg}}
  \item{covar}{a character vector of the name(s) of covariate(s)}
  \item{rug}{if \code{TRUE} (default), it adds a rug representation of times between time \code{s} and time \code{t}. }
  \item{main}{an overall title for the plot}
  \item{ylab}{a title for the y axis}
  \item{xlab}{a title for the x axis}
  \item{Ylim}{a list of y limits}
  \item{\dots}{arguments to be passed to methods.}
}

\seealso{
 \code{\link{TPreg}}, and generic functions \code{\link{print.TPreg}}, and \code{\link{summary.TPreg}}.
}

\examples{
data(colonTPreg)

co13 <- TPreg( ~ Age + Nodes + treatment, colonTPreg, link = "logit", s = 0, R = 9, t = 1200,
  by = 110, trans = "13", ncores = 1)
plot(co13, covar = "Nodes", Ylim = list(c(-0.5,0.5)))

 \donttest{
co11 <- TPreg( ~ Age + Nodes + treatment, colonTPreg, link = "logit", s = 0, by = 10,
  R = 199, t = 365*5, trans = "11")
plot(co11, covar = c("Age", "Nodes", "treatmentLev", "treatmentLev.5FU"),
Ylim = list(c(-0.1,0.1), c(-0.5,0.5), c(-2,2), c(-2,2)))
}

}