\name{iddata}
\alias{iddata}

\title{
Transforming a Simple Illness-Death Dataset into a Proper Data Frame to be Used Inside \code{TPreg} Function
}
\description{
Converts a raw illness-death dataset to a data frame which can
be passed into \code{TPreg} function.
}
\usage{
iddata(Stime, Sind, Iltime, Ilind, covar, ...)
}
\arguments{
  \item{Stime}{the possibly censored total survival time.}
  \item{Sind}{the indicator of uncensored total survival time.}
  \item{Iltime}{the arrival time to the intermediate state, for example diseased, recurrence, or relapse time.}
  \item{Ilind}{the indicator of visiting intermediate state.}
  \item{covar}{a vector of covariates.}
  \item{\dots}{arguments to be passed to methods.}
}
\value{
  \code{iddata} returns an object of class data frame, containing at least the following components:
  \item{id }{ identification variable}
  \item{Zt }{the possibly censored sojourn time in initial state }
  \item{delta1 }{the indicator of uncensored sojourn time in initial state }
  \item{Tt }{the possibly censored total survival time}
  \item{delta }{the indicator of uncensored total survival time }
 \item{covar }{a vector of covariates}

}
\examples{
Stime <- c(10, 20, 34, 12, 30)
Sind <- c(1, 0, 0, 1, 1)
Iltime <- c(7, 20, 28, 12, 30)
Ilind <- c(1, 0, 1, 0, 0)
covar <- c(20, 25, 57, 36, 46)
mydata <- iddata(Stime, Sind, Iltime, Ilind, covar)  
mydata
}
