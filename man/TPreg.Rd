\name{TPreg}
\alias{TPreg}
\title{
Regression Modeling of Transition Probabilities in a Progressive Illness-Death Model
}

\description{
 Fits a semi-parametric regression model to estimate the effects on
transition probabilities in a (possibly non-Markov) progressive illness-death model for a sequence of time.
}
\usage{
TPreg(formula, data, link, s = 0, t = NULL, R = 199, by = NULL, trans, ncores = NULL)
}

\arguments{
  \item{formula}{
an object of class formula which specifies the covariates. For example 
\code{formula}  \code{ = ~ age + sex}.
}
  \item{data}{
a data.frame of \code{iddata} class or a data.frame in which other than covariates five variables; \code{id}, \code{Zt}, \code{delta1}, \code{Tt}, \code{delta} are included. 
}
  \item{link}{
 a link function for \code{binomial} family which are \code{logit}, \code{probit} and \code{cauchit}, (corresponding to logistic, normal and Cauchy CDFs respectively). 
}
  \item{s}{
 the current time for the transition probabilities;  default is zero which reports the occupation probabilities. 
}
  \item{t}{
the Future time for the transition probabilities; default is \code{NULL} which is the largest uncensored sojourn time in the initial state. 
}
  \item{R}{
the number of bootstrap replicates. Default is 199. 
}
  \item{by}{
 number: increment of the sequence from time \code{s} to time \code{t}.The default is \code{NULL} which is
\eqn{\lfloor (\max({ Zt}) - \min({Zt}))/q_{0.01}({ Zt }) \rfloor} ,
where \eqn{q_{0.01}(.)} is the sample quantile corresponding to 0.01 probability and \eqn{\lfloor x\rfloor} gives
the largest integer less than or equal to \eqn{x}. A binomial regression at every \code{by}th time between \code{s} and 
\code{t} is performed. \code{by=1} reports all binomial regression results for each jump time,  corresponding to
 the specified transition(s),
between \code{s} and \code{t}.
 By increasing \code{by} we skip some times. In order to save the time, for a relatively large dataset
a relatively big \code{by} is recommended.
}
  \item{trans}{
the possible transition(s) for a progressive illness-death model. For \code{trans} argument there are five options:  \code{"11"},  \code{"12"},  \code{"13"}, \code{"23"}, and  \code{"all"}.
 }
  \item{ncores}{
the number of cores to use for parallel execution. Default is the number of CPU cores on the current host.
}
}

    \value{
\code{TPreg} returns an object of class {\code{TPreg}}.
 An object of class {\code{TPreg}} is a list containing at least the following components: 

  \item{co }{the list of:
    \itemize{
    \item \sQuote{transition} the transition,
    \item \sQuote{time} the jump times,
    \item \sQuote{coefficients} the estimated effects ,
    \item \sQuote{SD}   standard errors,
    \item \sQuote{LWL} lower confidence limits,
    \item \sQuote{UPL} upper confidence limits,
    \item \sQuote{p.value} p-values.}}
    \item{call }{the matched call.}
    \item{transition }{the  transition, this is equal to the \code{transition} in \code{co} unless for \code{trans="all"} }
    \item{s }{the current time for the transition probability.}
    \item{t }{the future time for the transition probability.}
    \item{n.misobs }{the number of missing observations.}
  
In addition, \code{trans="all"} will have four lists: \code{co11}, \code{co12}, \code{co13}, and \code{co23} instead of \code{co} and will give the information for all possible transitions in the progressive illness death model}

\references{
Azarang, L. Scheike, TH. and de Una-Alvarez, J. (2017) \emph{Direct modeling of regression effects for transition probabilities
in the progressive illness-death model}, Statistics in Medicine \bold{36}, \eqn{1964-1976}. 
}
\seealso{
 \code{\link{print.TPreg}}, \code{\link{summary.TPreg}}, and \code{\link{plot.TPreg}} as generic functions.
}

\examples{
data(colonTPreg)

co11 <- TPreg( ~ Age + Nodes + treatment, colonTPreg, link = "logit", s = 50, R = 19, t = 200,
trans = "11", ncores = 1)
co11
summary(co11)
}