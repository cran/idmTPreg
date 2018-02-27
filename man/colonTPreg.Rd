\name{colonTPreg}
\alias{colonTPreg}
\docType{data}
\title{
The colon cancer data
}
\description{
Colon cancer data  with 929 rows and 8 columns. There is one row per person.
The original data have 16 columns and two records per person and are available as part of \code{survival} Package.
These are data from one of the first successful trials of adjuvant chemotherapy for colon cancer.
Levamisole is a low-toxicity compound previously used to treat worm infestations in animals; 
5-FU is a moderately toxic (as these things go) chemotherapy agent.  
}
\usage{data("colonTPreg")}
\format{
  A data frame with 929 observations on the following 8 variables.
  \describe{
    \item{\code{id}}{identification variable.}
    \item{\code{Zt}}{recurrence free survival time, that is the minimum of recurrence, death, and censoring time.}
    \item{\code{Tt}}{total survival time, time to death or censoring}
    \item{\code{delta1}}{recurrence free survival indicator,  0=alive and recurrence free, 1=dead or relapsed.}
    \item{\code{delta}}{status indicator,  0=alive, 1=dead.}
    \item{\code{Nodes}}{ number of lymph nodes with detectable cancer.}
    \item{\code{Age}}{ age in years.}
    \item{\code{treatment}}{a factor with levels \code{Obs}, \code{Lev} and \code{Lev+5FU}}
  }
}