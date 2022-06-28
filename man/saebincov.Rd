\name{SAE EB Beta Binomial Covariates}
\alias{saebincov}
\title{SAE EB Beta Binomial Covariates}
\usage{
saebincov(y, non_y, n, x, N, area, unit, data)
}
\arguments{
  \item{y}{
  response variable contains the numbers of cases in each area.
}
  \item{non_y}{
  sample numbers in each area-numbers of cases in each area.
}
  \item{n}{
  sample numbers in each area.
}
  \item{x}{
  auxuliary variabel.
}
  \item{N}{
  populatio number in each area.
}
  \item{area}{
  id code area.
}
  \item{unit}{
  id code unit.
}
  \item{data}{
  data frame to use.
}
}
\description{
This function gives the area level EB and MSE estimator to binomial model respone with auxuliary variabel.
}
\examples{
\dontrun{
saebincov(y, non_y, n, x, N, area, unit, data)
}
}
