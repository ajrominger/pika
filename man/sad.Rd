% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sad.R
\name{sad}
\alias{sad}
\title{Create a species abundance distribution}
\usage{
sad(x, model, par, keepData = TRUE)
}
\arguments{
\item{x}{vector of integers representing a sample of species abundances}

\item{model}{character naming model to be fit (can only be length one)}

\item{par}{vector of model parameters}

\item{keepData}{logical, should the data be stored with the fitted \code{sad} object}
}
\value{
A list containing named objects of class \code{sad} with elements
\describe{
  \item{\code{MLE}}{The maximum likelihood estimate(s)}
  \item{\code{ll}}{The maximized log likelihood}
  \item{\code{df}}{The associated degrees of freedom}
  \item{\code{nobs}}{The associated number of observations}
  \item{\code{model}}{The name of the model fit}
  \item{\code{data}}{The data}
}
}
\description{
Function to create an object of class \code{sad} from data or theoretical predictions
}
\details{
An object of class \code{sad} can be created one three ways:
\describe{
 \item{Raw data}{specifying only \code{x} will produce an \code{sad} object without any parametric model fit to it}
 \item{Fit model}{specifying both \code{x} and \code{model} will fit the model to the data and return an \code{sad} object with both empirical and theoretical predictions}
 \item{Theoretical model}{specifying both \code{model} and \code{par} will produce an \code{sad} object with only the theoretical prediction but no data}
}

Note that this function behaves almost identically to \code{fitSAD} for the case when both \code{x} and \code{model} are spcified, but \code{sad} is more flexible allowing for unfit \code{sad} objects and specification of theoretical SADs without fitting of data by maximum likelihood.
}
\examples{
x <- rtpois(1000, 2)
sad(x, model='tpois')

}
\seealso{
logLik.sad, dlseries, dplnorm, dstick, dtnegb, dtpois
}
\author{
Andy Rominger <ajrominger@gmail.com>
}
\keyword{Maximum}
\keyword{SAD}
\keyword{abundance,}
\keyword{likelihood,}
\keyword{species}
