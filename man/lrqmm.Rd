\name{lrqmm}
\alias{lrqmm}

\title{
Fitting Linear Quantile Regression Mixed Models
}
\description{
Fit a quantile regression mixed model using a sparse implementation of the
  Frisch-Newton interior-point algorithm.
}
\usage{
lrqmm(id, sire, dam, X, Y, alpha = 0 , tau = 0.5 , Factor = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{ The number form animal record as vector or column matrix}
  \item{sire}{The number form father's animal record as vector or column matrix}
  \item{dam}{The number form mother's animal record as vector or column matrix}
  \item{X}{a fixed effect as vector or column matrix}
  \item{Y}{a response vector or column matrix}
  \item{alpha}{a parameter for raite error's varince to variance of random effects  }
  \item{tau}{desired quantile}
  \item{Factor}{type of fixed effect that "TRUE" as factor variable and "FALSE" as quantitative variable}
}
\details{
The function computes an estimate on the tau-th quantile effects of the linear mixed model. This is a sparse implementation of the Frisch-Newton algorithm for quantile regression described in Portnoy and Koenker (1997).

we used "GeneticsPed", "Matrix", "MasterBayes", "MCMCglmm", "MASS", "SparseM" and "quantreg" packages in this function. befor using "lrqmm" function be sure from installation this packages.
"GeneticsPed" available in <https://bioconductor.org/packages/release/bioc/src/contrib/GeneticsPed_1.46.0.tar.gz> or orders in <http://bioconductor.org/packages/release/bioc/html/GeneticsPed.html>. other packages are available in CRAN.
}
\value{
\item{Fixed effects}{estimate for fixed effect from linear quantile regression mixed models}
\item{Random effects}{estimate for random effect from linear quantile regression mixed models}
\item{summary}{reporting quantile for  effects estimation, number of observations and pedigree's length}
}
\references{
[1]Koenker, R. and S. Portnoy (1997) The Gaussian Hare and the Laplacean Tortoise: Computability of Squared-error vs Absolute Error Estimators, (with discussion). Statistical Science, 12, 279-300.
<https://www.jstor.org/stable/2246217>

[2]Koenker, R. W. (2005). Quantile Regression, Cambridge U. Press. ISBN: 0521608279.


[3]Mrode, R. A. (2005). Linear Models for the Prediction of Animal Breeding Values. 3rd edition.
CABI International. ISBN: 9781780643915.
}
\author{Sayyed Reza Alavian }



\examples{
data(Cow)
with(lrqmm(id,sire,dam,HERD,LACTLENGHT,2,0.5,Factor=TRUE),data=Cow)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{quantile regression}% use one of  RShowDoc("KEYWORDS")
\keyword{mixed models }% __ONLY ONE__ keyword per line
\keyword{Frisch-Newton interior-point }
