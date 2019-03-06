\name{SmoothECTest}
\alias{SmoothECTest}
\title{
Smooth Goodness-of-fit Test for Multivariate Elliptical Distributions
}
\description{
Smooth tests of goodness-of-fit for multivariate elliptical distributions with diagnostic (Dx) capabilities and
full invariance to affine-linear transformations. By increasing the value of the hyperparameter \code{K}, 
the test and the Dx become adaptively consistent against an increasing number of departures from 
the null model. The Dx pertains to elements \eqn{R} and \eqn{U} of the Cambanis, Huang & Simons 
stochastic representation of elliptical data. Note that p-values can be computed via an asymptotic chi-square approximation or by Monte Carlo.
}
\usage{
SmoothECTest(data, K = 7, family = "MVN", Est.Choice = "", Cpp = TRUE)
}
\arguments{
  \item{data}{
The data set to use. Cases with missing values are removed.
}
  \item{K}{
Integer. Hyperparameter controlling the size of the embedding
family. Should be greater than or equal to 3
for the Multivariate Normal Distribution. The computation time increases
with the size of the data frame and \code{K} . Please be patient.
}
  \item{family}{
The only family available in the current version of the package is the
Multivariate Normal Distribution.
}
\item{Est.Choice}{
Not used yet. Maximum Likelihood Estimation (MLE) or Method of moments. Currently, only the MLE is implemented.
}
  \item{Cpp}{
Logical. If \code{TRUE}, the faster C++ code is used. 
}
}
\value{List with components:
  \item{Q}{The global test statistic with hyperparameter \code{K}.
  }
  \item{dfQ}{Degrees of freedom of the asymptotic chi-square approximation.
  }
  \item{pval.asymp.Q}{Asymptotic p-value for \eqn{Q}.
  }
  \item{Uscaled}{Scaled component \eqn{U(s)_K} tests the uniformity of element \eqn{U}.
  }
  \item{dfU}{Degrees of freedom of the asymptotic chi-square approximation.
  }
  \item{pval.asymp.U}{Asymptotic p-value for \code{Uscaled}.
  }
  \item{Iscaled}{Scaled component \eqn{I(s)_K} tests the correlation between \eqn{R} and \eqn{U}.
  }
  \item{dfI}{Degrees of freedom of the asymptotic chi-square approximation.
  }
  \item{pval.asymp.I}{Asymptotic p-value for \code{Iscaled}.
  }
  \item{Rscaled}{Scaled component \eqn{R(s)_K} test the distribution of element \eqn{R}, radius of the data.
  }
  \item{dfR}{Degrees of freedom of the asymptotic chi-square approximation.
  }
  \item{pval.asymp.R}{Asymptotic p-value for \code{Rscaled}.
  }
}
\references{
Gilles R. Ducharme, Pierre Lafaye de Micheaux (2019). A Goodness-of-fit
Test for Elliptical Distributions with Diagnostic Capabilities.
ArXiv pre-print. \url{https://arxiv.org/abs/1902.03622}
}
\author{
G. R. Ducharme, P. Lafaye De Micheaux
}
\examples{
 # The famous (Fisher's or Anderson's) iris data set
 # Increase the value of K to K = 7 for better results.
  ressetosa <- SmoothECTest(iris[1:50, -5], K = 3)
  ressetosa

 # Examination marks (n = 88) in Vectors, Algebra and Statistics from the "Open
 # book-Closed book examination" data set (Mardia, Kent and Bibby, 1979,
 # p. 3-4).
 # Increase the value of K to K = 5 for better results.
  data <- scor[, c(2, 3, 5)]
  result <- SmoothECTest(data, K = 3)
  result
}
\keyword{distribution}
\keyword{multivariate}
\keyword{htest}
