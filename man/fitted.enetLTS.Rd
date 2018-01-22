\name{fitted.enetLTS}
\alias{fitted.enetLTS}
%-------------------------------------------------
\title{
 the fitted values from the \code{"enetLTS"} object.
}
%-------------------------------------------------
\description{
  A numeric vector which extract fitted values from the current model.
}
%-------------------------------------------------
\usage{
  \method{fitted}{enetLTS}(object,vers=c("reweighted","raw","both"),type=c("response","class"),...)
}
%-------------------------------------------------
\arguments{
\item{object}{the model fit from which to extract fitted values.}

\item{vers}{a character string specifying for which fit to make
  predictions.  Possible values are \code{"reweighted"} (the default) for
  predicting values from the reweighted fit, \code{"raw"} for predicting
  values from the raw fit, or \code{"both"} for predicting values from both
  fits.}

\item{type}{type of prediction required. \code{type="response"} gives the
  fitted probabilities for \code{"binomial"} and gives the fitted values for
  \code{"gaussian"}. \code{type="class"} is available only for \code{"binomial"} model,
  and produces the class label corresponding to the maximum probability.}

\item{\dots}{additional arguments from the \code{enetLTS} object if needed.}
}
%-------------------------------------------------
\value{
A numeric vector containing the requested fitted values.
}
%-------------------------------------------------
\seealso{
\code{\link{enetLTS}},
\code{\link{predict.enetLTS}},
\code{\link{residuals.enetLTS}}
}
%-------------------------------------------------
\examples{
## for gaussian

set.seed(86)
n <- 100; p <- 25                             # number of observations and variables
beta <- rep(0,p); beta[1:6] <- 1              # 10\% nonzero coefficients
sigma <- 0.5                                  # controls signal-to-noise ratio
x <- matrix(rnorm(n*p, sigma),nrow=n)
e <- rnorm(n,0,1)                             # error terms
eps <- 0.1                                    # contamination level
m <- ceiling(eps*n)                           # observations to be contaminated
eout <- e; eout[1:m] <- eout[1:m] + 10        # vertical outliers
yout <- c(x \%*\% beta + sigma * eout)          # response
xout <- x; xout[1:m,] <- xout[1:m,] + 10      # bad leverage points

\donttest{
fit1 <- enetLTS(xout,yout,alphas=0.5,lambdas=0.05,plot=FALSE)
fitted(fit1)
fitted(fit1,vers="raw")
fitted(fit1,vers="both")
fitted(fit1,vers="reweighted",type="response")}

## for binomial
eps <-0.05                                     # \%10 contamination to only class 0
m <- ceiling(eps*n)
y <- sample(0:1,n,replace=TRUE)
xout <- x
xout[y==0,][1:m,] <- xout[1:m,] + 10;          # class 0
yout <- y                                      # wrong classification for vertical outliers

\dontshow{
set.seed(86)
n <- 5; p <- 15                            
beta <- rep(0,p); beta[1:6] <- 1              
sigma <- 0.5                        
x <- matrix(rnorm(n*p, sigma),nrow=n)   
e <- rnorm(n,0,1)                             # error terms
eps <- 0.1                                    # contamination level
m <- ceiling(eps*n)                           # observations to be contaminated
eout <- e; eout[1:m] <- eout[1:m] + 10        # vertical outliers
yout <- c(x \%*\% beta + sigma * eout)          # response
xout <- x; xout[1:m,] <- xout[1:m,] + 10      # bad leverage points
fit2 <- enetLTS(xout,yout,alphas=0.5,lambdas=0.05,plot=FALSE)
fitted(fit2)
}
\donttest{
fit2 <- enetLTS(xout,yout,family="binomial",alphas=0.5,lambdas=0.05,plot=FALSE)
fitted(fit2)
fitted(fit2,vers="raw")
fitted(fit2,vers="both",type="class")
fitted(fit2,vers="both")
fitted(fit2,vers="reweighted",type="class")}
}
\author{
Fatma Sevinc KURNAZ, Irene HOFFMANN, Peter FILZMOSER
\cr Maintainer: Fatma Sevinc KURNAZ <fskurnaz@gmail.com>;<fskurnaz@yildiz.edu.tr>}

\keyword{regression}
\keyword{classification}


