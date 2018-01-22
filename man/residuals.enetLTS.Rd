\name{residuals.enetLTS}
\alias{residuals.enetLTS}
%-------------------------------------------------
\title{
the residuals from the \code{"enetLTS"} object
}
%-------------------------------------------------
\description{
A numeric vector which returns residuals from the enetLTS object.
}
%-------------------------------------------------
\usage{
\method{residuals}{enetLTS}(object,vers=c("reweighted","raw","both"),...)
}
%-------------------------------------------------
\arguments{
\item{object}{the model fit from which to extract residuals.}

\item{vers}{a character string specifying for which estimator to extract
outlier weights.  Possible values are \code{"reweighted"} (the default) for
weights indicating outliers from the reweighted fit, \code{"raw"} for
weights indicating outliers from the raw fit, or \code{"both"} for the
outlier weights from both estimators.}

\item{\dots}{additional arguments from the enetLTS object.}
}

%-------------------------------------------------
\value{
A numeric vector containing the requested residuals.
}
%-------------------------------------------------
\seealso{
\code{\link{enetLTS}},
\code{\link{fitted.enetLTS}},
\code{\link{predict.enetLTS}},
\code{\link{coef.enetLTS}}
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
residuals(fit1)
residuals(fit1,vers="raw")
residuals(fit1,vers="both")}


## for binomial

eps <-0.05                                    # \%10 contamination to only class 0
m <- ceiling(eps*n)
y <- sample(0:1,n,replace=TRUE)
xout <- x
xout[y==0,][1:m,] <- xout[1:m,] + 10;         # class 0
yout <- y                                     # wrong classification for vertical outliers

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
residuals(fit2)
}
\donttest{
fit2 <- enetLTS(xout,yout,family="binomial",alphas=0.5,lambdas=0.05,plot=FALSE)
residuals(fit2)
residuals(fit2,vers="raw")
residuals(fit2,vers="both")}
}
%-------------------------------------------------
\author{
Fatma Sevinc KURNAZ, Irene HOFFMANN, Peter FILZMOSER
\cr Maintainer: Fatma Sevinc KURNAZ <fatmasevinckurnaz@gmail.com>;<fskurnaz@yildiz.edu.tr>}
%-------------------------------------------------
\keyword{regression}
\keyword{classification}

