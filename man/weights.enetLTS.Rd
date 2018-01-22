\name{weights.enetLTS}
\alias{weights.enetLTS}
%-------------------------------------------------
\title{
 binary weights from the \code{"enetLTS"} object
}
%-------------------------------------------------
\description{
Extract binary weights that indicate outliers from the current model.
}
%-------------------------------------------------
\usage{
\method{weights}{enetLTS}(object,vers=c("reweighted","raw","both"),index=FALSE,...)
}
%-------------------------------------------------
\arguments{
 \item{object}{the model fit from which to extract outlier weights.}

\item{vers}{a character string specifying for which estimator to extract
outlier weights. Possible values are \code{"reweighted"} (the default) for
weights indicating outliers from the reweighted fit, \code{"raw"} for
weights indicating outliers from the raw fit, or \code{"both"} for the
outlier weights from both estimators.}

\item{index}{a logical indicating whether the indices of the weight vector should
  be included or not (the default is \code{FALSE}).}

\item{\dots}{additional arguments from the \code{enetLTS} object if needed.}
}
%-------------------------------------------------
\value{
A numeric vector containing the requested outlier weights.
}
%-------------------------------------------------
\note{
The weights are \eqn{1} for observations with reasonably small
residuals and \eqn{0} for observations with large residuals.
Here, residuals represent standardized residuals
for linear regression and Pearson residuals for logistic residuals.


Use weights with or without index is available.
}
%-------------------------------------------------
\seealso{
\code{\link{enetLTS}}
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
weights(fit1)
weights(fit1,vers="raw",index=TRUE)
weights(fit1,vers="both",index=TRUE)}


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
weights(fit2)
}
\donttest{
fit2 <- enetLTS(xout,yout,family="binomial",alphas=0.5,lambdas=0.05,plot=FALSE)
weights(fit2)
weights(fit2,vers="raw",index=TRUE)
weights(fit2,vers="both",index=TRUE)}
}
%-------------------------------------------------
\author{
Fatma Sevinc KURNAZ, Irene HOFFMANN, Peter FILZMOSER
\cr Maintainer: Fatma Sevinc KURNAZ <fatmasevinckurnaz@gmail.com>;<fskurnaz@yildiz.edu.tr>}
%-------------------------------------------------
\keyword{regression}
\keyword{classification}


