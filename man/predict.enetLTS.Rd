\name{predict.enetLTS}
\alias{predict.enetLTS}
%-------------------------------------------------
\title{
make predictions from the \code{"enetLTS"} object.
}
%-------------------------------------------------
\description{
Similar to other predict methods, this function predicts fitted values, logits,
coefficients and nonzero coefficients from a fitted \code{"enetLTS"} object.
}
%-------------------------------------------------
\usage{
\method{predict}{enetLTS}(object,newX,vers=c("reweighted","raw","both"),
    type=c("response","coefficients","nonzero","class"),...)
}
%------------------------------------
\arguments{
\item{object}{the model fit from which to make predictions.}

\item{newX}{new values for the predictor matrix \code{X}.
  Must be a matrix; can be sparse as in \code{Matrix} package.
  This argument is not used for \code{type=c("coefficients","nonzero")}.}

\item{vers}{a character string denoting which fit to use for the predictions.
  Possible values are \code{"reweighted"} (the default) for
  predicting values from the reweighted fit, \code{"raw"} for predicting
  values from the raw fit, or \code{"both"} for predicting values from both
  fits.}

\item{type}{type of prediction required. \code{type="response"} gives the
  fitted probabilities for \code{"binomial"} and gives the fitted values for
  \code{"gaussian"}. \code{type="coefficients"} computes the coefficients from the
  fitted model. \code{type="nonzero"} returns a list of the indices of the nonzero
  coefficients. \code{type="class"} is available only for \code{"binomial"} model,
  and produces the class label corresponding to the maximum probability.}

\item{\dots}{additional arguments from the \code{enetLTS} object if needed.}
}
%-------------------------------------------------
\details{
The \code{newdata} argument defaults to the matrix of predictors used to fit
the model such that the fitted values are computed.

\code{coef.enetLTS(...)} is equivalent to \code{predict.enetLTS(object,newX,type="coefficients",...)}, where newX argument is the matrix as in \code{enetLTS}.
}
%-------------------------------------------------
\value{
The requested predicted values are returned.
}
%-------------------------------------------------
\seealso{
\code{\link{enetLTS}},
\code{\link{coef.enetLTS}},
\code{\link{nonzeroCoef.enetLTS}}
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
predict(fit1,newX=xout)
predict(fit1,newX=xout,type="coefficients",vers="both")
predict(fit1,newX=xout,type="nonzero",vers="raw")
# provide new X matrix
newX <- matrix(rnorm(n*p, sigma),nrow=n)
predict(fit1,newX=newX,type="response",vers="both")
predict(fit1,newX=newX,type="coefficients")
predict(fit1,newX=newX,type="nonzero",vers="both")}


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
predict(fit2,newX=xout)
}
\donttest{
fit2 <- enetLTS(xout,yout,family="binomial",alphas=0.5,lambdas=0.05,plot=FALSE)
predict(fit2,newX=xout)
predict(fit2,newX=xout,type="coefficients",vers="both")
predict(fit2,newX=xout,type="nonzero",vers="raw")
predict(fit2,newX=newX,type="class",vers="both")
predict(fit2,newX=newX,type="coefficients",vers="raw")
predict(fit2,newX=newX,type="nonzero",vers="both")}

}
%-------------------------------------------------
\author{
Fatma Sevinc KURNAZ, Irene HOFFMANN, Peter FILZMOSER
\cr Maintainer: Fatma Sevinc KURNAZ <fatmasevinckurnaz@gmail.com>;<fskurnaz@yildiz.edu.tr>}
%-------------------------------------------------
\keyword{regression}
\keyword{classification}


