\name{enetLTS}
\alias{enetLTS}
%-------------------------------------------------
\title{
Robust and sparse estimation for linear and logistic regression}
%-------------------------------------------------
\description{
Compute fully robust versions of the elastic net estimator, which allows for sparse model estimates,
for linear and logistic regression.}
%-------------------------------------------------
\usage{
enetLTS(xx, yy, family=c("gaussian","binomial"),
alphas, lambdas, lambdaw, hsize=0.75,
intercept=TRUE, nsamp=500, s1=10, nCsteps=20, nfold=5,
seed=NULL, plot=TRUE, repl=5, para=TRUE, ncores=1,
del=0.0125, tol=-1e6, scal=TRUE, type=c("response","class"))}
%-------------------------------------------------
\arguments{
\item{xx}{a numeric matrix containing the predictor variables.}

\item{yy}{response variable. Quantitative for \code{family="gaussian"}. For
  \code{family="binomial"} should be a factor with two levels which is coded as \code{0} and {1}.}

\item{family}{a description of the error distribution and link function to be used
  in the model. \code{"gaussian"} and \code{"binomial"} options are available.}

\item{alphas}{a user supplied alpha sequence for the elastic net penalty, which is
  the mixing proportion of the ridge and lasso penalties and takes value in [0,1].
  \eqn{\alpha=1} is the lasso penalty, and \eqn{\alpha=0} the
  ridge penalty. If not provided a sequence, default is 41 equally spaced values.}

\item{lambdas}{a user supplied lambda sequence for the strength of the elastic net penalty.
  If not provided a sequence, default is chosen with steps of size -0.025 lambda0 with
  \eqn{0\le\lambda\le lambda0} for linear regression and
  -0.025 lambda00 with \eqn{0\le\lambda\le lambda00} for logistic regression. lambda0
  is determined based on the Pearson correlation between y and the jth predictor variable x_j
  on winsorized data for linear regression. In lambda00 for logistic regression, the Pearson
  correlation is replaced by a robustified point-biserial correlation.}

\item{lambdaw}{a user supplied lambda sequence for reweighting step. If not provided,
  default is computed by using k-fold cross-validation via \code{cv.glmnet} function.}

\item{hsize}{a user supplied numeric value giving the percentage of the residuals for
  which the elastic net penalized sum of squares for linear regression or for which the
  elastic net penalized sum of deviances for logistic regression should be minimized.
  The default is 0.75.}

\item{intercept}{a logical indicating whether a constant term should be
  included in the model (the default is \code{TRUE}).}

\item{scal}{a logical indicating whether scale the predictors by their arithmetic means
  and standard deviations. For \code{family="gaussian"}, it also indicates if
  mean-center the response variable or not. The default is \code{TRUE}. Note that scaling
  is performed on the subsamples rather than the full data set.}

\item{nsamp}{a numeric vector giving the number of subsamples to be used in
  the beginning of the algorithm, which gives the number of
  initial subsamples to be used. The default is to first perform C-steps on 500
  initial subsamples, and then to keep the \code{s1} subsamples with the lowest value (or highest   value based on which model is used - \code{"gaussian"} or \code{"binomial"}) of the objective    function   for additional C-steps until convergence.}

\item{s1}{a number of subsamples to keep after perform C-steps on \code{nsamp} initial subsets.
  For those remaining subsets, additional C-steps are performed until convergence. The
  default is 10.}

\item{nCsteps}{a positive integer giving the number of C-steps to perform on
  determined s1 subsamples. The default is 20.}

\item{nfold}{a user supplied numeric value for fold number of k-fold cross-validation which
  used in varied functions of the algorithm. The default is 5-fold cross-validation.}

\item{seed}{optional initial seed for the random number generator (see\code{\link{.Random.seed}}   ) when determine initial subsets at the beginning of the algorithm. The default is NULL.}

\item{plot}{a logical indicating if produces a plot for k-fold cross-validation based on alpha and lambda combinations. The default is TRUE.}

\item{repl}{a user supplied positive number for more stable results, repeat the k-fold CV        \code{repl} times and take the average of the corresponding evaluation measure. The default is 5.}

\item{para}{if TRUE, use parallel to fit each fold. Must register parallel before hand,
  such as doMC or others.}

\item{ncores}{a positive integer giving the number of processor cores to be
used for parallel computing (the default is 1 for no parallelization).  If
this is set to \code{NA}, all available processor cores are used.  For
prediction error estimation, parallel computing is implemented on the \R
level using package \pkg{parallel}.}

\item{del}{The default is 0.0125.}

\item{tol}{a small numeric value for convergence. The default is -1e6.}

\item{type}{type of prediction required. \code{type="response"} gives the
  fitted probabilities for \code{"binomial"} and  gives the fitted values for
  \code{"gaussian"}. \code{type="class"} is available only for \code{"binomial"} model,
  and produces the class label corresponding to the maximum probability.}
}
%-------------------------------------------------
\details{
The idea of repeatedly applying the non-robust classical elastic net estimators to data subsets
only is used for linear and logistic regression. The algorithm starts with 500 elemental subsets
only for one combination of \eqn{\alpha} and \eqn{\lambda}, and takes the warm start strategy
for subsequent combinations. This idea saves the computation time.
To choose the elastic net penalties, k-fold cross-validation is used and the replication option is
provided for more stable results.
Robustness has been achieved by using trimming idea, therefore a reweighting step is introduced
in order to improve the efficiency.
The outliers are identified according to current model.
For \code{family="gaussian"}, standardized residuals are used. For \code{family="binomial"}, the Pearson
residuals which are approximately standard normally distributed is used. Then the weights are defined by
the binary weight function using \code{del=0.0125}, which allows to be flagged as outliers of the
2.5\% of the observations in the normal model.
Therefore, binary weight function produces a clear distinction between the "good observations" and "outliers".}

%-------------------------------------------------
\value{
\item{objective}{a numeric vector giving the respective values of the
  enetLTS objective function, i.e., the elastic net penalized sums of
  the \eqn{h} smallest squared residuals from the raw fits for \code{family="gaussian"}
  and the elastic net penalized sums of the \eqn{h} deviances from the raw fits for
  \code{family="binomial"}.}

\item{best}{an integer vector containing the respective best
  subsets of \eqn{h} observations found and used for computing the raw
  estimates.}

\item{raw.wt}{an integer vector containing binary weights
  that indicate outliers from the respective raw fits, i.e., the weights used
  for the reweighted fits.}

\item{wt}{an integer vector containing binary weights that
  indicate outliers from the respective reweighted fits, i.e., the weights are
  \eqn{1} for observations with reasonably small reweighted residuals and
  \eqn{0} for observations with large reweighted residuals.}

\item{a00}{intercept term obtained from the raw fit.}

\item{raw.coefficients}{a numeric vector containing the
  respective coefficient estimates from the raw fit.}

\item{a0}{intercept term obtained from the reweighted fit.}

\item{coefficients}{a numeric vector containing the
  respective coefficient estimates from the reweighted fit.}

\item{alpha}{an optimal elastic net mixing parameter value obtained with
  k-fold cross-validation.}

\item{lambda}{an optimal value for the strength of the elastic net penalty
  obtained with k-fold cross-validation.}

\item{lambdaw}{an optimal value for the strength of the elastic net penalty
  re-obtained with k-fold cross-validation for reweighted fit.}

\item{num.nonzerocoef}{the number of the nonzero coefficients in the model.}

\item{h}{the number of observations used to compute the raw estimates.}

\item{raw.residuals}{a numeric vector containing the
  respective residuals from the raw fits.}

\item{residuals}{a numeric vector containing the
  respective residuals from the reweighted fits.}

\item{raw.fitted.values}{a numeric vector containing the
  respective fitted values of the response from the raw fits.}

\item{fitted.values}{a numeric vector containing the
  respective fitted values of the response from the reweighted fits.}

\item{raw.rmse}{root mean squared error for raw fit, which is available for only
\code{family="gaussian"}.}

\item{rmse}{root mean squared error for reweighted fit, which is available for only
\code{family  ="gaussian"}.}

\item{classnames}{class names for logistic model, which is available for only
\code{family="binomial"}.}

\item{classize}{class sizes for logisitic model, which is available for only
\code{family="binomial"}.}

\item{inputs}{all inputs used in the function \code{enetLTS.R}.}

\item{call}{the matched function call.}
}
%-------------------------------------------------
\references{
Kurnaz, F.S., Hoffmann, I. and Filzmoser, P. (2017) Robust and sparse estimation methods
for high dimensional linear and logistic regression. \emph{Chemometrics and Intelligent Laboratory Systems.}
}
%-------------------------------------------------
\author{
Fatma Sevinc KURNAZ, Irene HOFFMANN, Peter FILZMOSER
}
%---------------------------------------------------
\seealso{
\code{\link{print}},
\code{\link{predict}},
\code{\link{coef}},
\code{\link{nonzeroCoef.enetLTS}},
\code{\link{plot}},
\code{\link{plotCoef.enetLTS}},
\code{\link{plotResid.enetLTS}},
\code{\link{plotDiagnostic.enetLTS}},
\code{\link{residuals}},
\code{\link{fitted}},
\code{\link{weights}}
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
fit <- enetLTS(xout,yout,alphas=0.5,lambdas=0.05,plot=FALSE)
# determine user supplied alpha and lambda sequences
# alphas=seq(0,1,length=11)
# l0 <- robustHD::lambda0(xout,yout)          # use # lambda0 function from robustHD package
# lambdas <- seq(l0,0,by=-0.1*l0)
# fit <- enetLTS(xout,yout,alphas=alphas,lambdas=lambdas)}

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
}
\donttest{
	fit <- enetLTS(xout,yout,family="binomial",alphas=0.5,lambdas=0.05,plot=FALSE)
# determine user supplied alpha and lambda sequences
# l00 <- lambda00(xout,yout,normalize=TRUE,intercept=TRUE)
# lambdas <-  seq(l00,0,by=-0.01*l00)
# fit <- enetLTS(xout,yout,family="binomial",alphas=alphas,lambdas=lambdas)}
}
%-------------------------------------------------
\keyword{regression}
\keyword{classification}
\keyword{sparse}
\keyword{robust}



