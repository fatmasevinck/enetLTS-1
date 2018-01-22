
###################################################################################################################
###################################################################################################################
prepara <- function(x,y,family,index=NULL,robu=NULL){
  # default clasical scale, robust scale for robu=1

  if (is.null(robu)) robu=0
  if (is.null(index)) {
    if (robu>0){
      if(family=="binomial"){
        muy=y
      } else if(family=="gaussian"){
        muy <- median(y)
      }
      mux <- apply(x,2,median)   ## to go back original coef
      sigx <- apply(x,2,mad)
    } else {
      if(family=="binomial"){
        muy=y
      } else if(family=="gaussian"){
        muy <- mean(y)
      }
      mux <- apply(x,2,mean)   ## to go back original coef
      sigx <- apply(x,2,sd)
    }
  } else {
    if (robu>0){
      if(family=="binomial"){
        muy=y
      } else if(family=="gaussian"){
        muy <- median(y[index])
      }
      mux <- apply(x[index,],2,median)
      sigx <- apply(x[index,],2,mad)
    } else {
      if(family=="binomial"){
        muy=y
      } else if(family=="gaussian"){
        muy <- mean(y[index])
      }
      mux <- apply(x[index,],2,mean)
      sigx <- apply(x[index,],2,sd)
    }
  }
  xnor <- scale(x,mux,sigx)
  if(family=="binomial"){
    ycen <- y}
  else if(family=="gaussian"){
    ycen <- scale(y,muy,FALSE)}
  return(list(xnor=xnor,ycen=ycen,mux=mux,sigx=sigx,muy=muy))
}

################################################################################################################################
weight.gaussian <- function(resi,ind,del){
  if(is.logical(ind)){
    h <- length(which(ind==TRUE))
  }else{
    h <- length(ind)
  }
  n <- length(resi)
  mu <- mean(resi[ind])
  rc <- (resi - mu)
  qn <- qnorm((h+n)/ (2*n))                         # required quantile
  cdelta <- 1 / sqrt(1 - (2*n)/(h/qn) * dnorm(qn))
  s <- sqrt(mean(rc[ind]^2)) * cdelta
  we <- as.integer(abs(rc/s) <= qnorm(1-del))
  out <- list(we=we,mu=mu,s=s)

  return(out)
}


weight.binomial <- function(x,y,beta,intercept,del){
  if(intercept==TRUE){
    pi <- exp(cbind(1,x)%*%beta)/(1+exp(cbind(1,x)%*%beta))
    res <- (y - pi) / sqrt(pi*(1-pi))
  } else{
    pi <- exp(x%*%beta[-1])/(1+exp(x%*%beta[-1]))
    res <- (y - pi) / sqrt(pi*(1-pi))
  }
  we <- as.integer(abs(res) <= qnorm(1-del))
  return(we)
}


addColnames <- function(x) {
  # 'x' needs to be a matrix
  if(is.null(colnames(x))) colnames(x) <- paste("x", seq_len(ncol(x)), sep="")
  x
}

## add intercept column to design matrix
addIntercept <- function(x, check = FALSE) {
  if(!check || is.na(match("(Intercept)", colnames(x)))) {
    cbind("(Intercept)"=rep.int(1, nrow(x)), x)
  } else x
}


uptrimMSE<- function(x,trim=0.1){
  # computes trim% upper trimmed mean
  return(mean(x[x<quantile(x,1-trim)]))
}

## use in lambda00
winsorize.default <- function(x, standardized = FALSE, centerFun = median,
                              scaleFun = mad, const = 2,
                              return = c("data", "weights"), ...) {
   ## initializations
   standardized <- isTRUE(standardized)
   if(standardized) return <- match.arg(return)
   else {
      # standardize data
      x <- robStandardize(x, centerFun=centerFun, scaleFun=scaleFun, ...)
      center <- attr(x, "center")
      scale <- attr(x, "scale")
   }
   ## winsorize standardized data
   #   ind <- abs(x) > const           # observations in 'x' that need to be shrunken
   #   x[ind] <- const * sign(x[ind])  # winsorize
   weights <- pmin(const / abs(x), 1)
   if(standardized && return == "weights") return(weights)
   x <- weights * x
   ## finalizations
   if(!standardized) {
      # transform back to original scale and remove attributes
      x <- c(x * scale + center)
   }
   x
}



