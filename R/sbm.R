#' Comprehensive structural change estimation and testing
#' @description
#' `sbm()` is used to estimate linear models with structural breaks, including both pure and
#' partial model with pre-specified number of breaks given a pre-specified trimming level. It
#' can be used to carry out inference for each regime, testing and selecting number of
#' structural breaks using information criteria.
#' 
#' @usage sbm(formula, data, m = 5, eps = 0.15)
#' 
#' 
#' @param formula an object of class formula: a symbolic description of the pure/partial
#' structural break model to be fitted. The details of model specification are given 
#' under 'Details' 
#' @param data a data frame, or list containing the variables in the model. 
#' @param eps1 trimming level (in percentage) influencing the minimal segment length that is 
#' equivalent to \code{default} = int(\code{eps1}*T) (T is total sample size).
#  There are five options:
#' \itemize{
#' \item{`eps1=0.05` corresponds to maximum 9 possible breaks}
#' \item{`eps1=0.10` corresponds to maximum 8 possible breaks}
#' \item{`eps1=0.15` corresponds to maximum 5 possible break}
#' \item{`eps1=0.20` corresponds to maximum 3 possible break}
#' \item{`eps1=0.25` corresponds to maximum 2 possible break}
#' } If users want to specify their own combinations of number of breaks `m` and minimum segment length `h`, use [estimate()]
#' @param cov. There are six parameters users can specify for covariance matrix:
#' \itemize{
#'  \item{`prewhit` sets to \code{1} to apply AR(1) prewhitening prior to estimating
#' the long run covariance matrix. The method used is \emph{Andrews(1991)} automatic bandwidth with AR(1) approximation with quadratic
#' kernel. Note: Do not set to \code{1} if lagged dependent variables are
#' included as regressors.}
#' \item{`robust` sets to \code{1} to allow for heterogeneity and autocorrelation in the residuals, \code{0} otherwise.}
#' \item{`hetdat` is used for the construction of the F tests. Set to \code{1} if users want to allow different moment 
#' matrices of the regressors across segments. If \code{hetdat} = \code{0}, the same moment matrices are assumed for each segment
#' and estimated from the ful sample. It is recommended to set \code{hetdat}=\code{1} if number of regressors \code{x} > \code{0}.}
#' \item{`hetvar` is used for the construction of the F tests. Set to \code{1}
#' if users want to allow for the variance of the residuals to be different across segments. If \code{hetvar}=\code{0}, 
#' the variance of the residuals is assumed constant across segments and constructed from the full sample. 
#' \code{hetvar}=\code{1} when \code{robust} =\code{1}.}
#' \item{`hetomega` used in the construction of the confidence intervals for the break
#' dates. If \code{hetomega}=\code{0}, the long run covariance matrix of zu is
#' assumed identical across segments (the variance of the errors u if \code{robust}={0}).}
#' \item{`hetq` used in the construction of the confidence intervals for the break
#' dates. If \code{hetq}=\code{0}, the moment matrix of the data is assumed identical
#' across segments.}
#' }
#'
#' @details
#' Structural break models for `sbm` are specified symbolically similar to linear regression 
#' model specified in `lm`. A typical structural break model has the form `y ~ z | x` where 
#' response `y` is the (numeric) response vector and `x` and `z` terms is a series of terms which 
#' are classified in two distinct types. `x` terms have unchanged coefficients across regimes and `z` terms 
#' have varied coefficients across regimes. A terms specification of the form `z1 + z2 | x1 + x2`  indicates 
#' all the terms in `z1` with all the terms in `z2` with duplicates removed have coefficients varied 
#' across regimes, whereas all terms in `x1` with all terms in `x2` with duplicates removed have
#' coefficients fixed across regimes. If no `x` regressors are specified in the symbolic formula,
#' the model is pure structural break model, whereas if any `x` regressors are specified, the model is
#' partial structural break model. A formula has an implied intercept term that varies across regimes. 
#' If intercept is explicitly removed by using `y ~ z - 1` for pure structural model, or `y ~ z - 1 | x` for partial structural
#' change model, then the formula has an implied intercept that fixes across regimes. To remove intercept from
#' the model entirely, use `y ~ z - 1| x - 1`
#' 
#' @return `sbm` returns an object of class `sbm`. 
#' The functions summary is used to obtain and print a summary. The generic accessor
#' functions `` extract useful features of the value returned by `sbm`. An object of
#' class `sbm` is a list containing at least the following components:
#' 
#' @seealso [dotest()], [doseqtests()], [doorder()], [dosequa()], and [dofix()]
#' which are functions called by [mdl()].
#'  
#' @export
#' @examples
#' model1 <- sbm(inf~inflag+inffut|ygap,data=nkpc)
#' model2 <- sbm(rate~1,data=real)

sbm <- function(formula, data, eps1=0.15, cov. = list(),...){
  
  default.Opts <- list('prewhit'=1,'robust'=1,
                       'hetdat'=1,'hetvar'=1,'hetq'=1,'hetomega'=1)
  default.maxM = c(2,3,5,8,9)
  default.eps1 = c(0.05,0.10,0.15,0.20,0.25)
  
  #save all the values passed in arguments for later functions
  argAll <- c(as.list(environment()), list(...))
  
  #collect all arguments passed to covariance matrix specification
  #Opts = c('prewhit', 'robust', 'hetdat' , 'hetvar','hetq','hetomega')
  
  #copy down argument for covariance matrix specification
  default.Opts[names(cov.)] <- cov.
  argMatched <- default.Opts
  
  
  #extract structural break model from supplied formula
  mf.call <- match.call(expand.dots = FALSE)
  mf.args.matched <- match(c("formula", "data", "subset", "na.action"), names(mf.call), 0L)
  mf.call <- mf.call[c(1L, mf.args.matched)]
  #get model from symbolic expression/formula
  mf.f <- check_formula(formula)
  
  mf.call[[2]] <- mf.f #substitute reformatted formula
  mf.call$drop.unused.levels <- TRUE
  mf.call[[1]] <- quote(stats::model.frame) #set call to model.frame to set up data frame
  mf <- eval(mf.call, parent.frame()) #evaluate model.frame() from saved call mf.call, and saved data frame to mf
  
  #extract frames as y,z and x matrices
  y <- as.matrix(mf[,1,drop=FALSE])
  bigT = dim(y)[1]
  const = matrix(1,bigT,1)
  
  if (attr(mf.f,'z-intercept') == 1){
    zz <-mf[,attr(mf.f,'z-regs'),drop=FALSE]
    z <- cbind(const,zz)
    colnames(z) <- c('Intercept',attr(mf.f,'z-regs'))
    z <- as.matrix(z)
  } else { 
    z<-as.matrix(mf[,attr(mf.f,'z-regs'),drop=FALSE])}

  if (is.na(attr(mf.f,'x-regs'))) 
    {x = matrix(NA,bigT,0)}
  else{
    x <- as.matrix(mf[,attr(mf.f,'x-regs'),drop=FALSE])}
  if (attr(mf.f,'x-intercept') == 1){
    x <- cbind(const,x)
    colnames(x) <- c('Intercept', attr(mf.f,'x-regs'))
    x <- as.matrix(x)
  }

  #compute other parameters to fit in arguments of aux functions
  p = dim(x)[2]
  q = dim(z)[2]
  
  #check user specifications
  h = bigT*eps1
  idx.eps <- which(default.eps1 == eps1) 
  if (identical(integer(0),idx.eps)) {stop(paste('Invalid trimming level eps1 = ',eps1,'. Specify only 0.05, 0.10, 0.15, 0.20 or 0.25',sep=''))}
  m = default.maxM[idx.eps]
  conditions = check_trimming(bigT,eps1,m,h,p,q)
  
  h=conditions$h
  eps1=conditions$eps1
  m=conditions$m

  #estimate the break date with given number of break
  if (p == 0){
    est_mdl <- dating(as.matrix(y),as.matrix(z),h,m,q,bigT)
  }
  else{
    est_mdl <- nldat(y,z,x,h,m,p,q,bigT)
  }
  est_m <- m
  brdate <- as.matrix(est_mdl$datevec)
  colnames(brdate) <-paste('m=',c(1:5), sep='')
  gSSR <- cbind(c(1:m),est_mdl$glb)
  colnames(gSSR) <- c('numBreaks','SSR')
  vecSSR <- est_mdl$bigvec
  #estimate the sbm model parameters with estimated break date (for inference)
  #and covariance specification in params
  #construct a call to estim using arguments matched, and store this as a new class inside
  sbm_estCALL <- as.call(c(list(quote(estim),'m'=m,'q'=q,'z'=z,'y'=y,'b'=as.matrix(brdate[,m]),'x'=x,'p'=p),argMatched))
  #print(as.list(sbm_estCALL))
  sbm_est <- eval(sbm_estCALL)
  out = structure(list(),class='sbm')

  #group data into one separate field model.sbm
  model.sbm <- cbind(y,z,x)
  #attr(model.sbm, 'dep.var') <- colnames(y)
  attr(model.sbm, 'xregs') <- if(p == 0) NA else colnames(x)
  attr(model.sbm, 'zregs') <- colnames(z)
  #other attributes of the model
  attr(model.sbm, 'formula') <- formula
  attr(model.sbm, 'MAXbreaks') <- m
  attr(model.sbm, 'min.segment') <- h
  attr(model.sbm, 'trim.level') <- eps1
  

  # break date estimates and CI stored as bound (for inference); and also for further analysis
  breakdate <- c()
  breakdate$date <- brdate
  breakdate$bound <- sbm_est$CI
  attr(breakdate, 'numbreaks') <- est_m #this denotes differently
  #to reflect it is estimated by IC/test (for later more complex class)
  #could be more way of estimating confidence interval for break date later (for example bootstrap)
  attr(breakdate, 'procedure') <- 'Asymptotic'
  # regressors estimates and CI stored as stdev (for inference)
  
  #add backs all above fields
  out$model <- model.sbm
  out$breakdate <- breakdate
  out$maxM <- sbm_est  #mb class with respect to maximum number of breaks estimated specified by users via m

  #post-regression
  out$SSR <- gSSR     #global minimum SSR 
  out$vecSSR <-vecSSR #triangular SSR
  
  #keep tracks of options used
  out$opts <- as.list(argMatched)
  class(out) <- 'sbm'
  
  #select breaks based on IC
  ic.temp <- ic(out)
  out$IC <- ic.temp$IC
  out$SQ <- seqm(out)
  SQ.CALL <- as.call(c(list(quote(estim),'m'=out$SQ$SQ,'q'=q,'z'=z,'y'=y,'b'=as.matrix(out$SQ$date),'x'=x,'p'=p),argMatched))
  out$SQ$model <- eval(SQ.CALL)
  
  return(out)
}


#Lists of methods to do
#1) summary
#2) print (format as suggested)
#3) plot*** (important)
#4) coef (must be able to split out estimated coefs on z and x (optional))
#5) brdate (get the estimated break date out, in vector and other form for display)
#6) model.frame (need to rewrite this), which extracting model frame from a formula/fit
# (Unnecessary, can reuse base R)
#7) resid (return already computed residuals of the model)


#' Custom print method for MyClass
#'
#' This method provides a custom printing behavior for objects of class MyClass.
#' It displays the attribute of the object.
#'
#' @param x An object of class "sbm"
#'
#' @export
print.sbm <- function(x,digits = 2,...){
  if (is.na(attr(x$model,'xregs') ) ){p=0}else{p <- length(attr(x$model,'xregs'))}
  q <- length(attr(x$model,'zregs'))
  m <- attr(x$model,'MAXbreaks')
  if (p == 0) type <- 'Pure' else type <- 'Partial'
  cat('\n')
  cat(paste(type,'structural break model with maximum', m ,'breaks\n'))
  cat(paste('Model specification: ',format(attr(x$model,'formula')),'\n'))
  msel <- x$SQ$SQ 
  cat(paste('\n',msel,' breaks selected by sequential procedure',sep=''))
  if (msel> 0) {
    print(x$SQ$model,digits)
  } 
  #print(x$maxM,digits)
  invisible(x)
}

#' @export
summary.sbm <- function(x, digits = 2, ... ){
  if (is.na(attr(x$model,'xregs') ) ){p=0}else{p <- length(attr(x$model,'xregs'))}
  q <- length(attr(x$model,'zregs'))
  m <- attr(x$model,'MAXbreaks')
  if (p == 0) type <- 'Pure' else type <- 'Partial'
  cat('\n')
  cat(paste(type,'structural break model with maximum', m ,'breaks\n'))
  cat(paste('Model specification: ',format(attr(x$model,'formula')),'\n'))
  
  if (!is.null(x$IC)){
    cat('\nNumber of breaks selected by information criteria:\n')
    cat(paste(names(x$IC),x$IC,sep='='))
    cat('\n')
  }else{
    cat('\nNo break selection procedure via information criterion.\n')
  }
  
  if (!is.null(x$SQ)){
    msel <- a$SQ$SQ 
    cat(paste('\nNumber of breaks selected by sequential procedure:',msel,'\n'))
  }else{
    cat('No break selection via sequential procedure.')
  }
  cat('For diagnostic tests, run suptests() and seqtests()')
  
  #options for covariance
  #summarize SSRs
  #summarize break selection
  #summarize tests
  invisible(x)
}

#' @export
brdate.sbm <- function(x,...) x$breakdate$date

#' @export
coef.sbm <- function(x,...){ 
  out <- c()
  numbreaks = attr(x$breakdate, 'estbreaks')
  p = attr(x$coefficients, 'numx')
  q = attr(x$coefficients, 'numz')
  cname.z = c()
  coef.z = matrix(0L,q,numbreaks+1)
  rname.z = attr(x$coefficients, 'zregs')
  for (i in 1:(numbreaks+1)){
    cname.z = cbind(cname.z,paste('Regime',i))
  }
  for (j in 1:q){
    coef.z[j,] = x$coefficients$coef.z[((j-1)*(numbreaks+1)+1):(j*(numbreaks+1))]}
  colnames(coef.z) <- cname.z
  rownames(coef.z) <- rname.z
  if (p == 0){
    coef.x = NA
    colnames(coef.x) = NA
  }else{
  coef.x = x$coefficients$coef.x
  colnames(coef.x) <- attr(x$coefficients, 'xregs')
  }
  out$coefficients.z <- coef.z
  out$coefficients.x <- coef.x
  out
}


#' @export 
getCovOpts <- function(x,...){
  #default list of options for covariance matrix
  opts = x$opts
  cat(gettextf('robust = %s; hetdat = %s; hetvar = %s; hetomega = %s; hetq = %s',
               opts$robust,opts$hetdat,opts$hetvar,opts$hetomega,opts$hetq))
}

#' @export
residuals.sbm <- function(x,...) x$residuals

#' @export
fitted.sbm <- function(x,...) x$fitted.values





#' check object class
#' @noRd

is.sbm <- function(x) inherits(x, 'sbm')

##generics declaration
#summary <- function(x, ...){
#  UseMethod('summary',x)
#}
# print <- function(x, ...){
#  UseMethod('print',x)
# }
# print.default <- base::print
# 
# plot <- function(x, ...){
#   UseMethod('plot',x)
# }
# coef <- function(x, ...){
#   UseMethod('coef',x)
# }
# residuals <- function(x,...){
#   UseMethod('residuals',x)
# }
# fitted <- function(x,...){
#   UseMethod('fitted',x)
# }
# brdate <- function(x,...){
#   UseMethod('brdate',x)
# }
# getOption <- function(x,...){
#   UseMethod('getOption',x)
# }

#extra internal functions for analysis of a sbm class
# sbm.select <- function (x,maxm=5){
#   
# }


### Extra explanations
### NOTE ###
#create new S3 class (building block of everything):
#S3 class to capture the structural break model, will have the following
#properties:
#   numbr: number of breaks (not maximum, but a specify number of breaks); 
#other S3 class can inherit a different meaning of numbr later 
#which is different from the maximum number of break you want
# numbr is numeric
#   mform: special type of class formula, which can understand the formula syntax
# dedicated to mbreaks package. (More detailed on mformula.R)
###


#The model should be:
#Input is a formula or a data
#initialization of class sbm (only care about estimation of structural break model given
# number of breaks, not concerning with break selection. Break selection model should be
# subclass inherits many of the methods below from this block model/class)
# this class will only use specific trimming. For estimation with specific h, use doglob()
# or create a sub function to estimate specific h and m
# this return a totally new sbm class for printing as well

ic <- function(s, IC = c('BIC','LWZ','KT')){
  available.IC <- c('BIC','LWZ','KT')
  matched.IC <- match(available.IC,IC,0)
  IC <- IC[matched.IC]
  if (is.null(ic)){stop('Invalid information criteria. Please use either BIC, LWZ and/or KT')}
  if(is.na(attr(s$model,'xregs'))) p <- 0 else p <- length(attr(s$model,'xregs'))
  q <- length(attr(s$model,'zregs'))
  m <- attr(s$model,'MAXbreaks')
  z_name <- attr(s$model,'zregs')
  x_name <- attr(s$model,'xregs')
  y <- as.matrix(s$model[,1])
  z <- as.matrix(s$model[,z_name])
  x <- as.matrix(s$model[,x_name])
  bigT = dim(s$model)[1]
  if (p == 0){zz = z} else{zz = cbind(z,x)}
  glb = as.matrix(s$SSR[,'SSR'])
  ssr0 = nssr(y,zz)
  delta0 = 0.1 #optimal parameters in LWZ paper
  c0 = 0.299
  glob= matrix(0L, nrow = m+1, ncol=1)
  glob[1,1] = ssr0
  glob[seq(2,m+1),1] = glb
  datevec = s$breakdate$date
  bic = matrix(0L,nrow = m+1, ncol = 1)
  lwz = matrix(0L,nrow = m+1, ncol = 1)
  kt  = matrix(0L,nrow = m+1, ncol = 1)
  for (i in seq(1,m+1)){
    #BIC criterion
    bic [i,1] = log(glob[i,1]/bigT) + log(bigT)*(i-1)*(q+1)/bigT
    #LWZ criterion
    lwz[i,1] = log(glob[i,1]/(bigT-i*q-i+1)) +
    ((i-1)*(q+1)*c0*(log(bigT))^(2+delta0))/bigT
    #Kurozumi and Tuvaandori (2011)
    if (i==1){bd=c(0,bigT)}
    else{bd=c(0,datevec[1:i-1,i-1],bigT)}
    for (l in seq(1,i)){
      segy   = y[seq(bd[l]+1,bd[l+1],1),,drop=FALSE]
      segz   = z[seq(bd[l]+1,bd[l+1],1),,drop=FALSE]
      segres = segy-segz%*%solve(t(segz)%*%segz)%*%t(segz)%*%segy
      dt     = bd[l+1]-bd[l]
      kt[i,1]= kt[i,1]+(dt*log(t(segres)%*%segres/dt)+q*log(dt))}
      kt[i,1]    = kt[i,1]+2*i*log(bigT)}
    
  mBIC = which.min(bic) - 1
  mLWZ = which.min(lwz) - 1
  mKT = which.min(kt) - 1
  bic <-  cbind(seq(0,m,1),bic)
  lwz <-  cbind(seq(0,m,1),lwz)
  kt  <-  cbind(seq(0,m,1),kt)
  colnames(bic) <- c('m','penalized SSR')
  colnames(lwz) <- c('m','penalized SSR')
  colnames(kt) <- c('m','penalized SSR')
  mSEL = c('BIC' = mBIC, 'LWZ' = mLWZ, 'KT' = mKT)
  penSSR = list('BIC' = bic, 'LWZ' = lwz , 'KT' = kt)
  mSEL.out <- list()
  mSEL.out$IC <- mSEL[IC]
  mSEL.out$mSSR <- penSSR[IC]
  mSEL.out
}

#select breaks by sequential approach
#' @export
seqm <- function(s, siglev = 0.05){
  #needed options for sequential approach
  reqOpts <- c('robust','prewhit','hetdat','hetvar')
  opts <- s$opts[reqOpts]
  
  signifAll <- c(0.10,0.05,0.025,0.01)
  signif <- match(siglev,signifAll)
  if(is.na(signif)) stop('No significant level found. Please use only 0.10, 0.05, 0.025 or 0.01')
  
  if(is.na(attr(s$model,'xregs'))) p <- 0 else p <- length(attr(s$model,'xregs'))
  q <- length(attr(s$model,'zregs'))
  m <- attr(s$model,'MAXbreaks')
  z_name <- attr(s$model,'zregs')
  x_name <- attr(s$model,'xregs')
  y <- as.matrix(s$model[,1])
  z <- as.matrix(s$model[,z_name])
  x <- as.matrix(s$model[,x_name])
  bigT <- length(y)
  sequaCALL <- as.call(c(list(quote(sequa),
          'm'=m,'q'=q,'z'=z,'y'=y,'x'=x,'p'=p,'bigT'=bigT,'signif'=signif,'eps1'=attr(s$model,'trim.level'),'h'=attr(s$model,'min.segment')),opts))
  seq <- eval(sequaCALL)
  out <- list('SQ' = seq$nbreak, 'date' = as.matrix(seq$dv0,ncol=1))
  out
}


#select breaks by repartition approach
repart <- function(s, siglev = 0.05){
  
  if(is.na(attr(s$model,'xregs'))) p <- 0 else p <- length(attr(s$model,'xregs'))
  q <- length(attr(s$model,'zregs'))
  m <- attr(s$model,'MAXbreaks')
  z_name <- attr(s$model,'zregs')
  x_name <- attr(s$model,'xregs')
  y <- as.matrix(s$model[,1])
  z <- as.matrix(s$model[,z_name])
  x <- as.matrix(s$model[,x_name])
  h <- attr(s$model,'min.segment')
  bigT <- length(y)
  
  temp = seqm(s,siglev)
  if (temp$SQ == 0){
    stop('There are no breaks selected by sequential procedure and the repartition procedure is stopped')
  }
  else {
    repartda = preparti(y,z,temp$SQ,temp$date,h,x,p)
    reparv <- list('SQ' = temp$SQ, 'date' = repartda)
    reparv
  }
  
}

# all supF tests
suptests <- function(s){
  siglev=matrix(c(10,5,2.5,1),4,1)
  #needed options for sequential and sup F test
  reqOpts <- c('robust','prewhit','hetdat','hetvar')
  opts <- s$opts[reqOpts]
  
  if(is.na(attr(s$model,'xregs'))) p <- 0 else p <- length(attr(s$model,'xregs'))
  q <- length(attr(s$model,'zregs'))
  m <- attr(s$model,'MAXbreaks')
  z_name <- attr(s$model,'zregs')
  x_name <- attr(s$model,'xregs')
  y <- as.matrix(s$model[,1])
  z <- as.matrix(s$model[,z_name])
  x <- as.matrix(s$model[,x_name])
  h <- attr(s$model,'min.segment')
  eps1 <- attr(s$model,'trim.level')
  bigT <- length(y)
  
  ftest = matrix(0L, nrow = m, ncol = 1)
  wftest = matrix(0L, nrow = m, ncol = 1)
  for (i in 1:m){
    ftestCALL <- as.call(c(list(quote(pftest),
                               'y'=y,'q'=q,'i'=i,'x'=x,'p'=p,'bigT'=bigT,'datevec'=s$breakdate$date,'z'=z),opts))
    ftest[i,1] <- eval(ftestCALL)
  }
  
  cv_supF = matrix(0L,4,m)
  for (sl in 1:4){
    #critical values for supF test
    cv = getcv1(sl,eps1)
    #pad missing values with NA
    if(dim(cv)[2]<m){
      tmp_cv = dim(cv)[2]
      cv_supF[sl,1:tmp_cv] = cv[q,1:tmp_cv,drop=FALSE]
      cv_supF[sl,tmp_cv+1:m] = matrix(NA,1,m - tmp_cv)
    }else{
      cv_supF[sl,] = cv[q,1:m,drop=FALSE]}
  }
  
  #procedure for Dmax and UDmax test
  cv_Dmax = matrix(0L,4,1)
  for (sl in 1:4) {
    #critical values for Dmax test
    cvm = getdmax(sl,eps1)
    cv_Dmax[sl,1] = cvm[q,1]
  }
  
  for (sl in 1:4){
    #computation of WDmax test
    cv = getcv1(sl,eps1)
    for( i in 1:m){
      wftest[i,1] = cv[q,1] * ftest[i,1] / cv[q,i]
    }
  }
  rownames(cv_supF) = siglev
  rownames(cv_Dmax) = siglev
  if (m > 5){
    UDmax = max(ftest[1:5,])}
  else{
    UDmax = max(ftest)
  }
  
  out = list('ftest' = ftest, 'cv_supF' = cv_supF,
             'cv_Dmax' = cv_Dmax, 'UDmax' = UDmax)
  out$mbreak = m
  class(out) = 'sbtests'
  out = compile_sbtests(out)
  out
}

seqtests <- function(s){
  siglev=matrix(c(10,5,2.5,1),4,1)
  
  #needed options for sequential and sup F test
  reqOpts <- c('robust','prewhit','hetdat','hetvar')
  opts <- s$opts[reqOpts]
  
  if(is.na(attr(s$model,'xregs'))) p <- 0 else p <- length(attr(s$model,'xregs'))
  q <- length(attr(s$model,'zregs'))
  m <- attr(s$model,'MAXbreaks')
  z_name <- attr(s$model,'zregs')
  x_name <- attr(s$model,'xregs')
  y <- as.matrix(s$model[,1])
  z <- as.matrix(s$model[,z_name])
  x <- as.matrix(s$model[,x_name])
  h <- attr(s$model,'min.segment')
  eps1 <- attr(s$model,'trim.level')
  bigT <- length(y)
  
  supfl = matrix (0L,m,1)
  ndat = matrix (0L,m,1)
  pftestCALL <- as.call(c(list(quote(pftest),
                               'bigT'=bigT,'y'=y,'q'=q,'i'=1,'x'=x,'p'=p,'datevec'=s$breakdate$date[1:1,1,drop=FALSE],'z'=z),opts))
  supfl[1,1] = eval(pftestCALL)
  if (m > 1){
  for (i in seq(1,m-1,1)){
    seqFtestsCALL <- as.call(c(list(quote(spflp1),
                                    'y'=y,'bigvec'=s$vecSSR,'q'=q,'nseg'=i+1,'x'=x,'p'=p,'dt'=s$breakdate$date[1:i,i,drop=FALSE],'z'=z,'h'=h),opts))
    out1 <- eval(seqFtestsCALL)
    supfl[i+1,1] = out1$maxf}
  } 
  cv_supFl = matrix(0L,4,m)
  
  for (c in 1:4){
    cv = getcv2(c,eps1)
    cv_supFl[c,] = cv[q,1:m,drop=FALSE]
  }
  
  rownames(cv_supFl) = siglev
  
  out = list('supfl' = supfl, 'cv' = cv_supFl)
  out$mbreak = m
  class(out) = 'seqtests'
  out = compile_seqtests(out)
  return(out)
}
# additional estimation function for specific h and m.
estimate <- function(s,mnew,h){
  
  q <- length(attr(s$model,'zregs'))
  z_name <- attr(s$model,'zregs')
  z <- s$model[,z_name,drop=FALSE]
  y <- s$model[,1,drop=FALSE]
  bigT <- length(y)
  if(is.na(attr(s$model,'xregs'))) {
    p <- 0 
    x <- matrix(NA,nrow=bigT,ncol=0)
  }
  else {
    p <- length(attr(s$model,'xregs'))
    x_name <- attr(s$model,'xregs')
    x <- s$model[,x_name,drop=FALSE]
  }
  
  #check user specifications
  h.checked <- check_h_estim(bigT,mnew,h,p+q)
  
  
  #estimate the break date with given number of break
  if (p == 0){
    est_mdl <- dating(y,z,h.checked,mnew,q,bigT)
  }
  else{
    est_mdl <- nldat(y,z,x,h.checked,mnew,p,q,bigT)
  }
  
  #construct a call to estim using arguments matched, and store this as a new class inside
  sbm_estCALL <- as.call(c(list(quote(estim),'m'=mnew,'q'=q,'z'=z,'y'=y,'b'=as.matrix(est_mdl$datevec[,mnew]),'x'=x,'p'=p),s$opts))
  out <- eval(sbm_estCALL)
  out
  
}

