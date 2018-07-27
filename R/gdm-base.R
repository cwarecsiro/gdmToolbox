#' @title gdm model fit with base R glm
#' @rdname gdm.base
#' @name gdm.base
#' @description function for fitting a gdm model. The function differs to the \code{\link{gdm}} function
#' only in that it uses the base R glm function to construct and fit the model rather than using the existing 
#' C code. The advantages of fitting a gdm in this way may be limited for some use cases, and will likely incur a 
#' performance (speed) penalty. For other use cases, however, fitting a model using the glm framework provides a flexible
#' way of model fitting not possible using the existing C code (in which model options were hard-coded). As with 
#' \code{\link[stats]{glm.fit}}, gdm.base allows the specification of different link functions, 
#' families, opimizers, likelihood functions, and spline functions (possibly... to come...), in addition to 
#' providing a potentially more useful (or at least familiar) summary output.
#' @param data (data.frame) gdm site pair table of the format returned by \code{\link{formatsitepair}}.
#' @param link (character or function) link function. default is 'negexp' (e.g., Ferrier etal., 2007), can call other
#' link functions \code{\link[stats]{make.link}} or custom link function.
#' @param optim.meth (character or function) optimisation method (default nnls).
#' @param geo (boolean) Whether or not geographic distance should be included in the model 
#' @param geo.method (character) measure to use to calculate geographic distance. Currently supports 'euclidean' or 'vincenty' (requires SDMTools)
#' @param splines (vector) number of splines to transform model predictors. If splines has length of one, then the same number of splines is applied to each predictor.
#' Else, the length of splines must equal the number of predictors (as per the splines parameter in \code{\link{gdm}}).
#' @param quantiles (vector) as per knots parameter in \code{\link{gdm}}.
#' @param control model options to pass to glm.control(). Default glm.control().
#' @param output (character) output format. Options are one of 'glm' (default) or 'gdm.' If 'glm' model summary is that returned by the \code{\link{glm}}.
#' If 'trad' model summary is in the format returned by \code{\link{gdm}}.
#' @return gdm model fit as either glm object (default) of if output = 'gdm', returns 'traditional' output object. 
#' @note The present implementation differs in some key ways to the \code{\link{gdm}} function (e.g. different coefficients). 
#' The reason for this is not yet clear, but needs understanding before releasing the function into the wild... 
#' Also haven't included a weights arg into the model fit (yet).
#' Also also, documentation for the supporting functions is still to come...
#' Also also also, the functions here are probably best split over a few .R file, but are bundled for now. 
#' 
#' @importFrom nnls nnls
#' @importFrom SDMTools distance
#' 
#' @author Andrew Hoskins, Skipton Woolley, Chris Ware
#' 
#' @examples
#' ## create input table using gdm example  
#' load(system.file("./data/gdm.RData", package = 'gdm'))
#' sppData <- gdmExpData[, c(1,2,13,14)]
#' envTab <- gdmExpData[, c(2:ncol(gdmExpData))]
#' testData1a = reshape2::dcast(sppData, site~species)
#' coords <- unique(sppData[, 2:ncol(sppData)])
#' testData1b <- merge(testData1a, coords, by="site"))
#' exFormat1a <- suppressWarnings(formatsitepair(testData1a, 1, 
#'                                                siteColumn="site", 
#'                                                XColumn="Long", 
#'                                                YColumn="Lat",
#'                                                predData=envTab))
#' ## rm categorical (?) variable
#' data = exFormat1a[, -c(grep('bio18', names(exFormat1a)))]
#' ## fit the model
#' fitgdm = gdm.base(data, geo = TRUE, geo.method = 'vincenty', output = 'glm')
#' 
#' @export

gdm.base = function(data,
                    link = 'negexp',
                    optim.method = 'nnls.fit',
                    geo = FALSE,
                    geo.method = 'euclidean',
                    splines = 3,
                    quantiles = NULL,
                    control = glm.control(), 
                    output = 'glm'){ 
  
  ## checks
  if (is.character(optim.method))
    optim.method <- get(optim.method, mode = "function", envir = parent.frame())
  
  ## capture name for trad output
  dataname = deparse(substitute(data))

  ## do geo 
  if (geo){ 
    gdist = calc_geo(data[,3], data[,4], data[,5], data[,6], geo.method)
  }
  
  ## do splines 
  if (ncol(data) > 6) {
    splinedData = splineData(data[, -c(1:6)], splines, quantiles)  
    if(geo){
      gdist_splined  = splineData(gdist, splines, quantiles)  
      splinedData = list(data = cbind(gdist_splined$data, 
                                      splinedData$data),
                         quantiles = c(gdist_splined$quantiles, 
                                     splinedData$quantiles))
    }
  
  } else {
      ## assume here that if ncol(data) == 6, 
      ## the intention is to fit a distance only model,
      ## as per the current gdm implementation
    if(geo){
      gdist_splined  = splineData(gdist, splines, quantiles)  
      splinedData = list(data = gdist_splined$data,
                         quantiles = gdist_splined$quantiles)
      
    } else {
      stop(paste0('Input data table should have eight or more columns if',
                   'geo = FALSE, and be in the format returned by', 
                   'gdm::formatsitepair'))
    }
      
  }
  
  ## build data set
  data = data.frame(distance = data[,1], splinedData$data)  
  
  ## model config 
  formula = paste(colnames(data)[-1], collapse = "+")
  formula = as.formula(paste0(names(data)[1], "~", formula))
  
  ## fit model
  fit = tryCatch({
    ## The need to incldue this error handling should probably be 
    ## removed, but presently nnls struggles to fit a model in the case 
    ## where only geographic coordinates are supplied as predictors, 
    ## at least in the cases i've tried (i.e. distance only model 
    ## where data[, 1:6] is supplied as the input table) - maybe because
    ## it results in many 0 values once splined?
    ## Anayway, that's all this section deals with, over and above fitting 
    ## a model.
    if(link == 'negexp'){
      glm(formula, family = binomial(link = negexp()), 
          data = data, control = control,
          method = optim.method)
    
    } else {
      glm(formula, family = binomial(link = link), 
          data = data, control = control,
          method = optim.method)
      
    }
    }, error = function(e){e}
  )
  
  if(class(fit)[1] == 'simpleError'){
    if(length(grep('no valid set of coefficients has been found', 
                   as.character(fit))) > 0){
      ## smells like a distance only model with 0 spline coefficients - 
      ## fit a model to get starting params, and then try and re-fit model
      ## as requested. 
      warning(paste0('Model would not fit - ',
                     'trying again but supplying starting parameters'))
      start_coef <- glm(formula, family=binomial(link = negexp()), 
                        data = data)
      ## refit
      fit = tryCatch({
        if(link == 'negexp'){
          glm(formula, family = binomial(link = negexp()), 
              data = data, control = control,
              method = optim.method, start = start_coef$coefficients)
          
        } else {
          glm(formula, family = binomial(link = link), 
              data = data, control = control,
              method = optim.method, start = start_coef$coefficients)
          
        }
      }, error = function(e){e}
      )
      
      if(class(fit)[1] != 'glm'){
        stop('Could not fit a model')
      }
      
    }
  }
  
  ## do gdm output if requested
  if(output == 'gdm'){
    expl = 100 - fit$deviance / fit$null.deviance * 100
    predlist = lapply(strsplit(colnames(splinedData$data), '_'), '[')
    predlist = unique(unlist(lapply(predlist, function(x) x[-c(length(x))])))
    quant = splinedData$quantiles
    if(geo){
      predlist = c('geo', predlist)
      quant = c(gdist_splined$quantiles, quant)
    }
    
    if(length(splines) == 1) {
      splines = rep(splines, length(predlist))
    }
    
    gdmModOb <- structure(list(geo = geo, 
                               dataname = dataname,
                               sample = nrow(data), 
                               gdmdeviance = fit$deviance, 
                               nulldeviance = fit$null.deviance, 
                               explained = expl, 
                               intercept = unname(fit$coefficients[1]), 
                               predictors = predlist, 
                               coefficients = fit$coefficients, 
                               knots = quant, 
                               splines = splines, 
                               creationdate = date(), 
                               observed = fit$y, 
                               predicted = fit$fitted.values, 
                               ecological = fit$linear.predictors))
    
    class(gdmModOb) <- c("gdm", "list")
    return(gdmModOb)
    
  } else {
    return(fit)
  }

} 


calc_geo = function(x1, y1, x2, y2, method = 'euclidean'){
  if(method == 'euclidean'){
    dist = sqrt(((x2 - x1)^2 + (y2 - y1)^2))    
  } else if (method == 'vincenty'){
    dist = distance(y1, x1, y2, x2)
    dist = dist$distance
  } else {
    stop(paste0('Distance method not recognised - ',
                'only euclidean or vincenty allowed currently'))
  }
  return(dist)
}


ispline <- function(predVal, q1, q2, q3 ){
  
  outVal <- rep(NA,length(predVal))
  outVal[predVal <= q1] <- 0
  outVal[predVal >= q3] <- 1
  outVal[predVal > q1 & predVal <= q2] = 
    ( ( (predVal[predVal > q1 & predVal <= q2] - q1) * 
          (predVal[predVal > q1 & predVal <= q2] - q1) ) / 
        ( (q2-q1) * (q3-q1) ) )
  outVal[predVal > q2 & predVal < q3] = 
    ( 1.0 - ( ( (q3 - predVal[predVal > q2 & predVal < q3]) * 
                  (q3 - predVal[predVal > q2 & predVal < q3]) ) / 
                ( (q3-q2) * (q3-q1) ) ) )
  
  return(outVal)
  
}


splineData <- function(X, splines=NULL, quantiles=NULL){
  
  if(is.numeric(X)){ ## single vector case
    
    if(is.null(splines)){
      message("No splines specified. Using 3 splines (0%, 50%, 100% quantiles)")
      splines = 3
    }
    
    if(is.null(quantiles)){
      
      if(splines != 3){
        stop("Must specify quantile positions if all(splines) != 3")
      }
      
      quantiles = quantile(X, c(0,0.5,1))
      
    } else {
      
      quantiles = quantile(X, quantiles)
      
    }
    
    spl_geo = c()
    for(sp  in 1:splines){
      
      if(sp == 1){
        spl = ispline(X, quantiles[1], quantiles[1], quantiles[2])
              } 
      if(sp == splines){
        spl = ispline(X, quantiles[splines-1], quantiles[splines], 
                       quantiles[splines])
      }
      if(sp != 1 & sp != splines){
        spl = ispline(X, quantiles[sp-1], quantiles[sp], quantiles[sp+1])
      }
      
      spl_geo <- cbind(spl_geo, spl)
      
    }
    
    spl_geo = as.data.frame(spl_geo)
    names(spl_geo) = paste('geo_spl', 1:splines, sep = '_')
    
    return(list(data =  spl_geo, quantiles = quantiles))
    
  } else { ## table case
    
    X = as.matrix(X)
    
    ## fold X and create site vector
    nc <- ncol(X)
    nc2 <- nc/2
    if(nc %% 2 != 0){stop("X must be a matrix with even columns")}
    X1 <- X[,1:nc2]
    X2 <- X[,(nc2+1):nc]
    nms <- colnames(X1)
    ## assumes 's1. / s2. is the only site pair notation used...
    ## probably needs to be generalised
    nms = gsub('s1.', '', nms) 
    colnames(X2) <- nms
    ## site vector
    sv <- c(rep(1,nrow(X1)),rep(2,nrow(X2)))
    XX <- rbind(X1,X2)
    
    if(is.null(splines)){
      message("No splines specified. Using 3 splines (0%, 50%, 100% quantiles)")
      splines <- 	rep(3,ncol(XX))
    } else {
      if(length(splines) == 1){
        splines = rep(splines, ncol(XX))
      }
    }
    
    if(is.null(quantiles)){
      if(all(splines != 3)){
        stop("Must specify quantile positions if all(splines) != 3")
      }
      quantiles <- unlist(lapply(1:ncol(XX), function(x)
        quantile(XX[,x],c(0,0.5,1))))
    }
    
    if(length(quantiles) != sum(splines)){
      stop("Number of quantiles must equal number of splines")
    }
    
    ## spline data
    csp <- c(0,cumsum(splines))
    out.tab <- c()
    for(c in 1:ncol(XX)){
      # c = 1
      ns <- splines[c]
      predVal <- XX[,c]
      quan <- quantiles[(csp[c]+1):(csp[c]+ns)]
      
      for(sp  in 1:ns){
        
        if(sp == 1){
          spl = ispline(predVal,quan[1],quan[1],quan[2])
        }
        if(sp == ns){
          spl = ispline(predVal,quan[ns-1],quan[ns],quan[ns])
        }
        if(sp != 1 & sp != ns){
          spl = ispline(predVal,quan[sp-1],quan[sp],quan[sp+1])
        }
        out.tab <- cbind(out.tab, spl)
        
        if(anyNA(spl)){
          print(c);print(sp);break
        }
        
      }
    }
    
    NMS <- rep(nms,splines)
    SPNMS <- paste0("spl",unlist(lapply(splines,function(x) 1:x)))
    NMS <- paste(NMS, SPNMS, sep="_")
    colnames(out.tab) <- NMS
    
    XX1 <- out.tab[sv == 1,]
    XX2 <- out.tab[sv == 2,]
    
    xout = as.data.frame(abs(XX1 - XX2))
    return(list(data = xout, quantiles = quantiles))
  }
  
}


nnls.fit <- function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL, 
                      mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
                      control = list(), intercept = TRUE) {
  control <- do.call("glm.control", control)
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y)) 
    rownames(y)
  else names(y)
  conv <- FALSE
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  if (is.null(weights)) 
    weights <- rep.int(1, nobs)
  if (is.null(offset)) 
    offset <- rep.int(0, nobs)
  variance <- family$variance
  linkinv <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv)) 
    stop("'family' argument seems not to be a valid family object", 
         call. = FALSE)
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  unless.null <- function(x, if.null) if (is.null(x)) 
    if.null
  else x
  valideta <- unless.null(family$valideta, function(eta) TRUE)
  validmu <- unless.null(family$validmu, function(mu) TRUE)
  if (is.null(mustart)) {
    eval(family$initialize)
  }
  else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  if (EMPTY) {
    eta <- rep.int(0, nobs) + offset
    if (!valideta(eta)) 
      stop("invalid linear predictor values in empty model", 
           call. = FALSE)
    mu <- linkinv(eta)
    if (!validmu(mu)) 
      stop("invalid fitted means in empty model", call. = FALSE)
    dev <- sum(dev.resids(y, mu, weights))
    w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
    residuals <- (y - mu)/mu.eta(eta)
    good <- rep(TRUE, length(residuals))
    boundary <- conv <- TRUE
    coef <- numeric()
    iter <- 0L
  }
  else {
    coefold <- NULL
    eta <- if (!is.null(etastart)) 
      etastart
    else if (!is.null(start)) 
      if (length(start) != nvars) 
        stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", 
                      nvars, paste(deparse(xnames), collapse = ", ")), 
             domain = NA)
    else {
      coefold <- start
      offset + as.vector(if (NCOL(x) == 1L) 
        x * start
        else x %*% start)
    }
    else family$linkfun(mustart)
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta))) 
      stop("cannot find valid starting values: please specify some", 
           call. = FALSE)
    devold <- sum(dev.resids(y, mu, weights))
    boundary <- conv <- FALSE
    for (iter in 1L:control$maxit) {
      good <- weights > 0
      varmu <- variance(mu)[good]
      if (any(is.na(varmu))) 
        stop("NAs in V(mu)")
      if (any(varmu == 0)) 
        stop("0s in V(mu)")
      mu.eta.val <- mu.eta(eta)
      if (any(is.na(mu.eta.val[good]))) 
        stop("NAs in d(mu)/d(eta)")
      good <- (weights > 0) & (mu.eta.val != 0)
      if (all(!good)) {
        conv <- FALSE
        warning("no observations informative at iteration ", 
                iter)
        break
      }
      z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
      w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
      ngoodobs <- as.integer(nobs - sum(!good))
      
      ## USE NNLS instead of least squares
      require(nnls)
      fit <- nnnpls(x[good, , drop = FALSE] *  w, z * w,con=c(-1,rep(1,ncol(x)-1)))
      fit$coefficients <- fit$x
      ## calculate qr decomposition..
      QR <- qr(x[good, , drop = FALSE] *  w,tol=min(1e-07, control$epsilon/1000))
      fit$qr <- QR$qr
      fit$rank <- QR$rank
      fit$pivot <- QR$pivot
      fit$qraux <- QR$qraux	
      fit$effects <- fit$fitted
      ### END NNLS
      
      if (any(!is.finite(fit$coefficients))) {
        conv <- FALSE
        warning(gettextf("non-finite coefficients at iteration %d", 
                         iter), domain = NA)
        break
      }
      if (nobs < fit$rank) 
        stop(gettextf("X matrix has rank %d, but only %d observations", 
                      fit$rank, nobs), domain = NA)
      start[fit$pivot] <- fit$coefficients ## change here - remove $qr$
      eta <- drop(x %*% start)
      mu <- linkinv(eta <- eta + offset)
      dev <- sum(dev.resids(y, mu, weights))
      if (control$trace) 
        cat("Deviance =", dev, "Iterations -", iter, 
            "\n")
      boundary <- FALSE
      if (!is.finite(dev)) {
        if (is.null(coefold)) 
          stop("no valid set of coefficients has been found: please supply starting values", 
               call. = FALSE)
        warning("step size truncated due to divergence", 
                call. = FALSE)
        ii <- 1
        while (!is.finite(dev)) {
          if (ii > control$maxit) 
            stop("inner loop 1; cannot correct step size", 
                 call. = FALSE)
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
          dev <- sum(dev.resids(y, mu, weights))
        }
        boundary <- TRUE
        if (control$trace) 
          cat("Step halved: new deviance =", dev, "\n")
      }
      if (!(valideta(eta) && validmu(mu))) {
        if (is.null(coefold)) 
          stop("no valid set of coefficients has been found: please supply starting values", 
               call. = FALSE)
        warning("step size truncated: out of bounds", 
                call. = FALSE)
        ii <- 1
        while (!(valideta(eta) && validmu(mu))) {
          if (ii > control$maxit) 
            stop("inner loop 2; cannot correct step size", 
                 call. = FALSE)
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
        }
        boundary <- TRUE
        dev <- sum(dev.resids(y, mu, weights))
        if (control$trace) 
          cat("Step halved: new deviance =", dev, "\n")
      }
      if (((dev - devold)/(0.1 + abs(dev)) >= control$epsilon) & 
          (iter > 1)) {
        if (is.null(coefold)) 
          stop("no valid set of coefficients has been found: please supply starting values", 
               call. = FALSE)
        warning("step size truncated due to increasing deviance", 
                call. = FALSE)
        ii <- 1
        while ((dev - devold)/(0.1 + abs(dev)) > -control$epsilon) {
          if (ii > control$maxit) 
            break
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
          dev <- sum(dev.resids(y, mu, weights))
        }
        if (ii > control$maxit) 
          warning("inner loop 3; cannot correct step size")
        else if (control$trace) 
          cat("Step halved: new deviance =", dev, "\n")
      }
      if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
        conv <- TRUE
        coef <- start
        break
      }
      else {
        devold <- dev
        coef <- coefold <- start
      }
    }
    if (!conv) 
      warning("glm.fit2: algorithm did not converge. Try increasing the maximum iterations", 
              call. = FALSE)
    if (boundary) 
      warning("glm.fit2: algorithm stopped at boundary value", 
              call. = FALSE)
    eps <- 10 * .Machine$double.eps
    if (family$family == "binomial") {
      if (any(mu > 1 - eps) || any(mu < eps)) 
        warning("glm.fit2: fitted probabilities numerically 0 or 1 occurred", 
                call. = FALSE)
    }
    if (family$family == "poisson") {
      if (any(mu < eps)) 
        warning("glm.fit2: fitted rates numerically 0 occurred", 
                call. = FALSE)
    }
    if (fit$rank < nvars) 
      coef[fit$pivot][seq.int(fit$rank + 1, nvars)] <- NA ## change here - $qr$ removed
    xxnames <- xnames[fit$pivot] # change here - $qr$ removed
    residuals <- (y - mu)/mu.eta(eta)
    fit$qr <- as.matrix(fit$qr) ## change here $qr$ removed
    nr <- min(sum(good), nvars)
    if (nr < nvars) {
      Rmat <- diag(nvars)
      Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars] ## change here $qr$ removed
    }
    else Rmat <- fit$qr[1L:nvars, 1L:nvars] ## change here $qr$ removed
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    names(coef) <- xnames
    colnames(fit$qr) <- xxnames ## change here $qr$ removed
    dimnames(Rmat) <- list(xxnames, xxnames)
  }
  names(residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames
  wt <- rep.int(0, nobs)
  wt[good] <- w^2
  names(wt) <- ynames
  names(weights) <- ynames
  names(y) <- ynames
  if (!EMPTY) 
    names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("", 
                                                                sum(good) - fit$rank))
  wtdmu <- if (intercept) 
    sum(weights * y)/sum(weights)
  else linkinv(offset)
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)
  rank <- if (EMPTY) 
    0
  else fit$rank
  resdf <- n.ok - rank
  aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
  list(coefficients = coef, residuals = residuals, fitted.values = mu, 
       effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat, 
       rank = rank, qr = if (!EMPTY) structure(fit[c("qr", 
                                                     "rank", "qraux", "pivot", "tol")], class = "qr"), 
       family = family, linear.predictors = eta, deviance = dev, 
       aic = aic.model, null.deviance = nulldev, iter = iter, 
       weights = wt, prior.weights = weights, df.residual = resdf, 
       df.null = nulldf, y = y, converged = conv, boundary = boundary)  ## change in here $qr removed from fit$qr[c('qr....)] now fit[c(...)]
}


negexp <- function(){
  ## link
  linkfun <- function(mu) -log(1-mu)
  ## inverse link
  linkinv <- function(eta) 1-exp(-eta)
  ## derivative of inverse link wrt eta
  mu.eta <- function(eta) exp(-eta)
  valideta <- function(eta) all(is.finite(eta))
  link <- paste0("negexp")
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, name = link),
            class = "link-glm")
}

