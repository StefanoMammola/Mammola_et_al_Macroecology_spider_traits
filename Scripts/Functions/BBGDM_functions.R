#Functions needed for null models of BBGDM ----------------------------------------------------------------------------

## SCRIPT FOR BBGDM

###### LOAD GDM PACKAGE
library(gdm)

#### LOAD CUSTOM BBGDM FUNCTIONS FROM SKIP
############
#' @title Running a bbgdm within the gdm modelling framework.
#' @rdname bbgdm
#' @name bbgdm
#' @description creates a \code{bbgdm} model, comprising multiple
#'   gdms. \code{bbgdm} models are core of \code{bbgdm}.
#' @param \dots for \code{bbgdm()}: one or more \code{plot()}, \code{print()} and
#'   \code{predict()}: further arguments passed to or from other methods
#' @param gdmTab formula for bbgdm model
#' @param geo presence absence matrix, sp as columns sites as rows.
#' @param splines environmental or spatial covariates at each site.
#' @param knots
#' @param bootstraps
#' @param cores
#' 
#' @author Skipton Woolley
#' 
#' @importFrom surveillance plapply
#' @importFrom plyr ldply
#' 
#' @examples
#' library(gdm)
#' load(system.file("./data/gdm.RData", package="gdm"))
#' sppTab <- gdmExpData[, c("species", "site", "Lat", "Long")]
#' envTab <- gdmExpData[, c(2:ncol(gdmExpData))]
#' gdmTab <- formatsitepair(sppTab, bioFormat=2, XColumn="Long", YColumn="Lat",
#'                          sppColumn="species", siteColumn="site", predData=envTab)
#' fm100 <- bbgdm(gdmTab, bootstraps=100, ncores=3)

bbgdm <- function(data, geo=FALSE, splines=NULL, knots=NULL, bootstraps=10, ncores=1){
  
  ## generate the number of sites for bayesian bootstrap.
  sitecols <- c('s1.xCoord','s1.yCoord','s2.xCoord','s2.yCoord')
  sitedat <- rbind(as.matrix(data[,sitecols[1:2]]),as.matrix(data[,sitecols[3:4]]))
  nsites <- nrow(sitedat[!duplicated(sitedat[,1:2]),])
  
  mods <- surveillance::plapply(1:bootstraps, bb_apply, nsites, data, geo, splines, knots, .parallel = ncores)
  
  median.intercept <- apply(plyr::ldply(mods, function(x) c(x$intercept)),2,median,na.rm=TRUE)
  quantiles.intercept <- apply(plyr::ldply(mods, function(x) c(x$intercept)),2, function(x) quantile(x,c(.05,.25,.5,.75, .95),na.rm=TRUE))  
  
  median.coefs <- apply(plyr::ldply(mods, function(x) c(x$coefficients)),2,median,na.rm=TRUE)
  quantiles.coefs <- apply(plyr::ldply(mods, function(x) c(x$coefficients)),2, function(x) quantile(x,c(.05,.25,.5,.75, .95),na.rm=TRUE))
  results <- list(gdms=mods,
                  median.intercept=median.intercept,
                  quantiles.intercept=quantiles.intercept,
                  median.coefs=median.coefs,
                  quantiles.coefs=quantiles.coefs)
  class(results) <- c("bbgdm", "list")
  return(results)
}

bb_apply <- function(x,nsites,data,geo,splines,knots){
  
  # creates the weights for bayes boot
  w <- gtools::rdirichlet(nsites, rep(1/nsites,nsites))
  wij <- w%*%t(w)
  wij <- wij[upper.tri(wij)]
  tmp <- data
  tmp$weights <- wij
  x <- gdm::gdm(tmp, geo=geo, splines=splines, knots=knots)
  return(x)
}

#' @title Running a bbgdm within the gdm modelling framework.
#' @rdname bbgdm
#' @name bbgdm.predict
#' @description creates a \code{bbgdm} model, comprising multiple
#'   gdms. \code{bbgdm} models are core of \code{bbgdm}.
#' @param mod bbgdm model object
#' @param newobs new sitepair table observations, can calculate this using the gdm functions.
#' @param bbgdm is
#' 
#' @author Skipton Woolley
#' 
#' @importFrom surveillance plapply
#' @importFrom plyr ldply
#' 
#' @examples
#' preds <- predict.bbgdm(fm100,gdmTab)
#' ## sanity check, should be pretty correlated
#' plot(gdmTab$distance,preds)


predict.bbgdm <- function(mod, newobs){
  
  predicted <- rep(0,times=nrow(newobs))
  object <- mod$gdms[[1]]
  knots <- mod$gdms[[1]]$knots
  splineNo <- mod$gdms[[1]]$splines
  intercept <- mod$median.intercept 
  coefs <- mod$median.coefs
  
  pred <- .C( "GDM_PredictFromTable",
              as.matrix(newobs),
              as.integer(FALSE),
              as.integer(length(object$predictors)), 
              as.integer(nrow(newobs)), 
              as.double(knots),
              as.integer(splineNo),
              as.double(c(intercept,coefs)),
              preddata = as.double(predicted),
              PACKAGE = "gdm")
  
  pred$preddata
  
}


#' @rdname bbgdm
#'
#' @method plot bbgdm
#' @param add_coefs if TRUE plot coefficent estimates on prediction plot
#' @param plot_derivate if TRUE plot derivated (numberically estimated) of spline predictions
#'
#' @export
#'
#' @examples
#' #plot bbgdm fit
#'  \dontrun{ 
#' x11(width = 7,height = 7)
#' par(mfrow=c(3,3),oma = c(5,4,0,0) + 0.1, mar = c(0.8,4,4,4) + 0.1)
#' plot.bbgdm(fm100)
#' 
#' x11(width = 7,height = 7)
#' par(mfrow=c(4,3),oma = c(5,4,0,0) + 0.1, mar = c(0.8,4,4,4) + 0.1)
#' plot.bbgdm(fm100,add_coefs = TRUE)
#' 
#' x11(width = 12,height = 7)
#' par(mfrow=c(3,4),oma = c(5,4,0,0) + 0.1, mar = c(0.8,4,4,4) + 0.1)
#' plot.bbgdm(mods,plot_derivate = TRUE)}

plot.bbgdm <- function(mods,add_coefs=FALSE,plot_derivate=FALSE, ...){

  if(plot_derivate==FALSE){
    PSAMPLE <- 200
    preddata <- rep(0,times=PSAMPLE)
    preds <- length(gdms[[1]]$predictors)
    predmax <- 0
    splineindex <- 1
    numsplines <- 3
    
    for(i in 1:preds){  
      i=1
      predplotdat <- lapply(gdms,function(x).C( "GetPredictorPlotData", 
                                                pdata = as.double(preddata),
                                                as.integer(PSAMPLE),
                                                as.double(x$coefficients[splineindex:(splineindex+numsplines-1)]),
                                                as.double(x$knots[splineindex:(splineindex+numsplines-1)]),
                                                as.integer(numsplines),
                                                PACKAGE = "gdm"))
      
      ## generating posterior stats from runs.
      quant.preds <- apply(plyr::ldply(predplotdat, function(x) c(x$pdata)),2,
                           function(x)quantile(x,c(.05,.5,.95),na.rm=T))
      varNam <- gdms[[1]]$predictors
      v <- max(quant.preds)
      if (v > predmax) predmax <- v
      
      splineindex <- splineindex + numsplines
    }
    
    splineindex <- 1
    numsplines <- 3
    for(i in 1:preds){  
      predplotdat <- lapply(gdms,function(x).C( "GetPredictorPlotData", 
                                                pdata = as.double(preddata),
                                                as.integer(PSAMPLE),
                                                as.double(x$coefficients[splineindex:(splineindex+numsplines-1)]),
                                                as.double(x$knots[splineindex:(splineindex+numsplines-1)]),
                                                as.integer(numsplines),
                                                PACKAGE = "gdm"))
      
      ## generating posterior stats from runs.
      quant.preds <- apply(plyr::ldply(predplotdat, function(x) c(x$pdata)),2,
                           function(x)quantile(x,c(.05,.5,.95),na.rm=T))
      varNam <- gdms[[1]]$predictors[i]
      plot(gdms[[1]]$knots))
      env_grad <- seq(from=gdms[[1]]$knots[splineindex], to=gdms[[1]]$knots[(splineindex+numsplines-1)], length=PSAMPLE)
      plot(env_grad, quant.preds[2,], 
           xlab=varNam, ylab=paste("f(", varNam, ")", sep="" ), ylim=c(0,predmax), type="n")
      polygon(c(rev(env_grad), env_grad), c(rev(quant.preds[3,]), quant.preds[1,]),
              col = 'grey90', border = NA)
      lines(env_grad, quant.preds[2,],col='dodgerblue',lwd=2)
      
      if(add_coefs==TRUE){
        
        quantcoefs <- apply(plyr::ldply(x, function(y) c(y$coefficients)),2,function(y)quantile(y,c(.05,.5,.95),na.rm=T))
        
        knots <- gdms[[1]]$knots[splineindex:(splineindex+numsplines-1)]
        # hack: we draw arrows but with very special "arrowheads"
        par(new = T)
        plot(knots, quantcoefs[2,splineindex:(splineindex+numsplines-1)],pch=16,col='tomato', axes=F, xlab=NA, ylab=NA,ylim=c(0,max(quantcoefs)))
        arrows(knots, quantcoefs[1,splineindex:(splineindex+numsplines-1)],
               knots, quantcoefs[3,splineindex:(splineindex+numsplines-1)], col='tomato',length=0, angle=90, code=3,lwd=2)
        axis(side = 4, col="tomato",col.axis="tomato",las=1)
        mtext("Spline coefficent estimates",side=4, col="tomato",line=-1.5, outer = TRUE)
      }
      splineindex <- splineindex + numsplines
    }  
    
    
    
    
    
    
  } else {
    PSAMPLE <- 200
    preddat <- rep(0,times=PSAMPLE)
    preds <- length(gdms[[1]]$predictors)
    predmin <- 0
    predmax <- 0
    splineindex <- 1
    numsplines <- 3
    
    for(i in 1:preds){  
      predplotdat <- lapply(gdms,function(x).C( "GetPredictorPlotData", 
                                                pdata = as.double(preddata),
                                                as.integer(PSAMPLE),
                                                as.double(x$coefficients[splineindex:(splineindex+numsplines-1)]),
                                                as.double(x$knots[splineindex:(splineindex+numsplines-1)]),
                                                as.integer(numsplines),
                                                PACKAGE = "gdm"))
      
      
      diffquants <- apply(apply(plyr::ldply(predplotdat, function(x) c(x$pdata)),1,diff),1,function(x)quantile(x,c(.1,.5,.9)))
      
      varNam <- gdms[[1]]$predictors[i]
      vmin <- min(diffquants)
      vmax <- max(diffquants)
      if (vmin < predmin) predmin <- vmin
      if (vmax > predmax) predmax <- vmax
      splineindex <- splineindex + numsplines
    }
    
    splineindex <- 1
    numsplines <- 3
    for(i in 1:preds){  
      predplotdat <- lapply(gdms,function(x).C( "GetPredictorPlotData", 
                                                pdata = as.double(preddata),
                                                as.integer(PSAMPLE),
                                                as.double(x$coefficients[splineindex:(splineindex+numsplines-1)]),
                                                as.double(x$knots[splineindex:(splineindex+numsplines-1)]),
                                                as.integer(numsplines),
                                                PACKAGE = "gdm"))
      
      
      diffquants <- apply(apply(plyr::ldply(predplotdat, function(x) c(x$pdata)),1,diff),1,function(x)quantile(x,c(.1,.5,.9)))
      
      varNam <- gdms[[1]]$predictors[i]
      
      env_grad <- seq(from=gdms[[1]]$knots[splineindex], to=gdms[[1]]$knots[(splineindex+numsplines-1)], length=PSAMPLE)
      
      #calc approximate derivates.
      dX <- diff(env_grad)
      dYl <- diffquants[1,]/dX
      dYm <- diffquants[2,]/dX
      dYu <- diffquants[3,]/dX
      dX <- rowMeans(embed(env_grad,2))
      plot(dX, dYm, 
           xlab=varNam, ylab=paste("f'(", varNam, ")", sep="" ), type="n",
           ylim=c(min(dYl),max(dYu)))
      polygon(c(rev(dX), dX), c(rev(dYu), dYl),
              col = 'grey90', border = NA)
      lines(dX, dYm,col='tomato',lwd=2)
      splineindex <- splineindex + numsplines
    }
  }
}


bbgdm.transform <- function(model, data){
  #################
  ##lines used to quickly test function
  #model <- gdmModel 
  #data <- climCurrExt
  #data <- cropRasts[[3:nlayers(cropRasts)]]
  #################
  options(warn.FPU = FALSE)
  rastDat <- NULL
  dataCheck <- class(data)
  
  ##error checking of inputs
  ##checks to make sure a gdm model is given
  if(class(model)[1]!="bbgdm"){
    stop("model argument must be a gdm model object")
  }
  ##checks to make sure data is a correct format
  if(!(dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick" | dataCheck=="data.frame")){
    stop("Data to be transformed must be either a raster object or data frame")
  }
  
  ##checks rather geo was T or F in the model object
  geo <- model$gdms[[1]]$geo
  
  ##turns raster data into dataframe
  if(dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick"){
    ##converts the raster object into a dataframe, for the gdm transformation
    rastDat <- data
    data <- rasterToPoints(rastDat)
    ##determines the cell number of the xy coordinates
    rastCells <- cellFromXY(rastDat, xy=data[,1:2]) 
    
    ##checks for NA in the 
    checkNAs <- as.data.frame(which(is.na(data), arr.ind=T))
    if(nrow(checkNAs)>0){
      warning("After extracting raster data, NAs found from one or more layers. Removing NAs from data object to be transformed.")
      data <- na.omit(data)
      rastCells <- rastCells[-c(checkNAs$row)]
    }
    
    ##if geo was not T in the model, removes the coordinates from the data frame
    if(geo==FALSE){
      data <- data[,3:ncol(data)]
    }
  }
  
  sizeVal <- 10000000
  ##sets up the data to be transformed into pieces to be transformed
  holdData <- data
  fullTrans <- matrix(0,nrow(holdData),ncol(holdData))
  rows <- nrow(holdData)
  istart <- 1
  iend <- min(sizeVal,rows)
  ##to prevent errors in the transformation of the x and y values when geo is a predictor,
  ##extracts the rows with the minimum and maximum x and y values, these rows will be added
  ##onto the "chuck" given to transform, and then immediately removed after the transformation,
  ##this makes sure that the c++ code will always have access to the minimum and maximum 
  ##x and y values
  if(geo==TRUE){
    if(dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick"){
      xMaxRow <- holdData[which.max(holdData[,"x"]),]
      xMinRow <- holdData[which.min(holdData[,"x"]),]
      yMaxRow <- holdData[which.max(holdData[,"y"]),]
      yMinRow <- holdData[which.min(holdData[,"y"]),]
    }
  }
  
  ##transform the data based on the gdm
  ##part of a loop to prevent memory errors 
  while(istart < rows){
    ##Call the dll function
    data <- holdData[istart:iend,]
    ##adds coordinate rows to data to be transformed
    if((dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick") & geo==TRUE){
      data <- rbind(xMaxRow, xMinRow, yMaxRow, yMinRow, data)
    }
    transformed <- matrix(0,nrow(data),ncol(data))
    z <- .C( "GDM_TransformFromTable",
             as.integer(nrow(data)), 
             as.integer(ncol(data)),
             as.integer(model$gdms[[1]]$geo),
             as.integer(length(model$gdms[[1]]$predictors)), 
             as.integer(model$gdms[[1]]$splines),             
             as.double(model$gdms[[1]]$knots),             
             as.double(model$median.coefs),
             as.matrix(data),
             trandata = as.double(transformed),
             PACKAGE = "gdm")
    
    ## Convert transformed from a vector into a dataframe before returning...
    nRows <- nrow(data)
    nCols <- ncol(data)
    
    ## z$trandata is the transformed data vector created
    myVec <- z$trandata
    pos <- 1
    ##fills out dataframe with transformed values
    for (i in seq(from = 1, to = nCols, by = 1)) {
      tmp <- myVec[seq(from=pos, to=pos+nRows-1)]
      transformed[,i] <- tmp
      pos <- pos + nRows
    }
    
    ##remove the coordinate rows before doing anything else
    if((dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick") & geo==TRUE){
      transformed <- transformed[-c(1:4),]
    }
    
    ##places the transformed values into the readied data frame 
    fullTrans[istart:iend,] <- transformed
    istart <- iend + 1
    iend <- min(istart + (sizeVal-1), rows)
  }
  
  ##if wanted output data as raster, provides maps raster, or output table
  if(dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick"){
    ##maps the transformed data back to the input rasters
    rastLay <- rastDat[[1]]
    rastLay[] <- NA
    outputRasts <- stack()
    for(nn in 1:ncol(fullTrans)){
      #print(nn)
      #nn=1
      holdLay <- rastLay
      holdLay[rastCells] <- fullTrans[,nn]
      #holdLay[rastCells] <- holdData[,nn]
      
      outputRasts <- stack(outputRasts, holdLay)
    }
    ##renames raster layers to be the same as the input
    if(geo){
      names(outputRasts) <- c("xCoord", "yCoord", names(rastDat))
    } else {
      names(outputRasts) <- names(rastDat)
    }
    
    ##get the predictors with non-zero sum of coefficients      
    splineindex <- 1
    predInd <- NULL
    for(i in 1:length(model$gdms[[1]]$predictors)){  
      #i <- 1
      ##only if the sum of the coefficients associated with this predictor is > 0.....
      numsplines <- model$gdms[[1]]$splines[i]
      if(sum(model$median.coefs[splineindex:(splineindex+numsplines-1)])>0){
        predInd <- c(predInd, i)
      }
      splineindex <- splineindex + numsplines
    }
    if(geo){
      predInd <- c(1,2,predInd[-1]+1)
    }
    
    outputRasts <- outputRasts[[predInd]]
    
    ##returns rasters
    return(outputRasts)
  }else{
    if(is.null(rastDat)){
      ##if not raster data, sends back the transformed data
      colnames(fullTrans) <- colnames(data)
      return(fullTrans)
    }else{
      ##returns only the transformed variable data as a table, and the cells with which to map to
      colnames(fullTrans) <- colnames(data)
      return(list(fullTrans, rastCells))
    }
  }
}

residuals.bbgdm <- function (object){
  y <- object$gdms[[1]]$observed  # offset <- object$offset
  preds<- t(plyr::ldply(object$gdms,function(x)x$predicted))
  pi <- apply(preds,1,mean)
  a <- pbinom(y-1, 1, pi)#-1
  b <- pbinom(y, 1, pi)
  u <- runif(n = length(y), min = a, max = b)
  # u[u==1] <- u[u==1]-1e-5
  # u[u==0] <- u[u==0]+1e-5
  res <- qnorm(u)
  structure(list(res=res,pi=pi,y=y))
}

plot.residuals <- function(x,...){
  qqnorm(x$res,ylab='Random Quantile Residuals',main = "")
  qqline(x$res, col = 'red')
  hist(x$res, xlab = "Random Quantile Residuals",main = "",...)
  plot(x$res,x$pi,xlab="Predicted Dissimilarity",ylab="Random Quantile Residuals",...)
  plot(x$pi,x$y,xlab="Predicted Dissimilarity",ylab="Observed Dissimilarity",...)
}

bbgdm.wald.test <- function(object,H0=0){
  # H0: Hypothesis test = 0
  # IM: Identiy Matrix
  # beta: parameter estimates from model
  # vcov: Variance-covariance matrix estimated from BB
  
  A <- plyr::ldply(object$gdms,function(x)c(x$intercept,x$coefficients))# all.coefs.se #matrix of B bootstrap coeficient estimates.
  esti.var <- var(A) #make sure this is a matrix
  beta <- c(object$median.intercept,object$median.coefs) #medians of coef estimates
  
  #Intercept
  intercept_IM <- matrix(c(1,rep(0,length(beta)-1)),nrow=1)
  wd_inter <- t(intercept_IM%*%beta-H0) %*% solve(intercept_IM%*%esti.var%*%t(intercept_IM))%*%(intercept_IM%*%beta-H0)
  pval_i = 1-pchisq(wd_inter,1)
  
  #Splines
  # beta <- c(object$median.coefs) #medians of coef estimates
  splineLength <- object$gdms[[1]]$splines
  val1 <- seq(2,length(beta),splineLength[1])
  val2 <- seq(1+splineLength[1],length(beta),splineLength[1])
  wd_vals <- matrix(NA,length(val1)+1,3)
  wd_vals[1,]<- c(wd_inter,1,pval_i)
  for(i in 1:length(splineLength)){
    w <- splineLength[1]
    L <- matrix(rep(0, length(beta) * w), ncol = length(beta))
    Terms <- seq(val1[i],val2[i],1)
    for (ii in 1:w) L[ii, Terms[ii]] <- 1
    vcov2 <- try(solve(L %*% esti.var %*% t(L)),silent = TRUE)
    if ((class(vcov2)=="try-error")||(class(vcov2)=="try-error"))
    {
      wd <- 0
    }else{
      wd <- t(L %*% beta - H0) %*% vcov2 %*% (L %*% beta - H0)
    }
    pv <- 1 - pchisq(wd, df = w)
    wd_vals[1+i,]<- c(wd,w,pv)
  }
  colnames(wd_vals) <- c("bbgdm_W","bbgdm_df","bbgdm_p-value")
  # if(object$gdm[[1]]$geo){ 
  rownames(wd_vals)<-c('intercept',object$gdms[[1]]$predictors)
  # } else { 
  # rownames(wd_vals)<-c('intercept',object$gdms[[1]]$predictors)
  # }
  wd_vals
}

negexp<- function(){
  linkfun <- function(mu) -log(1-mu)
  linkinv <- function(eta) 1-exp(-eta)
  mu.eta <- function(eta) exp(-eta)
  valideta <- function(eta) all(is.finite(eta))
  link <- paste0("negexp")
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, name = link),
            class = "link-glm")
}


appendList <- function (x, val) {
  stopifnot(is.list(x), is.list(val))
  xnames <- names(x)
  for (v in names(val)) {
    x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
      appendList(x[[v]], val[[v]])
    else c(x[[v]], val[[v]])
  }
  x
}

beta_fd<- function(comm= comm_parsed, trait = trait_shuffled){

    ntasks<-nrow(comm)
  
    cores <- 6
    cl <- makeSOCKcluster(cores)
    registerDoSNOW(cl)
    
  alpha.FD <- foreach(
    i = 1:ntasks,
    .combine = hypervolume::hypervolume_join,
    .multicombine = TRUE,
    .errorhandling = 'stop'
  ) %dopar% {
    .libPaths(c("/projappl/project_2005062/project_rpackages", .libPaths()))
    if(!require("hypervolume")) {install.packages("hypervolume")}
    abun <- comm[i, which(comm[i, ] > 0)]
    w <- as.numeric(abun/sum(abun))
    hv <- hypervolume::hypervolume_gaussian(trait[names(comm[i, which(comm[i, ] > 0)]), ],
                                            weight= w,
                                            verbose = FALSE,
                                            name = rownames(comm)[i])
  }
  
  name_sites<-sapply(seq_along(alpha.FD@HVList), function(i) alpha.FD@HVList[[i]]@Name)
  
  pairwise_beta <-  foreach(i=1:ntasks,
                            .combine=appendList) %:%
    foreach(j=i:ntasks, .combine=appendList) %do% {
      if(!require("hypervolume")) {install.packages("hypervolume")}
      hyperSet <- hypervolume::hypervolume_set(
        alpha.FD@HVList[[i]], alpha.FD@HVList[[j]],      
        check.memory = FALSE, verbose = FALSE, num.points.max = 10000)
      union <- hyperSet[[4]]@Volume
      unique1 <- hyperSet[[5]]@Volume
      unique2 <- hyperSet[[6]]@Volume
      
      union <- 2 * union - unique1 - unique2
      Btotal <- (unique1 + unique2)/union
      Brepl <- 2 * min(unique1, unique2)/union
      Brich <- abs(unique1 - unique2)/union
      output<-list(Btotal=Btotal,Brepl=Brepl,Brich=Brich)
      return(output)
    }
  
  output_beta<-list()
  output_beta$Btotal <- matrix(NA,ntasks,ntasks)
  output_beta$Brepl <- matrix(NA,ntasks,ntasks)
  output_beta$Brich <- matrix(NA,ntasks,ntasks)
  
  output_beta$Btotal[lower.tri(output_beta$Btotal,diag=TRUE)] <-round(pairwise_beta$Btotal,3)
  output_beta$Brepl [lower.tri(output_beta$Btotal,diag=TRUE)] <- round(pairwise_beta$Brich,3)
  output_beta$Brich [lower.tri(output_beta$Btotal,diag=TRUE)] <- round(pairwise_beta$Brepl,3)
  Fbeta <- lapply(output_beta, as.dist)
  names(Fbeta$Btotal) <- names(Fbeta$Brepl) <- names(Fbeta$Brich) <- name_sites
  
  return(Fbeta)
}

run_bbgdm <- function(x, pred = predictors){

  ID <- names(x)
  df<- data.frame(as.matrix(x))
  df<- cbind(ID,df)
  
  ###### EXCLUDE SITES IN PREDICTORS BUT NOT IN METRIC
  pred <- pred[pred$ID %in% names(x),]
  pred <- droplevels(pred)
  
  ###### GET PREDICTORS
  predTab <- data.frame(pred)
  
  ###### FORMAT TO GDM TABLE
  gdmTab <- formatsitepair(df,
                           bioFormat=3,
                           XColumn="decimalLongitude", 
                           YColumn="decimalLatitude",
                           predData=predTab, 
                           siteColumn="ID")

  ##### REMOVE NA
  gdmTab <- gdmTab[complete.cases(gdmTab),]
  
  
  ###### FIT BBGDM
  
  # fit bbgdm
  options(warn = 1) # turn off warnings
  mod1 <- bbgdm(data= gdmTab, geo=TRUE, bootstraps=1000, ncores=1, knots=NULL, splines= NULL) 
  
  ###### EXTRACT OUTPUTS
  
  ###### Get Wald Test results and signs of predictors
  waldtest <- bbgdm.wald.test(mod1)
  
  results <- list(beta_metric = x,
                  wald_value = waldtest,# waldtest results
                  model = mod1) # actual model objects
  
  return(results)
}


combine_custom<- function(list1,list2) {
  
  if (is.null(names(list1))|is.null(names(list2))) { 
    ls<- append(list1, list(list2))
  } else { ls<- list(list1,list2)}
  return(ls)
}
results$Trial_2$Btotal$model$median.coefs


gdm)
#######################################################################################################################