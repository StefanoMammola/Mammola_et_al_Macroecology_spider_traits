make_tab <- function(x, pred) {
  
  gdmTab <- x %>%
    as.matrix() %>%
    reshape2::melt() %>%
    set_names(c("s2","s1","distance")) %>%
    mutate(across(c(s2,s1),as.character)) %>% 
    left_join(pred %>%
                set_names(paste0("s1.",names(pred))), by= c("s1" = "s1.site")) %>%
    left_join(pred %>%
                set_names(paste0("s2.",names(pred))), by =c("s2" = "s2.site")) %>%
    mutate(grp = paste0(pmin(s1,s2),"-",pmax(s1,s2))) %>%
    distinct(grp, .keep_all=TRUE) %>%
    filter(s1 != s2) %>%
    mutate(weights = 1) %>%
    relocate(s2,s1,grp,distance,weights,s1.x,s1.y,s2.x,s2.y) %>%
    rename(s1.xCoord = "s1.x",
           s2.xCoord = "s2.x",
           s1.yCoord = "s1.y",
           s2.yCoord = "s2.y") %>%
    column_to_rownames("grp") %>%
    dplyr::select(-s2,-s1) %>% 
    drop_na()
  class(gdmTab) <- c("gdmData", "data.frame")
  return(gdmTab)
}

run_bbgdm <- function(x, pred, boot=999){
  
  gdmTab <- make_tab(x, pred)
  
  ###### FIT BBGDM
  
  # fit bbgdm
  
  options(warn = 1) # turn off warnings
  mod1 <- bbgdm(data= gdmTab, geo=TRUE, bootstraps=boot, ncores=1, knots=NULL, splines= NULL) 
  
  ###### EXTRACT OUTPUTS
  
  ###### Get Wald Test results and signs of predictors
  waldtest <- bbgdm.wald.test(mod1)
  ##### Get spline values 
  index<- which(!sapply(mod1$gdms,is.null, simplify=TRUE))[1] #retrieve index for baseline model 
  gdms <- mod1$gdms
  PSAMPLE <- 200
  preddata <- rep(0,times=PSAMPLE)
  preds <- length(gdms[[index]]$predictors)
  pred.names <- gdms[[index]]$predictors
  predmax <- 0
  splineindex <- 1
  numsplines <- 3
  
  curves <- list()
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
    
    env_grad <- matrix(seq(from=gdms[[index]]$knots[splineindex], to=gdms[[index]]$knots[(splineindex+numsplines-1)], length=PSAMPLE),ncol=1)
    colnames(env_grad) <- pred.names[i]
    #I added some modifications to this piece of code so we can make a table ready for ggplot
    curves[[i]] <- data.frame(env_grad,t(quant.preds))
    curves[[i]]$variable <- rep(colnames(curves[[i]])[1],PSAMPLE)
    colnames(curves[[i]])[1]<- "value"
    splineindex <- splineindex + numsplines
  }  
  
  r2_all_models<- sapply(gdms, function(x) x$explained)
  r2_median <-median(unlist(r2_all_models))
  
  results <- list(curves = curves, #save models 
                  wald_value = waldtest,
                  r2 = r2_median,
                  mod = mod1)  #waldtest results
  return(results)
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
  nsites <-  data %>% 
    rownames_to_column("sites") %>% 
    dplyr::select(sites) %>% 
    separate(sites, into=c("s1","s2"), sep="-") %>% 
    distinct(s1) %>% 
    pull() %>%
    length %>% 
    {.+1}
  
  mods <- surveillance::plapply(1:bootstraps, bb_apply, nsites, data, geo, splines, knots, .parallel = ncores, .verbose=FALSE)
  
  median.intercept <- apply(plyr::ldply(mods, function(x) c(x$intercept), .progress = "none"),2,median,na.rm=TRUE)
  quantiles.intercept <- apply(plyr::ldply(mods, function(x) c(x$intercept), .progress = "none"),2, function(x) quantile(x,c(.05,.25,.5,.75, .95),na.rm=TRUE))  
  
  median.coefs <- apply(plyr::ldply(mods, function(x) c(x$coefficients), .progress = "none"),2,median,na.rm=TRUE)
  quantiles.coefs <- apply(plyr::ldply(mods, function(x) c(x$coefficients), .progress = "none"),2, function(x) quantile(x,c(.05,.25,.5,.75, .95),na.rm=TRUE))
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
  
  gdms <- mods$gdms
  if(plot_derivate==FALSE){
    PSAMPLE <- 200
    preddata <- rep(0,times=PSAMPLE)
    preds <- length(gdms[[1]]$predictors)
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
  preds<- t(plyr::ldply(object$gdms,function(x){x$predicted}, .progress = "none"))
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
  index<- which(!sapply(object$gdms,is.null, simplify=TRUE))[1] #retrieve index for baseline model
  A <- plyr::ldply(object$gdms,function(x){c(x$intercept,x$coefficients)}, .progress = "none")# all.coefs.se #matrix of B bootstrap coeficient estimates.
  esti.var <- var(A) #make sure this is a matrix
  beta <- c(object$median.intercept,object$median.coefs) #medians of coef estimates
  
  
  
  #Intercept
  intercept_IM <- matrix(c(1,rep(0,length(beta)-1)),nrow=1)
  wd_inter <- t(intercept_IM%*%beta-H0) %*% solve(intercept_IM%*%esti.var%*%t(intercept_IM))%*%(intercept_IM%*%beta-H0)
  pval_i = 1-pchisq(wd_inter,1)
  
  #Splines
  # beta <- c(object$median.coefs) #medians of coef estimates
  splineLength <- object$gdms[[index]]$splines
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
  rownames(wd_vals)<-c('intercept',object$gdms[[index]]$predictors)
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