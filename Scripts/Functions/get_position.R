get_position<-function (x, choices = c(1, 2), labels, arrow.mul, at = c(0, 
                                                          0), axis = FALSE, p.max = NULL, col = "blue", bg, add = TRUE, 
          ...) {
  formals(arrows) <- c(formals(arrows), alist(... = ))
  labs <- list(v = rownames(x$vectors$arrows), f = rownames(x$factors$centroids))
  if (!missing(labels)) {
    if (is.list(labels)) {
      if (!is.null(labs$v) && !is.null(labels$vectors)) 
        labs$v <- labels$vectors
      if (!is.null(labs$f) && !is.null(labels$factors)) 
        labs$f <- labels$factors
    }
    else {
      if (!is.null(x$vectors) && !is.null(x$factors)) 
        stop("needs a list with both 'vectors' and 'factors' labels")
      if (!is.null(x$factors)) 
        labs$f <- labels
      else labs$v <- labels
    }
  }
  vect <- NULL
  if (!is.null(p.max)) {
    if (!is.null(x$vectors)) {
      take <- x$vectors$pvals <= p.max
      x$vectors$arrows <- x$vectors$arrows[take, , drop = FALSE]
      labs$v <- labs$v[take]
      x$vectors$r <- x$vectors$r[take]
      if (nrow(x$vectors$arrows) == 0) 
        x$vectors <- vect <- NULL
    }
    if (!is.null(x$factors)) {
      tmp <- x$factors$pvals <= p.max
      nam <- names(tmp)[tmp]
      take <- x$factors$var.id %in% nam
      x$factors$centroids <- x$factors$centroids[take, 
                                                 , drop = FALSE]
      labs$f <- labs$f[take]
      if (nrow(x$factors$centroids) == 0) 
        x$factors <- NULL
    }
  }
  if (!is.null(x$vectors)) {
    vect <- sqrt(x$vectors$r) * x$vectors$arrows[, choices, 
                                                 drop = FALSE]
    if (missing(arrow.mul)) {
      if (!add) 
        arrow.mul <- 1
      else arrow.mul <- ordiArrowMul(vect, at = at)
    }
    if (axis) {
      maxarr <- round(sqrt(max(x$vectors$r)), 1)
      ax <- -c(-1, 0, 1) * arrow.mul * maxarr
    }
    vect <- arrow.mul * vect
    vect <- sweep(vect, 2, at, "+")
    if (add) {
      vtext <- ordiArrowTextXY(vect, labs$v, rescale = FALSE, 
                               at = at)
    }
  }
  if (!add) {
    plot.new()
    if (is.null(vect) || is.null(x$factors)) {
      xstack <- rbind(vect, x$factors$centroids)
      plot.window(xlim = range(xstack[, 1], at[1]), ylim = range(xstack[, 
                                                                        2], at[2]), asp = 1, ...)
    }
    else {
      plot.window(xlim = range(x$factors$centroids[, 1], 
                               at[1]), ylim = range(x$factors$centroids[, 2], 
                                                    at[2]), asp = 1)
      vfill <- 0.75
      arrow.mul <- ordiArrowMul(vect, at = at, fill = 1)
      vect <- arrow.mul * vect
    }
    sw <- strwidth(c(labs$v, labs$f), ...)/2
    sh <- strheight(c(labs$v, labs$f), ...)
    xstack <- rbind(x$factors$centroids, vect)
    xlim <- range(xstack[, 1] + sw, xstack[, 2] - sw)
    ylim <- range(xstack[, 2] + sh, xstack[, 2] - sh)
    plot.window(xlim = xlim, ylim = ylim, asp = 1, ...)
    if (!is.null(vect)) {
      arrow.mul <- ordiArrowMul(vect, at = at, fill = 1)
      vect <- arrow.mul * vect
      vtext <- ordiArrowTextXY(vect, labs$v, at = at, 
                               rescale = FALSE, ...)
      sw <- strwidth(labs$v, ...)/2
      sh <- strheight(labs$v, ...)
      xlim <- range(xlim, vtext[, 1] + sw, vtext[, 1] - 
                      sw)
      ylim <- range(xlim, vtext[, 2] + sh, vtext[, 2] - 
                      sh)
      plot.window(xlim = xlim, ylim = ylim, asp = 1, ...)
    }
    axis(side = 1, ...)
    axis(side = 2, ...)
    box(...)
    alabs <- colnames(vect)
    title(..., ylab = alabs[2], xlab = alabs[1])
  }
  if (!is.null(vect)) {
    arrows(at[1], at[2], vect[, 1], vect[, 2], len = 0.05, 
           col = col)
    if (missing(bg)) 
      text(vtext, labs$v, col = col, ...)
    else ordilabel(vtext, labels = labs$v, col = col, fill = bg, 
                   ...)
  }
  if (!is.null(x$factors)) {
    if (missing(bg)) 
      text(x$factors$centroids[, choices, drop = FALSE], 
           labs$f, col = col, ...)
    else ordilabel(x$factors$centroids[, choices, drop = FALSE], 
                   labels = labs$f, col = col, fill = bg, ...)
  }
  if (axis && !is.null(vect)) {
    axis(3, at = ax + at[1], labels = c(maxarr, 0, maxarr), 
         col = col)
    axis(4, at = ax + at[2], labels = c(maxarr, 0, maxarr), 
         col = col)
  }
  return(rbind(vect,x$factors$centroids))
  
}
