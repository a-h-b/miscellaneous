annaSclass <- function (dfxy, fac, wt = rep(1, length(fac)), xax = 1, yax = 2, 
                        cstar = 1, cellipse = 1.5, axesell = TRUE, label = levels(fac), 
                        clabel = 1, cpoint = 1, pch = 20, col = rep(1, length(levels(fac))), 
                        xlim = NULL, ylim = NULL, grid = TRUE, addaxes = TRUE, origin = c(0, 
                                                                                          0), include.origin = TRUE, sub = "", csub = 1, possub = "bottomleft", 
                        colline = NULL, colfra = NULL,
                        cgrid = 1, pixmap = NULL, contour = NULL, area = NULL, add.plot = FALSE, mar=c(0.1,0.1,0.1,0.1), lwd=1,xlab="",ylab="",mgp=c(0.4,0.5,0),cex.lab=1.5,
                        main = "",cex.main=1.6,singleLab = NULL,boxcol="black") 
{
  f1 <- function(cl) {
    n <- length(cl)
    cl <- as.factor(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1:n) + n * (unclass(cl) - 1)] <- 1
    dimnames(x) <- list(names(cl), levels(cl))
    data.frame(x)
  }
  scatterutil.base <- function (dfxy, xax, yax, xlim, ylim, grid, addaxes, cgrid, include.origin, 
                                origin, sub, csub, possub, pixmap, contour, area, add.plot,xlab,ylab,mgp,cex.lab,main,cex.main) 
  {
    df <- data.frame(dfxy)
    if (!is.data.frame(df)) 
      stop("Non convenient selection for df")
    if ((xax < 1) || (xax > ncol(df))) 
      stop("Non convenient selection for xax")
    if ((yax < 1) || (yax > ncol(df))) 
      stop("Non convenient selection for yax")
    x <- df[, xax]
    y <- df[, yax]
    if (is.null(xlim)) {
      x1 <- x
      if (include.origin) 
        x1 <- c(x1, origin[1])
      x1 <- c(x1 - diff(range(x1)/10), x1 + diff(range(x1))/10)
      xlim <- range(x1)
    }
    if (is.null(ylim)) {
      y1 <- y
      if (include.origin) 
        y1 <- c(y1, origin[2])
      y1 <- c(y1 - diff(range(y1)/10), y1 + diff(range(y1))/10)
      ylim <- range(y1)
    }
    if (!is.null(pixmap)) {
      if (is.null(class(pixmap))) 
        pixmap <- NULL
      if (is.na(charmatch("pixmap", class(pixmap)))) 
        pixmap <- NULL
    }
    if (!is.null(contour)) {
      if (!is.data.frame(contour)) 
        contour <- NULL
      if (ncol(contour) != 4) 
        contour <- NULL
    }
    if (!is.null(area)) {
      if (!is.data.frame(area)) 
        area <- NULL
      if (!is.factor(area[, 1])) 
        area <- NULL
      if (ncol(area) < 3) 
        area <- NULL
    }
    if (!add.plot) 
      plot.default(0, 0, type = "n", asp = 1, xlab = xlab, ylab = ylab, 
                   xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim, 
                   xaxs = "i", yaxs = "i", frame.plot = FALSE,mgp=mgp,cex.lab=cex.lab,main=main,cex.main=cex.main,font.main=1)
    if (!is.null(pixmap)) {
      plot(pixmap, add = TRUE)
    }
    if (!is.null(contour)) {
      apply(contour, 1, function(x) segments(x[1], x[2], x[3], 
                                             x[4], lwd = 1))
    }
    if (grid & !add.plot) 
      scatterutil.grid(cgrid)
    if (addaxes & !add.plot) 
      abline(h = 0, v = 0, lty = 1)
    if (!is.null(area)) {
      nlev <- nlevels(area[, 1])
      x1 <- area[, 2]
      x2 <- area[, 3]
      for (i in 1:nlev) {
        lev <- levels(area[, 1])[i]
        a1 <- x1[area[, 1] == lev]
        a2 <- x2[area[, 1] == lev]
        polygon(a1, a2)
      }
    }
    if (csub > 0) 
      scatterutil.sub(sub, csub, possub)
    return(list(x = x, y = y))
  }
  
  opar <- par(mar = par("mar"))
  par(mar=mar)
  on.exit(par(opar))
  dfxy <- data.frame(dfxy)
  if (!is.data.frame(dfxy)) 
    stop("Non convenient selection for dfxy")
  if (any(is.na(dfxy))) 
    stop("NA non implemented")
  if (!is.factor(fac)) 
    stop("factor expected for fac")
  dfdistri <- f1(fac) * wt
  coul = col
  if(is.null(colline)){
    colline <- coul
  }
  w1 <- unlist(lapply(dfdistri, sum))
  dfdistri <- t(t(dfdistri)/w1)
  coox <- as.matrix(t(dfdistri)) %*% dfxy[, xax]
  cooy <- as.matrix(t(dfdistri)) %*% dfxy[, yax]
  if (nrow(dfxy) != nrow(dfdistri)) 
    stop(paste("Non equal row numbers", nrow(dfxy), nrow(dfdistri)))
  coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax, 
                          xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, 
                          cgrid = cgrid, include.origin = include.origin, origin = origin, 
                          sub = sub, csub = csub, possub = possub, pixmap = pixmap, 
                          contour = contour, area = area, add.plot = add.plot,xlab=xlab,ylab=ylab,mgp=mgp,cex.lab=cex.lab,
                          main=main,cex.main=cex.main)
  if (cpoint > 0) 
    for (i in 1:ncol(dfdistri)) {
      pch <- rep(pch, length = nrow(dfxy))
      if(any(pch %in% c(21:25))&!is.null(colfra)){
        points(coo$x[dfdistri[, i] > 0], coo$y[dfdistri[, 
                                                        i] > 0], pch = pch[dfdistri[, i] > 0], cex = par("cex") * 
                 cpoint, bg = coul[i],col=colfra[i])
      }else{
        points(coo$x[dfdistri[, i] > 0], coo$y[dfdistri[, 
                                                        i] > 0], pch = pch[dfdistri[, i] > 0], cex = par("cex") * 
                 cpoint, bg = coul[i],col=coul[i])
      }
    }
  if(!is.null(singleLab)){
    text(coo$x,coo$y,singleLab,pos=3,cex=0.6)
  }
  if (cstar > 0) 
    for (i in 1:ncol(dfdistri)) {
      scatterutil.star(coo$x, coo$y, dfdistri[, i], cstar = cstar, 
                       colline[i])
    }
  if (cellipse > 0) 
    for (i in 1:ncol(dfdistri)) {
      scatterutil.ellipse(coo$x, coo$y, dfdistri[, i], 
                          cellipse = cellipse, axesell = axesell, colline[i])
    }
  if (clabel > 0) 
    scatterutil.eti(coox, cooy, label, clabel, coul = col)
  box(lwd=lwd,col=boxcol)
  invisible(match.call())
}

library(gtools)
heatmap.2a <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
                        distfun = dist, hclustfun = hclust, dendrogram = c("both", 
                                                                           "row", "column", "none"), symm = FALSE, 
                        scale = c("none", "row", "column"), na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks, 
                        symbreaks = min(x < 0, na.rm = TRUE) || scale != "none", col = "heat.colors", colsep, rowsep, 
                        sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
                        notecol = "cyan", na.color = par("bg"), trace = c("column", "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
                        vline = median(breaks), linecol = tracecol, margins = c(5,5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
                        cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0, NA), adjCol = c(NA,0),
                        offsetRow = 0.5, offsetCol = 0.5, key = TRUE, keysize = 1.5, 
                        density.info = c("histogram", "density", "none"), denscol = tracecol, 
                        symkey = min(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25, 
                        main = NULL, xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, 
                        lwid = NULL, extrafun = NULL, font.col=1, font.row=1, keyName="Value", ...) 
{
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale)) 
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col)) 
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none")) 
    warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv)) 
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv)) 
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv)) 
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) 
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2) 
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote)) 
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv)) 
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv)) 
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc) 
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm) 
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm) 
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow)) 
    labRow <- if (is.null(rownames(x))) 
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol)) 
    labCol <- if (is.null(colnames(x))) 
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 
      1) {
    if (missing(col) || is.function(col)) 
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks) 
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function") 
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei)) 
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid)) 
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    if (!missing(ColSideColors)) {
      if (is.vector(ColSideColors)) ColSideColors <- as.matrix(t(ColSideColors))
      if (!is.character(ColSideColors) || ncol(ColSideColors) != nc) 
        stop("'ColSideColors' must be a character vector of length ncol(x) or character matrix with ncol ncol(x)")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      lhei <- c(lhei[1], 0.12*nrow(ColSideColors), lhei[2])
    }
    if (!missing(RowSideColors)) {
      if (is.vector(RowSideColors)) RowSideColors <- as.matrix(RowSideColors)
      if (!is.character(RowSideColors) || nrow(RowSideColors) != nr) 
        stop("'RowSideColors' must be a character vector of length nrow(x) or character matrix with nrow nrow(x)")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -  1), 1), lmat[, 2] + 1)
      lwid <- c(lwid[1], 0.12*ncol(RowSideColors), lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  if (length(lhei) != nrow(lmat)) 
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat)) 
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    image(matrix(1:(ncol(RowSideColors)*nr),ncol=nr,nrow=ncol(RowSideColors),byrow=F), col = t(RowSideColors[rowInd,]), axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    image(matrix(1:(nrow(ColSideColors)*nc),nrow=nc,ncol=nrow(ColSideColors),byrow=F),
          col=t(ColSideColors[,colInd]),axes = FALSE)
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr")) 
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
          c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr")) 
    retval$rowDendrogram <- ddr
  if (exists("ddc")) 
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
          col = na.color, add = TRUE)
  }
  if (is.null(srtCol)) 
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 + 
           offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1], 
         padj = adjCol[2], font=font.col)
  else {
    if (is.numeric(srtCol)) {
      if (missing(adjCol) || is.null(adjCol)) 
        adjCol = c(1, NA)
      xpd.orig <- par("xpd")
      par(xpd = NA)
      xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, 
                   tick = 0)
      text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * 
             strheight("M"), labels = labCol, adj = adjCol, 
           cex = cexCol, srt = srtCol, font=font.col)
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtCol ignored.")
  }
  if (is.null(srtRow)) {
    axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow, 
         tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2], font=font.row)
  }
  else {
    if (is.numeric(srtRow)) {
      xpd.orig <- par("xpd")
      par(xpd = NA)
      ypos <- axis(4, iy, labels = rep("", nr), las = 2, 
                   line = -0.5, tick = 0)
      text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"), 
           y = ypos, labels = labRow, adj = adjRow, cex = cexRow, 
           srt = srtRow, font=font.row)
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtRow ignored.")
  }
  if (!is.null(xlab)) 
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  if (!is.null(ylab)) 
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr)) 
    eval(substitute(add.expr))
  if (!missing(colsep)) 
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
                              xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep)) 
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
                                                      1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
                                                                                                       1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
                              col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol, 
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i - 0.5 + hline.vals, col = linecol, 
               lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote)) 
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main)) 
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row") 
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column") 
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, keyName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
    if (trace %in% c("both", "column")) {
      vline.vals <- scale01(vline, min.raw, max.raw)
      if (!is.null(vline)) {
        abline(v = vline.vals, col = linecol, lty = 2)
      }
    }
    if (trace %in% c("both", "row")) {
      hline.vals <- scale01(hline, min.raw, max.raw)
      if (!is.null(hline)) {
        abline(v = hline.vals, col = linecol, lty = 2)
      }
    }
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
                                  high = retval$breaks[-1], color = retval$col)
  if (!is.null(extrafun)) 
    extrafun()
  invisible(retval)
}

heatmap.2mn <- function(x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
                        distfun = dist, hclustfun = hclust, dendrogram = c("both", 
                                                                           "row", "column", "none"), symm = FALSE, 
                        scale = c("none", "row", "column"), na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks, 
                        symbreaks = min(x < 0, na.rm = TRUE) || scale != "none", col = "heat.colors", colsep, rowsep, 
                        sepcolor = c("white","white"), sepwidth = c(0.02, 0.02), cellnote, notecex = 1, 
                        notecol = "cyan", na.color = par("bg"), trace = c("column", "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
                        vline = median(breaks), linecol = tracecol, margins = c(5,5), ColSideColors, RowSideColors, cexRow = 0.32 + 1.6/log10(nr), 
                        cexCol = 0.32 + 1.6/log10(nc), labRow = NULL, labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0, NA), adjCol = c(NA,0.5),
                        offsetRow = 0.5, offsetCol = -0.6, key = TRUE, keysize = 27, keyheight = 100,
                        density.info = c("histogram", "density", "none"), denscol = tracecol, ColSideFac=1, RowSideFac=1,
                        symkey = min(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25, 
                        main = NULL, xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, 
                        lwid = NULL, extrafun = NULL, font.col=1, font.row=1, keyName="Value", ...) 
{
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale)) 
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  keyhigh <- keysize * keyheight / 100
  if (length(col) == 1 && is.character(col)) 
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none")) 
    warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv)) 
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv)) 
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv)) 
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) 
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2) 
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote)) 
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv)) 
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv)) 
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  } else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  } else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  }else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  } else if (identical(Colv, "Rowv")) {
    if (nr != nc) 
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  } else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm) 
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm) 
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow)) 
    labRow <- if (is.null(rownames(x))) 
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol)) 
    labCol <- if (is.null(colnames(x))) 
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  } else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col)) 
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks) 
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function") 
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei)) 
    lhei <- c(keyhigh, 100-keyhigh)
  if (missing(lwid) || is.null(lwid)) 
    lwid <- c(keysize, 100-keysize)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    if (!missing(ColSideColors)) {
      if (is.vector(ColSideColors)) ColSideColors <- as.matrix(t(ColSideColors))
      if (!is.character(ColSideColors) || ncol(ColSideColors) != nc) 
        stop("'ColSideColors' must be a character vector of length ncol(x) or character matrix with ncol ncol(x)")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      lhei <- c(lhei[1], 3*nrow(ColSideColors)*ColSideFac, lhei[2])
    }
    if (!missing(RowSideColors)) {
      if (is.vector(RowSideColors)) RowSideColors <- as.matrix(RowSideColors)
      if (!is.character(RowSideColors) || nrow(RowSideColors) != nr) 
        stop("'RowSideColors' must be a character vector of length nrow(x) or character matrix with nrow nrow(x)")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -  1), 1), lmat[, 2] + 1)
      lwid <- c(lwid[1], 2*ncol(RowSideColors)*RowSideFac, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  if (length(lhei) != nrow(lmat)) 
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat)) 
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(0.5, 0, margins[1], 0.5))
    image(matrix(1:(ncol(RowSideColors)*nr),ncol=nr,nrow=ncol(RowSideColors),byrow=F), col = t(RowSideColors[rowInd,]), axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    image(matrix(1:(nrow(ColSideColors)*nc),nrow=nc,ncol=nrow(ColSideColors),byrow=F),
          col=t(ColSideColors[,colInd]),axes = FALSE)
  }
  par(mar = c(0.5, 0, margins[1], margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr")) 
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
          c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr")) 
    retval$rowDendrogram <- ddr
  if (exists("ddc")) 
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
          col = na.color, add = TRUE)
  }
  if (is.null(srtCol)) 
    axis(3, 1:nc, labels = labCol, las = 2, line = offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1], 
         padj = adjCol[2], font=font.col)
  else {
    if (is.numeric(srtCol)) {
      if (missing(adjCol) || is.null(adjCol)) 
        adjCol = c(1, NA)
      xpd.orig <- par("xpd")
      par(xpd = NA)
      xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, 
                   tick = 0)
      text(x = xpos, y = par("usr")[3] + (offsetCol) * 
             strheight("M"), labels = labCol, adj = adjCol, 
           cex = cexCol, srt = srtCol, font=font.col)
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtCol ignored.")
  }
  if (is.null(srtRow)) {
    axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow, 
         tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2], font=font.row)
  }
  else {
    if (is.numeric(srtRow)) {
      xpd.orig <- par("xpd")
      par(xpd = NA)
      ypos <- axis(4, iy, labels = rep("", nr), las = 2, 
                   line = -0.5, tick = 0)
      text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"), 
           y = ypos, labels = labRow, adj = adjRow, cex = cexRow, 
           srt = srtRow, font=font.row)
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtRow ignored.")
  }
  if (!is.null(xlab)) 
    mtext(xlab, side = 3, line = margins[1] - 1.25)
  if (!is.null(ylab)) 
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr)) 
    eval(substitute(add.expr))
  
  if (!missing(rowsep)) 
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
                                                      1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
                                                                                                       1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
                              col = sepcolor[2], border = sepcolor[2])
  if (!missing(colsep)) 
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
                              xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                                1, lty = 1, lwd = 0.5, col = sepcolor[1], border = sepcolor[1])
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol, 
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i - 0.5 + hline.vals, col = linecol, 
               lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote)) 
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  } else plot.new()
  if (!is.null(main)) 
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    cex1=0.75*keysize/27
    par(mar = c(3, 4.7, 2.8, 0.5), cex = cex1)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    } else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv, mgp=c(3,0.6,0)*keysize/27)
    if (scale == "row") 
      mtext(side = 1, "Row Z-Score", line = 1.7, cex=cex1)
    else if (scale == "column") 
      mtext(side = 1, "Column Z-Score", line = 1.7, cex=cex1)
    else mtext(side = 1, keyName, line = 1.7, cex=cex1)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y),las=1)
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
            col = denscol)
      axis(2, at = pretty(hy,3)/max(hy) * 0.95, pretty(hy,3),las=1,mgp=c(3,0.9,0)*keysize/27)
      mtext("Color Key\nand Histogram",3,line=0.3,cex=cex1,adj=0,font=2)
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 3.5,cex=cex1)
    }
    else mtext("Color Key",3,line=0.3,cex=cex1,adj=0,font=2)
    if (trace %in% c("both", "column")) {
      vline.vals <- scale01(vline, min.raw, max.raw)
      if (!is.null(vline)) {
        abline(v = vline.vals, col = linecol, lty = 2)
      }
    }
    if (trace %in% c("both", "row")) {
      hline.vals <- scale01(hline, min.raw, max.raw)
      if (!is.null(hline)) {
        abline(v = hline.vals, col = linecol, lty = 2)
      }
    }
  }else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
                                  high = retval$breaks[-1], color = retval$col)
  if (!is.null(extrafun)) 
    extrafun()
  invisible(retval)
}

heatmap.2m <- function(x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
                       distfun = dist, hclustfun = hclust, dendrogram = c("both", 
                                                                          "row", "column", "none"), symm = FALSE, 
                       scale = c("none", "row", "column"), na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks, 
                       symbreaks = min(x < 0, na.rm = TRUE) || scale != "none", col = "heat.colors", colsep, rowsep, 
                       sepcolor = c("white","white"), sepwidth = c(0.02, 0.02), cellnote, notecex = 1, 
                       notecol = "cyan", na.color = par("bg"), trace = c("column", "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
                       vline = median(breaks), linecol = tracecol, margins = c(5,5), ColSideColors, RowSideColors, cexRow = 0.32 + 1.6/log10(nr), 
                       cexCol = 0.32 + 1.6/log10(nc), labRow = NULL, labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0, NA), adjCol = c(NA,0.5),
                       offsetRow = 0.5, offsetCol = 0.5, key = TRUE, keysize = 27, keyheight = 100,
                       density.info = c("histogram", "density", "none"), denscol = tracecol, ColSideFac=1, RowSideFac=1,
                       symkey = min(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25, 
                       main = NULL, xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, 
                       lwid = NULL, extrafun = NULL, font.col=1, font.row=1, keyName="Value", ...) 
{
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale)) 
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  keyhigh <- keysize * keyheight / 100
  if (length(col) == 1 && is.character(col)) 
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none")) 
    warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv)) 
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv)) 
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv)) 
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) 
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2) 
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote)) 
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv)) 
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv)) 
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  } else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  } else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  }else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  } else if (identical(Colv, "Rowv")) {
    if (nr != nc) 
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  } else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm) 
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm) 
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow)) 
    labRow <- if (is.null(rownames(x))) 
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol)) 
    labCol <- if (is.null(colnames(x))) 
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  } else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col)) 
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks) 
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function") 
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei)) 
    lhei <- c(keyhigh, 100-keyhigh)
  if (missing(lwid) || is.null(lwid)) 
    lwid <- c(keysize, 100-keysize)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    if (!missing(ColSideColors)) {
      if (is.vector(ColSideColors)) ColSideColors <- as.matrix(t(ColSideColors))
      if (!is.character(ColSideColors) || ncol(ColSideColors) != nc) 
        stop("'ColSideColors' must be a character vector of length ncol(x) or character matrix with ncol ncol(x)")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      lhei <- c(lhei[1], 3*nrow(ColSideColors)*ColSideFac, lhei[2])
    }
    if (!missing(RowSideColors)) {
      if (is.vector(RowSideColors)) RowSideColors <- as.matrix(RowSideColors)
      if (!is.character(RowSideColors) || nrow(RowSideColors) != nr) 
        stop("'RowSideColors' must be a character vector of length nrow(x) or character matrix with nrow nrow(x)")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -  1), 1), lmat[, 2] + 1)
      lwid <- c(lwid[1], 2*ncol(RowSideColors)*RowSideFac, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  if (length(lhei) != nrow(lmat)) 
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat)) 
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    image(matrix(1:(ncol(RowSideColors)*nr),ncol=nr,nrow=ncol(RowSideColors),byrow=F), col = t(RowSideColors[rowInd,]), axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    image(matrix(1:(nrow(ColSideColors)*nc),nrow=nc,ncol=nrow(ColSideColors),byrow=F),
          col=t(ColSideColors[,colInd]),axes = FALSE)
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr")) 
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
          c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr")) 
    retval$rowDendrogram <- ddr
  if (exists("ddc")) 
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
          col = na.color, add = TRUE)
  }
  if (is.null(srtCol)) 
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 + 
           offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1], 
         padj = adjCol[2], font=font.col)
  else {
    if (is.numeric(srtCol)) {
      if (missing(adjCol) || is.null(adjCol)) 
        adjCol = c(1, NA)
      xpd.orig <- par("xpd")
      par(xpd = NA)
      xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, 
                   tick = 0)
      text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * 
             strheight("M"), labels = labCol, adj = adjCol, 
           cex = cexCol, srt = srtCol, font=font.col)
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtCol ignored.")
  }
  if (is.null(srtRow)) {
    axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow, 
         tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2], font=font.row)
  }
  else {
    if (is.numeric(srtRow)) {
      xpd.orig <- par("xpd")
      par(xpd = NA)
      ypos <- axis(4, iy, labels = rep("", nr), las = 2, 
                   line = -0.5, tick = 0)
      text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"), 
           y = ypos, labels = labRow, adj = adjRow, cex = cexRow, 
           srt = srtRow, font=font.row)
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtRow ignored.")
  }
  if (!is.null(xlab)) 
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  if (!is.null(ylab)) 
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr)) 
    eval(substitute(add.expr))
  
  if (!missing(rowsep)) 
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
                                                      1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
                                                                                                       1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
                              col = sepcolor[2], border = sepcolor[2])
  if (!missing(colsep)) 
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
                              xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                                1, lty = 1, lwd = 0.5, col = sepcolor[1], border = sepcolor[1])
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol, 
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i - 0.5 + hline.vals, col = linecol, 
               lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote)) 
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  } else plot.new()
  if (!is.null(main)) 
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    cex1=0.75*keysize/27
    par(mar = c(3, 4.7, 2.8, 0.5), cex = cex1)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    } else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv, mgp=c(3,0.6,0)*keysize/27)
    if (scale == "row") 
      mtext(side = 1, "Row Z-Score", line = 1.7, cex=cex1)
    else if (scale == "column") 
      mtext(side = 1, "Column Z-Score", line = 1.7, cex=cex1)
    else mtext(side = 1, keyName, line = 1.7, cex=cex1)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y),las=1)
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
            col = denscol)
      axis(2, at = pretty(hy,3)/max(hy) * 0.95, pretty(hy,3),las=1,mgp=c(3,0.9,0)*keysize/27)
      mtext("Color Key\nand Histogram",3,line=0.3,cex=cex1,adj=0,font=2)
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 3.5,cex=cex1)
    }
    else mtext("Color Key",3,line=0.3,cex=cex1,adj=0,font=2)
    if (trace %in% c("both", "column")) {
      vline.vals <- scale01(vline, min.raw, max.raw)
      if (!is.null(vline)) {
        abline(v = vline.vals, col = linecol, lty = 2)
      }
    }
    if (trace %in% c("both", "row")) {
      hline.vals <- scale01(hline, min.raw, max.raw)
      if (!is.null(hline)) {
        abline(v = hline.vals, col = linecol, lty = 2)
      }
    }
  }else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
                                  high = retval$breaks[-1], color = retval$col)
  if (!is.null(extrafun)) 
    extrafun()
  invisible(retval)
}

heatmap.2pn <- function(x, PW, bgList, nw=F,pval=0.05, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
                        distfun = dist, hclustfun = hclust, dendrogram = c("both", 
                                                                           "row", "column", "none"), symm = FALSE, 
                        scale = c("none", "row", "column"), na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks, 
                        symbreaks = min(x < 0, na.rm = TRUE) || scale != "none", col = "heat.colors", colsep, rowsep, 
                        sepcolor = "white", sepwidth = c(0.02, 0.02), cellnote, notecex = 1, 
                        notecol = "cyan", na.color = par("bg"), trace = c("column", "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
                        vline = median(breaks), linecol = tracecol, margins = c(5,5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
                        cexCol = 0.32 + 1.6/log10(nc), labRow = NULL, labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0, NA), adjCol = c(NA,0.5),
                        offsetRow = 0.1, offsetCol = 0.1, key = TRUE, keysize = 27, keyheight = 100, ColSideFac=1, RowSideFac=1,
                        density.info = c("histogram", "density", "none"), denscol = tracecol, pathcut=3,
                        symkey = min(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25, writeRows=F,
                        main = NULL, xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, adjsca = 10,
                        lwid = NULL, extrafun = NULL, font.col=1, font.row=1, keyName="Value", pathFac= 1, ...) 
{
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale)) 
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  keyhigh <- keysize * keyheight / 100
  if (length(col) == 1 && is.character(col)) 
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none")) 
    warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv)) 
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv)) 
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv)) 
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) 
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2) 
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote)) 
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv)) 
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv)) 
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  } else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  } else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  }else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  } else if (identical(Colv, "Rowv")) {
    if (nr != nc) 
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  } else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm) 
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm) 
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow)) 
    labRow <- if (is.null(rownames(x))) 
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol)) 
    labCol <- if (is.null(colnames(x))) 
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  } else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col)) 
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks) 
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function") 
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei)) 
    lhei <- c(keyhigh, 100-keyhigh)
  if (missing(lwid) || is.null(lwid)) 
    lwid <- c(keysize, 100-keysize)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    if (!missing(ColSideColors)) {
      if (is.vector(ColSideColors)) ColSideColors <- as.matrix(t(ColSideColors))
      if (!is.character(ColSideColors) || ncol(ColSideColors) != nc) 
        stop("'ColSideColors' must be a character vector of length ncol(x) or character matrix with ncol ncol(x)")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      lhei <- c(lhei[1], nrow(ColSideColors)*2*ColSideFac, lhei[2])
      lmat <- cbind(lmat,c(0,0,max(lmat+1,na.rm=T)))
    }else{
      lmat <- cbind(lmat,c(0,max(lmat+1,na.rm=T)))
    }
    if (!missing(RowSideColors)) {
      if (is.vector(RowSideColors)) RowSideColors <- as.matrix(RowSideColors)
      if (!is.character(RowSideColors) || nrow(RowSideColors) != nr) 
        stop("'RowSideColors' must be a character vector of length nrow(x) or character matrix with nrow nrow(x)")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -  1), 1), lmat[, 2] + 1)
      lwid <- c(lwid[1], 2*ncol(RowSideColors)*RowSideFac, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    lwid <- c(lwid,lwid[length(lwid)]*pathFac)
  }
  if (length(lhei) != nrow(lmat)) 
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat)) 
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(0.5, 0, margins[1], 0.5))
    image(matrix(1:(ncol(RowSideColors)*nr),ncol=nr,nrow=ncol(RowSideColors),byrow=F), col = t(RowSideColors[rowInd,]), axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    image(matrix(1:(nrow(ColSideColors)*nc),nrow=nc,ncol=nrow(ColSideColors),byrow=F),
          col=t(ColSideColors[,colInd]),axes = FALSE)
  }
  par(mar = c(0.5, 0, margins[1], margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr")) 
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
          c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
  imai <- par("mai")
  retval$carpet <- x
  if (exists("ddr")) 
    retval$rowDendrogram <- ddr
  if (exists("ddc")) 
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
          col = na.color, add = TRUE)
    #imai <- par("mai")
  }
  tempParUsr <- par("usr")
  if (is.null(srtCol)) 
    axis(3, 1:nc, labels = labCol, las = 2, line = offsetCol - 0.8, 
         tick = 0, cex.axis = cexCol, hadj = adjCol[1], 
         padj = adjCol[2], font=font.col)
  else {
    if (is.numeric(srtCol)) {
      if (missing(adjCol) || is.null(adjCol)) 
        adjCol = c(1, NA)
      xpd.orig <- par("xpd")
      par(xpd = NA)
      xpos <- axis(3, 1:nc, labels = rep("", nc), las = 2, 
                   tick = 0)
      text(x = xpos, y = par("usr")[4] + (0.5 + offsetCol) * 
             strheight("M"), labels = labCol, adj = adjCol, 
           cex = cexCol, srt = srtCol, font=font.col)
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtCol ignored.")
  }
  if (!is.null(srtRow)) {
    if (is.numeric(srtRow)) {
      xpd.orig <- par("xpd")
      par(xpd = NA)
      ypos <- axis(4, iy, labels = rep("", nr), line = 0,las = 2, 
                   line = -0.5, tick = 0)
      if(writeRows){
        text(x = xpos[length(xpos)], 
             y = ypos, labels = labRow, adj = adjRow, cex = cexRow, 
             srt = srtRow, font=font.row,xpd=T)
      }
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtRow ignored.")
  }
  if (!is.null(xlab)) 
    mtext(xlab, side = 3, line = margins[1] - 1.25)
  if (!is.null(ylab)) 
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr)) 
    eval(substitute(add.expr))
  if (!missing(colsep)) 
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
                              xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                                1, lty = 1, lwd = 0.5, col = sepcolor, border = sepcolor)
  if (!missing(rowsep)) 
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
                                                      1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
                                                                                                       1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
                              col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol, 
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i - 0.5 + hline.vals, col = linecol, 
               lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote)) 
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
         col = notecol, cex = notecex)
  par(mar = c(0.5, 0, margins[1], 0))
  imai <- par("mai")
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  } else plot.new()
  if (!is.null(main)) 
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    cexO <- par("cex")
    cex1=0.75*keysize/27
    par(mar = c(3, 4.7, 2.8, 0.5), cex = cex1)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    } else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    par("cex"=cexO)
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv, mgp=c(3,0.6,0)*keysize/27,cex.axis=1.2*cexCol*27/keysize)
    if (scale == "row") 
      mtext(side = 1, "Row Z-Score", line = 1.7, cex=cex1*27/keysize)
    else if (scale == "column") 
      mtext(side = 1, "Column Z-Score", line = 1.7, cex=cex1*27/keysize)
    else mtext(side = 1, keyName, line = 1.7, cex=cex1*27/keysize)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y),las=1)
      title("Color Key\nand Density Plot",cex=cex1*27/keysize)
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
      par(cex=cexO)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
            col = denscol)
      axis(2, at = pretty(hy,3)/max(hy) * 0.95, pretty(hy,3),las=1,mgp=c(3,0.9,0)*keysize/27)
      mtext("Color Key\nand Histogram",3,line=0.3,cex=cexCol*keysize/27,adj=0,font=2)
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 3.5,cex=cex1)
      par(cex=cexO)
    }
    else mtext("Color Key",3,line=0.3,cex=cex1*27/keysize,adj=0,font=2)
    if (trace %in% c("both", "column")) {
      vline.vals <- scale01(vline, min.raw, max.raw)
      if (!is.null(vline)) {
        abline(v = vline.vals, col = linecol, lty = 2)
      }
    }
    if (trace %in% c("both", "row")) {
      hline.vals <- scale01(hline, min.raw, max.raw)
      if (!is.null(hline)) {
        abline(v = hline.vals, col = linecol, lty = 2)
      }
    }
  }else plot.new()
  if(!nw){
    enrPW <- annaEnrichPath(data.frame("koID"=labRow,stringsAsFactors=F),bgList,pval=pval)
    print(enrPW$numberInSig)
    enrPW <- enrPW[as.numeric(enrPW$numberInSig) >= pathcut,]
    print(pathcut)
  }else{
    enrPW <- annaEnrichPathNW(labRow,pval=pval)
    enrPW <- enrPW[as.numeric(enrPW$numberInSig) >= pathcut,]
  }
  print(enrPW$pathway)
  if(nrow(enrPW)>8){
    enrPW <- enrPW[order(enrPW$adjpval,decreasing=F)[1:8],]
    warning("enriched pathways truncated to 8")
  }
  enrPW <- enrPW[order(sapply(enrPW$KOs,function(x) median(which(labRow %in% unlist(strsplit(x,split=";"))))),decreasing=F),]
  
  par(mai=imai)
  plot(c(0,1),c(1,nr),type="n",yaxs='i',ann=F,axes=F)
  heiwi <- par("pin")/par("din")
  print(heiwi)
  tohei <- par("fig")[3]
  pathposyA <- seq(from=1,to=nr,length.out=length(enrPW$pwDesc)+2)
  text(0.27,pathposyA,label=c("",enrPW$pwDesc,""),adj=0,cex=cexCol)
  pathposy <- pathposyA[-c(1,length(pathposyA))]
  print(par("usr"))
  xsca <- par("usr")[1:2]
  ysca <- par("usr")[3:4]
  pushViewport(viewport(xscale=xsca, yscale=ysca,x=1,width=unit(heiwi[1],"npc"),
                        y=tohei,just=c("right","bottom"),
                        height=unit(heiwi[2],"npc")))
  # grid.points(x=rep(0,nr),
  #              y=adjsca+nr*(1:nr)/(nr+1),
  #              default.units="native",pch=".")
  for(i in 1:nrow(enrPW)){
    concol <- alpha(brewer.pal(nrow(enrPW),"Dark2")[i],0.5)
    for(ko in unlist(strsplit(enrPW$KOs[i],split=";"))){
      lry <- which(labRow==ko)
      grid.bezier(x=c(xsca[1],0.2,0.1,0.25),
                  y=c(adjsca+nr*(lry)/(nr+1),adjsca+nr*(lry)/(nr+1),
                      #pathposy[i],pathposy[i]),
                      adjsca+nr*(pathposy[i])/(nr+1),adjsca+nr*(pathposy[i])/(nr+1)),
                  default.units="native",gp=gpar(col=concol,lwd=0.7))
    }
  }
  popViewport()  
  
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
                                  high = retval$breaks[-1], color = retval$col)
  if (!is.null(extrafun)) 
    extrafun()
  invisible(retval)
}

heatmap.2p <- function(x, PW, bgList, pwn=pwn,nw=F,pval=0.05, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
                       distfun = dist, hclustfun = hclust, dendrogram = c("both", 
                                                                          "row", "column", "none"), symm = FALSE, 
                       scale = c("none", "row", "column"), na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks, 
                       symbreaks = min(x < 0, na.rm = TRUE) || scale != "none", col = "heat.colors", colsep, rowsep, 
                       sepcolor = "white", sepwidth = c(0.02, 0.02), cellnote, notecex = 1, 
                       notecol = "cyan", na.color = par("bg"), trace = c("column", "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
                       vline = median(breaks), linecol = tracecol, margins = c(5,5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
                       cexCol = 0.32 + 1.6/log10(nc), labRow = NULL, labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0, NA), adjCol = c(NA,0.5),
                       offsetRow = 0.1, offsetCol = 0.5, key = TRUE, keysize = 27, keyheight = 100, ColSideFac=1, RowSideFac=1,
                       density.info = c("histogram", "density", "none"), denscol = tracecol, pathcut=3,
                       symkey = min(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25, writeRows=F,
                       main = NULL, xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, 
                       lwid = NULL, extrafun = NULL, font.col=1, font.row=1, keyName="Value", pathFac= 1, ...) 
{
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale)) 
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  keyhigh <- keysize * keyheight / 100
  if (length(col) == 1 && is.character(col)) 
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none")) 
    warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv)) 
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv)) 
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv)) 
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) 
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2) 
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote)) 
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv)) 
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv)) 
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  } else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  } else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  }else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  } else if (identical(Colv, "Rowv")) {
    if (nr != nc) 
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  } else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm) 
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm) 
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow)) 
    labRow <- if (is.null(rownames(x))) 
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol)) 
    labCol <- if (is.null(colnames(x))) 
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  } else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col)) 
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks) 
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function") 
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei)) 
    lhei <- c(keyhigh, 100-keyhigh)
  if (missing(lwid) || is.null(lwid)) 
    lwid <- c(keysize, 100-keysize)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    if (!missing(ColSideColors)) {
      if (is.vector(ColSideColors)) ColSideColors <- as.matrix(t(ColSideColors))
      if (!is.character(ColSideColors) || ncol(ColSideColors) != nc) 
        stop("'ColSideColors' must be a character vector of length ncol(x) or character matrix with ncol ncol(x)")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      lhei <- c(lhei[1], nrow(ColSideColors)*2*ColSideFac, lhei[2])
      lmat <- cbind(lmat,c(0,0,max(lmat+1,na.rm=T)))
    }else{
      lmat <- cbind(lmat,c(0,max(lmat+1,na.rm=T)))
    }
    if (!missing(RowSideColors)) {
      if (is.vector(RowSideColors)) RowSideColors <- as.matrix(RowSideColors)
      if (!is.character(RowSideColors) || nrow(RowSideColors) != nr) 
        stop("'RowSideColors' must be a character vector of length nrow(x) or character matrix with nrow nrow(x)")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -  1), 1), lmat[, 2] + 1)
      lwid <- c(lwid[1], 2*ncol(RowSideColors)*RowSideFac, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    lwid <- c(lwid,lwid[length(lwid)]*pathFac)
  }
  if (length(lhei) != nrow(lmat)) 
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat)) 
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    image(matrix(1:(ncol(RowSideColors)*nr),ncol=nr,nrow=ncol(RowSideColors),byrow=F), col = t(RowSideColors[rowInd,]), axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    image(matrix(1:(nrow(ColSideColors)*nc),nrow=nc,ncol=nrow(ColSideColors),byrow=F),
          col=t(ColSideColors[,colInd]),axes = FALSE)
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr")) 
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
          c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
  imai <- par("mai")
  retval$carpet <- x
  if (exists("ddr")) 
    retval$rowDendrogram <- ddr
  if (exists("ddc")) 
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
          col = na.color, add = TRUE)
    imai <- par("mai")
  }
  tempParUsr <- par("usr")
  if (is.null(srtCol)) 
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 + 
           offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1], 
         padj = adjCol[2], font=font.col)
  else {
    if (is.numeric(srtCol)) {
      if (missing(adjCol) || is.null(adjCol)) 
        adjCol = c(1, NA)
      xpd.orig <- par("xpd")
      par(xpd = NA)
      xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, 
                   tick = 0)
      text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * 
             strheight("M"), labels = labCol, adj = adjCol, 
           cex = cexCol, srt = srtCol, font=font.col)
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtCol ignored.")
  }
  if (!is.null(srtRow)) {
    if (is.numeric(srtRow)) {
      xpd.orig <- par("xpd")
      par(xpd = NA)
      ypos <- axis(4, iy, labels = rep("", nr), las = 2, 
                   line = -0.5, tick = 0)
      if(writeRows){
        text(x = xpos[length(xpos)], 
             y = ypos, labels = labRow, adj = adjRow, cex = cexRow, 
             srt = srtRow, font=font.row,xpd=T)
      }
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtRow ignored.")
  }
  if (!is.null(xlab)) 
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  if (!is.null(ylab)) 
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr)) 
    eval(substitute(add.expr))
  if (!missing(colsep)) 
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
                              xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                                1, lty = 1, lwd = 0.5, col = sepcolor, border = sepcolor)
  if (!missing(rowsep)) 
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
                                                      1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
                                                                                                       1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
                              col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol, 
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i - 0.5 + hline.vals, col = linecol, 
               lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote)) 
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  } else plot.new()
  if (!is.null(main)) 
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    cexO <- par("cex")
    cex1=0.75*keysize/27
    par(mar = c(3, 4.7, 2.8, 0.5), cex = cex1)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    } else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    par("cex"=cexO)
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv, mgp=c(3,0.6,0)*keysize/27,cex.axis=1.2*cexCol*27/keysize)
    if (scale == "row") 
      mtext(side = 1, "Row Z-Score", line = 1.7, cex=cex1*27/keysize)
    else if (scale == "column") 
      mtext(side = 1, "Column Z-Score", line = 1.7, cex=cex1*27/keysize)
    else mtext(side = 1, keyName, line = 1.7, cex=cex1*27/keysize)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y),las=1)
      title("Color Key\nand Density Plot",cex=cex1*27/keysize)
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
      par(cex=cexO)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
            col = denscol)
      axis(2, at = pretty(hy,3)/max(hy) * 0.95, pretty(hy,3),las=1,mgp=c(3,0.9,0)*keysize/27)
      mtext("Color Key\nand Histogram",3,line=0.3,cex=cexCol*keysize/27,adj=0,font=2)
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 3.5,cex=cex1)
      par(cex=cexO)
    }
    else mtext("Color Key",3,line=0.3,cex=cex1*27/keysize,adj=0,font=2)
    if (trace %in% c("both", "column")) {
      vline.vals <- scale01(vline, min.raw, max.raw)
      if (!is.null(vline)) {
        abline(v = vline.vals, col = linecol, lty = 2)
      }
    }
    if (trace %in% c("both", "row")) {
      hline.vals <- scale01(hline, min.raw, max.raw)
      if (!is.null(hline)) {
        abline(v = hline.vals, col = linecol, lty = 2)
      }
    }
  }else plot.new()
  if(!nw){
    enrPW <- annaEnrichPath(data.frame("koID"=labRow,stringsAsFactors=F),bgList,pw=PW,pwn=pwn,pval=pval)
    enrPW <- enrPW[as.numeric(enrPW$numberInSig) >= pathcut,]
  }else{
    enrPW <- annaEnrichPathNW(labRow,pval=pval)
    enrPW <- enrPW[as.numeric(enrPW$numberInSig) >= pathcut,]
  }
  if(nrow(enrPW)>8){
    enrPW <- enrPW[order(enrPW$adjpval,decreasing=F)[1:8],]
    warning("enriched pathways truncated to 8")
  }
  enrPW <- enrPW[order(sapply(enrPW$KOs,function(x) median(which(labRow %in% unlist(strsplit(x,split=";"))))),decreasing=F),]
  print(enrPW)
  par(mai=imai)
  plot(c(0,1),c(1,nr),type="n",yaxs='i',ann=F,axes=F)
  heiwi <- par("pin")/par("din")
  tohei <- par("fig")[4]
  pathposyA <- seq(from=1,to=nr,length.out=length(enrPW$pwDesc)+2)
  text(0.3,pathposyA,label=c("",enrPW$pwDesc,""),adj=0,cex=cexCol)
  pathposy <- pathposyA[-c(1,length(pathposyA))]
  xsca <- par("usr")[1:2]
  ysca <- par("usr")[3:4]
  pushViewport(viewport(xscale=xsca, yscale=ysca,x=1,width=unit(heiwi[1],"npc"),y=tohei,just=c("right","top"),
                        height=unit(heiwi[2],"npc")))
  for(i in 1:nrow(enrPW)){
    concol <- alpha(brewer.pal(nrow(enrPW),"Dark2")[i],0.5)
    for(ko in unlist(strsplit(enrPW$KOs[i],split=";"))){
      lry <- which(labRow==ko)
      grid.bezier(x=c(xsca[1],0.2,0.1,0.25),y=c(lry,lry,pathposy[i],pathposy[i]),default.units="native",gp=gpar(col=concol,lwd=0.7))
    }
  }
  popViewport()  
  
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
                                  high = retval$breaks[-1], color = retval$col)
  if (!is.null(extrafun)) 
    extrafun()
  invisible(retval)
}

annaPrcompBiplot <- function(x, fac, choices = 1L:2L, scale = 1, pc.biplot=FALSE,  var.axes = TRUE, col="black", pch=20,
                             cex = rep(par("cex"), 2), xlabs = NULL, ylabs = NULL, expand=0.7, xlim = NULL, ylim = NULL,
                             arrowLen = 0.1, arrowNum=5,main = NULL, sub = NULL, xlab = NULL, ylab = NULL, font.loading=3, 
                             mar=c(1.5,1.9,0.2,0.2),clabel=0)
{
  if(length(choices) != 2L) stop("length of choices must be 2")
  if(!length(scores <- x$x)) stop(gettextf("object '%s' has no scores", deparse(substitute(x))),domain = NA)
  if(is.complex(scores)) stop("biplots are not defined for complex PCA")
  lam <- x$sdev[choices]
  n <- NROW(scores)
  lam <- lam * sqrt(n)
  if(scale < 0 || scale > 1) warning("'scale' is outside [0, 1]")
  if(scale != 0) lam <- lam^scale else lam <- 1
  if(pc.biplot) lam <- lam / sqrt(n)
  orx <- x
  x <- t(t(scores[, choices]) / lam)
  y <- t(t(orx$rotation[, choices]) * lam)
  n <- nrow(x)
  p <- nrow(y)
  if(missing(xlabs)) {
    xlabs <- dimnames(x)[[1L]]
    if(is.null(xlabs)) xlabs <- 1L:n
  }
  xlabs <- as.character(xlabs)
  dimnames(x) <- list(xlabs, dimnames(x)[[2L]])
  if(missing(ylabs)) {
    ylabs <- dimnames(y)[[1L]]
    if(is.null(ylabs)) ylabs <- paste("Var", 1L:p)
  }
  ylabs <- as.character(ylabs)
  dimnames(y) <- list(ylabs, dimnames(y)[[2L]])
  
  if(length(cex) == 1L) cex <- c(cex, cex)
  
  on.exit(par(op))
  op <- par(pty = "s")
  
  labels = c(paste("PC 1 (",sprintf("%.1f",100*summary(orx)$importance[2,choices[1]]),"%)",sep=""),
             paste("PC 2 (",sprintf("%.1f",100*summary(orx)$importance[2,choices[2]]),"%)",sep=""))
  
  annaSclass(orx$x[, choices], as.factor(fac),col=col, pch=pch,
             clabel=clabel, grid=0,mar=mar,lwd=3,xlab=labels[1],
             ylab=labels[2])
  a <- axis(1,line=10,ann=F,outer=T)
  b <- axis(2,line=10,ann=F,outer=T)
  xfac <- (max(a)-min(a))*expand/(max(y[,1])-min(y[,1]))
  yfac <- (max(b)-min(b))*expand/(max(y[,2])-min(y[,2]))
  arrowOrd <- order(apply(y,1,function(x) x[1]^2+x[2]^2),decreasing=T)[1:arrowNum]
  par(mar=mar)
  text(y[arrowOrd,1L]*xfac,
       y[arrowOrd,2L]*yfac, 
       labels=ylabs[arrowOrd], cex = cex[2L]*0.8, col = "red",
       font=font.loading,mar=c(1.5,1.9,0.2,0.2))
  arrows(0, 0, y[arrowOrd,1L] * 0.8 * xfac, 
         y[arrowOrd,2L] * 0.8 * yfac, col = "red", 
         length=arrowLen)
}

violinplot <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
                        horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
                        lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
                        at, add = FALSE, wex = 1, drawRect = TRUE, las=2 ,ann="xy") 
{
  datas <- list(x, ...)
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  if(length(col)!=n) col <- rep(col,length.out=n)
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h))) 
    args <- c(args, h = h)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      if(grepl("y",ann)) axis(2,las=las)
      if(grepl("x",ann))  axis(1, at = at, label = label,las=las)
    }
    box()
    for (i in 1:n) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col[i], border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
             q3[i], col = rectCol)
        points(at[i], med[i], pch = pchMed, col = colMed)
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), col = col[i], border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}

annaMat2Ind <- function(mat,individuals=indiv,round=F){
  matI <- matrix(0,nrow=nrow(mat),ncol=length(individuals),dimnames=list(rownames(mat),individuals))
  for(i in individuals){
    hits <- grep(i,colnames(mat))
    if(length(hits)>1) matI[,i] <- apply(mat[,hits],1,median) else matI[,i] <- mat[,hits]
    if(round) matI[,i] <- round(matI[,i])
  }
  matI <- matI[which(rowSums(matI)>0),]
  return(matI)
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# barplots
annaMustBarplot <- function(mat,name="barplot",xlab=colnames(mat),ylab="abundance [%]",col="",type="pdf",
                            space=Mspace,silent=TRUE,font.leg=1,sumNorm=T,log="",boco="black"){
  if(col==""){
    if(nrow(mat)%%2==1) {
      col <- colorRampPalette(brewer.pal(11,"Spectral"))(nrow(mat))[c(nrow(mat),c(rbind(1:(nrow(mat)/2),((nrow(mat)/2)+0.5):(nrow(mat)-1))))]
    }else {
      col <- colorRampPalette(brewer.pal(11,"Spectral"))(nrow(mat))[c(rbind(1:(nrow(mat)/2),((nrow(mat)/2)+1):nrow(mat)))]
    }
  }
  plotmat <- mat[order(rowSums(mat),decreasing=T),]
  if(sumNorm) plotmat <- decostand(plotmat,method="total",MARGIN=2)*100
  if(!silent){
    layout(matrix(c(1,2),nrow=1,ncol=2),c(0.8,0.2))
    par(mar=c(4.7,3.2,0.4,0.2),mgp=c(2.2,0.5,0))
    barplot(plotmat,las=2,names.arg=xlab,col=col,ylab=ylab,space=space,log=log,border=boco)
    par(mar=c(0.4,0.2,0.4,0.2))
    plot(1,1,ann=F,axes=F,type="n")
    legend("topleft",legend=rownames(plotmat),fill=col,bty="n",cex=0.8,xpd=T,
           y.intersp=0.75,text.font=font.leg)
  }
  if(type=="pdf"){
    pdf(paste(gsub("-","",sub("20","",Sys.Date())),"-",name,".pdf",sep=""),width=7,height=2.5,pointsize=8)
    layout(matrix(c(1,2),nrow=1,ncol=2),c(0.8,0.2))
    par(mar=c(4.7,3.2,0.4,0.2),mgp=c(2.2,0.5,0))
    barplot(plotmat,las=2,names.arg=xlab,col=col,ylab=ylab,space=space,log=log,border=boco)
    par(mar=c(0.4,0.2,0.4,0.2))
    plot(1,1,ann=F,axes=F,type="n")
    legend("topleft",legend=rownames(plotmat),fill=col,bty="n",cex=0.8,xpd=T,
           y.intersp=0.75,text.font=font.leg)
    dev.off()
  }else{
    tiff(paste(gsub("-","",sub("20","",Sys.Date())),"-",name,".tif",sep=""),width=7,height=2.5,pointsize=8,
         units="in",res=300)
    layout(matrix(c(1,2),nrow=1,ncol=2),c(0.8,0.2))
    par(mar=c(4.7,3.2,0.4,0.2),mgp=c(2.2,0.5,0))
    barplot(plotmat,las=2,names.arg=xlab,col=col,ylab=ylab,space=space,log=log,border=boco)
    par(mar=c(0.4,0.2,0.4,0.2))
    plot(1,1,ann=F,axes=F,type="n")
    legend("topleft",legend=rownames(plotmat),fill=col,bty="n",cex=0.8,xpd=T,
           y.intersp=0.75,text.font=font.leg)
    dev.off()
  }
}
# scatterplots
annaMustPlotI <- function(x,y,name="plot",xlab="",ylab="",col=Mcol,pch=Mchar,type="pdf",legend=Mleg,legpch=Mchar,
                          legcol=Mcol,main="",fmain=1,cmain=0.8,silent=TRUE,secy="",font.secy=1,font.leg=1,log="",cex=1){
  xlim <- c(0.9*min(x),1.1*max(x))
  ylim <- c(0.9*min(y),1.1*max(y))
  legx <- 1.14*max(x)
  legy <- 1.06*max(y)
  if(grepl("x",log)){
    xlim[1] <- min(x[x>0])
  }
  if(grepl("y",log)){
    ylim[1] <- min(y[y>0])
  }
  if(!silent){
    par(mar=c(3.4,4.5,0.6,4),mgp=c(3.5,0.6,0))
    plot(x,y,col=col,pch=pch,xlab="",ylab=ylab,las=1,bty="l",xlim=xlim,ylim=ylim,
         main=main,font.main=fmain,cex.main=cmain,log=log,cex=cex)
    mtext(xlab,1,2)
    if(length(secy)>0) mtext(strtrim(secy,40),2,2.5,font=font.secy)
    legend(legx,legy,legend=legend,bty="n",pch=legpch,col=legcol,cex=0.8,xpd=T,y.intersp=0.8,
           text.font=font.leg)
  }
  if(type=="pdf"){
    pdf(paste(gsub("-","",sub("20","",Sys.Date())),"-",name,".pdf",sep=""),width=3.5,height=2.3,pointsize=8)
    par(mar=c(3.4,4.5,0.6,4),mgp=c(3.5,0.6,0))
    plot(x,y,col=col,pch=pch,xlab="",ylab=ylab,las=1,bty="l",xlim=xlim,ylim=ylim,
         main=main,font.main=fmain,cex.main=cmain,log=log,cex=cex)
    mtext(xlab,1,2)
    if(length(secy)>0) mtext(strtrim(secy,40),2,2.5,font=font.secy)
    legend(legx,legy,legend=legend,bty="n",pch=legpch,col=legcol,cex=0.8,xpd=T,y.intersp=0.8,
           text.font=font.leg)
    dev.off()
  }else{
    tiff(paste(gsub("-","",sub("20","",Sys.Date())),"-",name,".tif",sep=""),width=3.5,height=2.2,pointsize=8,units="in",res=300)
    par(mar=c(3.3,4.5,0.4,3.8),mgp=c(3.5,0.6,0))
    plot(x,y,col=col,pch=pch,xlab="",ylab=ylab,las=1,bty="l",xlim=xlim,ylim=ylim,
         main=main,font.main=fmain,cex.main=cmain,log=1,cex=cex)
    legend(legx,legy,legend=legend,bty="n",pch=legpch,col=legcol,cex=0.8,xpd=T,y.intersp=0.8,
           text.font=font.leg)
    if(length(secy)>0) mtext(strtrim(secy,40),2,2.5,font=font.secy)
    mtext(xlab,1,2)
    dev.off()
  }
}
# of one parameter over visits (with line between points)
annaMustPlotV <- function(vec,visvec=visits,indvec=indiv[visFac],name="plot",xlab="visit",ylab="",col=Mcol,pch=Mchar,type="pdf",legend=Mleg,
                          legpch=Mchar,legcol=Mcol,main="",fmain=1,cmain=0.8,silent=TRUE,secy="",font.secy=1){
  ylim <- c(0.9*min(vec,na.rm=T),1.1*max(vec,na.rm=T))
  mat <- tapply(vec,list(indvec,visvec),sum)
  if(!silent){
    par(mar=c(3.4,4.5,0.6,4),mgp=c(3.5,0.6,0))
    for(i in 1:nrow(mat)){
      plot(mat[i,],ylim=ylim,ann=F,axes=F,type="b",col=col[i],pch=pch[i])
      par(new=T)
    }
    axis(1,at=1:3)
    axis(2,las=1,lwd=0.5)
    mtext(xlab,1,2)
    mtext(ylab,2,3.5)
    if(length(secy)>0) mtext(strtrim(secy,40),2,2.5,font=font.secy)
    legend(3.1,0.96*max(ylim),legend=legend,bty="n",pch=legpch,col=legcol,cex=0.8,xpd=T,y.intersp=0.8)
    box(bty="l")
  }
  if(type=="pdf"){
    pdf(paste(gsub("-","",sub("20","",Sys.Date())),"-",name,".pdf",sep=""),width=3.5,height=2.3,pointsize=8)
    par(mar=c(3.4,4.5,0.6,4),mgp=c(3.5,0.6,0))
    for(i in 1:nrow(mat)){
      plot(mat[i,],ylim=ylim,ann=F,axes=F,type="b",col=col[i],pch=pch[i])
      par(new=T)
    }
    axis(1,at=1:3)
    axis(2,las=1,lwd=0.5)
    mtext(xlab,1,2)
    mtext(ylab,2,3.5)
    if(length(secy)>0) mtext(strtrim(secy,40),2,2.5,font=font.secy)
    legend(3.1,0.96*max(ylim),legend=legend,bty="n",pch=legpch,col=legcol,cex=0.8,xpd=T,y.intersp=0.8)
    box(bty="l")
    dev.off()
  }else{
    tiff(paste(gsub("-","",sub("20","",Sys.Date())),"-",name,".tif",sep=""),width=3.5,height=2.2,pointsize=8,units="in",res=300)
    par(mar=c(3.4,4.5,0.6,4),mgp=c(3.5,0.6,0))
    for(i in 1:nrow(mat)){
      plot(mat[i,],ylim=ylim,ann=F,axes=F,type="b",col=col[i],pch=pch[i])
      par(new=T)
    }
    axis(1,at=1:3)
    axis(2,las=1,lwd=0.5)
    mtext(xlab,1,2)
    mtext(ylab,2,3.5)
    if(length(secy)>0) mtext(strtrim(secy,40),2,2.5,font=font.secy)
    legend(3.1,0.96*max(ylim),legend=legend,bty="n",pch=legpch,col=legcol,cex=0.8,xpd=T,y.intersp=0.8)
    box(bty="l")
    dev.off()
  }
}
# rank-abundance barplots
annaMustRAplotMetaT <- function(mat,name="RAplot",xlab=colnames(mat),ylab="abundance [%]",col="",type="pdf",maxlen=0,lev="other",silent=TRUE) {
  dir.create(paste("./",paste(gsub("-","",sub("20","",Sys.Date())),name,sep="_"),sep=""))
  setwd(paste("./",paste(gsub("-","",sub("20","",Sys.Date())),name,sep="_"),sep=""))
  if(maxlen==0) maxlen <- nrow(mat)
  mat <- mat[order(rowSums(mat),decreasing=T)[1:maxlen],]
  if(col=="") col <- colorRampPalette(brewer.pal(11,"Spectral"))(nrow(mat))
  if(lev %in% c("species","genus")) font <- 3 else font <- 1
  
  pdf(paste(gsub("-","",sub("20","",Sys.Date())),name,"BySample.pdf",sep="-"),width=3.5,height=3.5,pointsize=8)
  par(mar=c(6,3.2,1.4,0.4),mgp=c(2.2,0.6,0))
  for(i in 1:ncol(mat)){
    o2 <- order(mat[,i],decreasing=T)
    b <- barplot(decostand(mat[o2,i],method="total",MARGIN=2)*100,ylim=c(0,100),beside=T,las=2,names.arg="",
                 col=col[o2],ylab=ylab,main=xlab[i])
    text(b, -0.1, labels = rownames(mat)[o2], srt = 45, adj = c(1,1), xpd = TRUE, cex=.75,font=font)
  }
  dev.off()
  
  heigh <- decostand(apply(mat,1,mean),method="total",MARGIN=2)*100
  pdf(paste(gsub("-","",sub("20","",Sys.Date())),name,"Ave.pdf",sep="-"),width=3.5,height=3.5,pointsize=8)
  par(mar=c(6,3.2,0.4,0.4),mgp=c(2.2,0.6,0))
  b <- barplot2(heigh,ylim=c(0,100),beside=T,las=2,names.arg="",col=col,ylab=ylab,plot.ci=T,ci.l=heigh,
                ci.u=heigh+apply(decostand(mat,method="total",MARGIN=2),1,sd)*100)
  text(b, -0.1, labels = rownames(mat)[o2], srt = 45, adj = c(1,1), xpd = TRUE, cex=.75,font=font)
  dev.off()
  
  pdf(paste(gsub("-","",sub("20","",Sys.Date())),name,"All.pdf",sep="-"),width=7,height=3.5,pointsize=8)
  par(mar=c(6,3.2,0.4,0.4),mgp=c(2.2,0.6,0))
  b <- barplot(t(decostand(mat,method="total",MARGIN=2)*100),ylim=c(0,100),beside=T,density=40,
               angle=c(0,135,45)[as.numeric(gsub(".+V","",colnames(mat)))],las=2,names.arg=rep("",nrow(mat)),col=Mcol[visFac],ylab=ylab)
  text(apply(b,2,mean), -0.1, labels = rownames(mat)[o2], srt = 45, adj = c(1,1), xpd = TRUE, cex=.75,font=font)
  legend("topright",legend=xlab,density=40,angle=c(0,135,45)[as.numeric(gsub(".+V","",colnames(mat)))],fill=Mcol[visFac],bty="n",cex=0.8,xpd=T)
  dev.off()
  setwd("../")
  if(!silent){
    par(mar=c(6,3.2,0.4,0.4),mgp=c(2.2,0.6,0))
    b <- barplot2(heigh,ylim=c(0,100),beside=T,las=2,names.arg="",col=col,ylab=ylab,plot.ci=T,ci.l=heigh,
                  ci.u=heigh+apply(decostand(mat,method="total",MARGIN=2),1,sd)*100)
    text(b, -0.1, labels = rownames(mat)[o2], srt = 45, adj = c(1,1), xpd = TRUE, cex=.75,font=font)
  }
}

#####
#boxplot
annaMustPlotX <- function(y,x=miMeta$DIABETESTY1,name="boxplot",xlab="Diabetes type 1",ylab="",type="pdf",main="",pval=""
                          ,pch=Mchar[visFac],col=Mcol[visFac],secy="",font.secy=1,silent=TRUE){
  miny <- 0
  maxy <- 0
  if(min(y)>=0) miny <- min(y)*0.9 else miny <- min(y)*1.1
  if(max(y)>=0) maxy <- max(y)*1.1 else maxy <- max(y)*0.9
  if(!silent){
    par(mar=c(3.3,4.5,0.3,0.2),mgp=c(3.5,0.6,0))
    if(pval !="") par(mar=c(4.3,par("mar")[2:4]))
    if(main !="") par(mar=c(par("mar")[1:2],2,par("mar")[4]))
    boxplot(y ~ x, ylab=ylab,bty="l",las=1,outline=F,ylim=c(miny,maxy))
    points(jitter(as.numeric(as.factor(x)),0.5),y,pch=pch,col=col,cex=0.8)
    mtext(xlab,1,2)
    mtext(pval,1,3)
    if(length(secy)>0) mtext(strtrim(secy,40),2,2.5,font=font.secy)
    mtext(main,3,1,cex=0.6)
  } 
  if(type=="pdf"){
    pdf(paste(gsub("-","",sub("20","",Sys.Date())),"-",name,".pdf",sep=""),width=1.75,height=2.5,pointsize=8)
    par(mar=c(3.3,4.5,0.3,0.2),mgp=c(3.5,0.6,0))
    if(pval !="") par(mar=c(4.3,par("mar")[2:4]))
    if(main !="") par(mar=c(par("mar")[1:2],2,par("mar")[4]))
    boxplot(y ~ x, ylab=ylab,bty="l",las=1,outline=F,ylim=c(miny,maxy))
    points(jitter(as.numeric(as.factor(x)),0.5),y,pch=pch,col=col,cex=0.8)
    mtext(xlab,1,2)
    mtext(pval,1,3)
    if(length(secy)>0) mtext(strtrim(secy,40),2,2.5,font=font.secy)
    mtext(main,3,1,cex=0.6)
    dev.off()
  }else{
    tiff(paste(gsub("-","",sub("20","",Sys.Date())),"-",name,".tif",sep=""),width=1.75,height=2.5,pointsize=8,units="in",
         res=300)
    par(mar=c(3.3,4.5,0.3,0.2),mgp=c(3.5,0.6,0))
    if(pval !="") par(mar=c(4.3,par("mar")[2:4]))
    if(main !="") par(mar=c(par("mar")[1:2],2,par("mar")[4]))
    boxplot(y ~ x,ylab=ylab,bty="l",las=1,outline=F,ylim=c(miny,maxy))
    points(jitter(as.numeric(as.factor(x)),0.3),y,pch=pch,col=col,cex=0.8)
    mtext(xlab,1,2)
    mtext(pval,1,3)
    if(length(secy)>0) mtext(strtrim(secy,40),2,2.5,font=font.secy)
    mtext(main,3,1,cex=0.6)
    dev.off()
  }
}
# Wicoxon rank sum test and boxplot of significant observations (with multiple testing adjustment)
annaMustWilcPlotMT <- function(mat,x,grpName,sigMT,outName,pch=Mchar[visFac],col=Mcol[visFac],unit="[%]",
                               secy=matrix("",ncol=2,nrow=1),font.secy=1,silent=TRUE,res=F){
  wilcVec <- vector(length=nrow(mat))
  names(wilcVec) <- rownames(mat)
  wilcVecE <- vector(length=nrow(mat))
  names(wilcVecE) <- rownames(mat)
  for(i in 1:nrow(mat)){
    if(length(unique(mat[i,x]))>1){
      tmpw <- wilcox.test(as.matrix(mat)[i,x],as.matrix(mat)[i,!x],exact=F,conf.int = T)
      wilcVec[i] <- tmpw$p.value
      wilcVecE[i] <- tmpw$estimate
    }else{
      tmpw <- wilcox.test(as.matrix(mat)[i,x],as.matrix(mat)[i,!x],exact=F)
      wilcVec[i] <- tmpw$p.value
      wilcVecE[i] <- 0
    }
  }
  wilcVec[is.na(wilcVec)] <- 1
  wilcVecE[is.na(wilcVec)] <- 0
  if(length(wilcVec[p.adjust(wilcVec,"fdr")<=sigMT])>0){
    mt <- data.frame("group1median","greaterMedian",,"diff"=wilcVecE[p.adjust(wilcVec,"fdr")<=sigMT],"pval"=wilcVec[p.adjust(wilcVec,"fdr")<=sigMT],
                     "adjpval"=p.adjust(wilcVec,"fdr")[p.adjust(wilcVec,"fdr")<=sigMT],stringsAsFactors=F)
    colnames(mt) <- c("group1median","greaterMedian","diff","pval","adjpval")
    k=1
    for(i in which(p.adjust(wilcVec,"fdr")<=sigMT)){
      mt[k,2] <- grpName[1+as.numeric(median(as.matrix(mat)[i,x])>median(as.matrix(mat)[i,!x]))]
      mt[k,1] <- median(as.matrix(mat)[i,!x])
      k=k+1   
    }
    write.table(mt,paste(gsub("-","",sub("20","",Sys.Date())),"-",outName,"MTlt",sigMT,".tsv",sep=""),sep="\t")
    
    system(paste("mkdir ",gsub("-","",sub("20","",Sys.Date())),"_",outName,"MTlt",sigMT,"_Boxplots",sep=""))
    if(unit=="[%]") perc <- 100 else perc <- 1
    setwd(paste(gsub("-","",sub("20","",Sys.Date())),"_",outName,"MTlt",sigMT,"_Boxplots",sep=""))
    for(i in which(p.adjust(wilcVec,"fdr")<=sigMT)) {
      annaMustPlotX(as.matrix(mat)[i,]*perc,x,name=gsub("/","",paste("boxplot",rownames(mat)[i],grpName[length(grpName)],
                                                                     sep="")),
                    ylab=paste(rownames(mat)[i], unit),xlab=grpName[2],
                    pval=paste("adj p-val =",round(p.adjust(wilcVec,"fdr")[i],3)),pch=pch,col=col,
                    secy=secy[rownames(mat)[i]==secy[,1],2],font.secy=font.secy,silent=silent)
    }
    setwd("../")
    if(res) rownames(mat)[wilcVec<=sig]
  }else warning("No significant observations!")
}
# Wicoxon rank sum test and boxplot of significant observations (without multiple testing adjustment)
annaMustWilcPlotWOMT <- function(mat,x,grpName,sig,outName,pch=Mchar[visFac],col=Mcol[visFac],unit="[%]",
                                 secy=matrix("",ncol=2,nrow=1),font.secy=1,silent=TRUE,res=F){
  wilcVec <- vector(length=nrow(mat))
  names(wilcVec) <- rownames(mat)
  wilcVecE <- vector(length=nrow(mat))
  names(wilcVecE) <- rownames(mat)
  for(i in 1:nrow(mat)){
    if(length(unique(mat[i,x]))>1){
      tmpw <- wilcox.test(as.matrix(mat)[i,x],as.matrix(mat)[i,!x],exact=F,conf.int = T)
      wilcVec[i] <- tmpw$p.value
      wilcVecE[i] <- tmpw$estimate
    }else{
      tmpw <- wilcox.test(as.matrix(mat)[i,x],as.matrix(mat)[i,!x],exact=F)
      wilcVec[i] <- tmpw$p.value
      wilcVecE[i] <- 0
    }
  }
  wilcVec[is.na(wilcVec)] <- 1
  wilcVecE[is.na(wilcVec)] <- 0
  if(length(wilcVec[wilcVec<=sig])>0){
    mt <- data.frame("group1median","greaterMedian","diff"=wilcVecE[wilcVec<=sig],"pval"=wilcVec[wilcVec<=sig],
                     "adjpval"=p.adjust(wilcVec,"fdr")[wilcVec<=sig],stringsAsFactors=F)
    colnames(mt) <- c("group1median","greaterMedian","diff","pval","adjpval")
    k=1
    for(i in which(wilcVec<=sig)){
      mt[k,2] <- grpName[1+as.numeric(median(as.matrix(mat)[i,x])>median(as.matrix(mat)[i,!x]))]
      mt[k,1] <- median(as.matrix(mat)[i,!x])
      k=k+1   
    }
    write.table(mt,paste(gsub("-","",sub("20","",Sys.Date())),"-",outName,"pvallt",sig,".tsv",sep=""),sep="\t")
    
    system(paste("mkdir ",gsub("-","",sub("20","",Sys.Date())),"_",outName,"pvallt",sig,"_Boxplots",sep=""))
    if(unit=="[%]") perc <- 100 else perc <- 1
    setwd(paste(gsub("-","",sub("20","",Sys.Date())),"_",outName,"pvallt",sig,"_Boxplots",sep=""))
    for(i in which(wilcVec<=sig)) {
      annaMustPlotX(as.matrix(mat)[i,]*perc,x,name=gsub("/","",paste("boxplot",rownames(mat)[i],grpName[length(grpName)],
                                                                     sep="")),
                    ylab=paste(rownames(mat)[i], unit),xlab=grpName[2],font.secy=font.secy,
                    pval=paste("p-val =",round(wilcVec[i],3)),pch=pch,col=col,secy=secy[rownames(mat)[i]==secy[,1],2],silent=silent)
    }
    setwd("../")
    if(res) rownames(mat)[wilcVec<=sig]
  }else warning("No significant observations!")
}

# dataframe to matrix
annaSdf2Mat <- function(df){
  mat <- as.matrix(df[,2:ncol(df)])
  row.names(mat) <- df[,1]
  return(mat)
}
#### pathway enrichment analysis ####
annaEnrichPath <- function(koTab,bgList,pw=PW,pwn=pwn,koIDname="koID",pval=0.05){
  if(class(koTab)=="data.frame") sigPW <- merge(koTab,pw,by.x=koIDname,by.y="ko")
  if(class(koTab)=="matrix") sigPW <- merge(koTab,pw,by.x=0,by.y="ko")
  if(class(koTab)=="character" ) sigPW <- merge(koTab,pw,by.x=1,by.y="ko")
  sigPW <- unique(sigPW)
  colnames(sigPW)[1] <- "koID"
  
  bgSet <- merge(bgList,pw,by.x=1,by.y="ko")
  
  colnames(bgSet)[1] <- "ko" 
  bgSet <- unique(bgSet)
  pwInf <- data.frame("pathway"=unique(sigPW$pw),"numberInSig","KOs","numberTotal","pval",stringsAsFactors=F)
  colnames(pwInf) <- gsub(".","",gsub("X.","",colnames(pwInf)),fixed=T)
  
  for(i in pwInf$pathway){
    q <- length(unique(sigPW$koID[sigPW$pw==i]))
    m <- length(unique(bgSet$ko[bgSet$pw==i]))
    n <- length(unique(bgSet$ko)) - m
    k <- length(unique(sigPW$koID))
    pwInf$numberInSig[pwInf$pathway==i] <- q
    pwInf$KOs[pwInf$pathway==i] <- paste(sigPW$koID[which(sigPW$pw==i)],sep=";",collapse=";")
    pwInf$numberTotal[pwInf$pathway==i] <- m
    pwInf$pval[pwInf$pathway==i] <- phyper(q-1,m,n,k,lower.tail=F) 
  }
  pwInf[,6] <- p.adjust(pwInf[,5],"fdr")
  colnames(pwInf)[6] <- "adjpval"
  
  pwEn <- pwInf[ pwInf$adjpval<=pval,]
  pwEn <- merge(pwEn,pwn,by.x="pathway",by.y="pw",all.x=T)
  return(pwEn)
}

annaEnrichPathNW <- function(koVec,pval){
  sigPW <- merge(unlist(sapply(koVec,function(x) unlist(strsplit(x,split="-")))),pw,by.x=1,by.y="ko")
  sigPW <- unique(sigPW)
  colnames(sigPW)[1] <- "koID"
  bgSet <- merge(unlist(sapply(nwNodes$V1,function(x) unlist(strsplit(x,split="-")))),pw,by.x=1,by.y="ko")
  colnames(bgSet)[1] <- "ko" 
  
  bgSet <- unique(bgSet)
  pwInf <- data.frame("pathway"=unique(sigPW$pw),"numberInSig","KOs","numberTotal","pval",stringsAsFactors=F)
  colnames(pwInf) <- gsub(".","",gsub("X.","",colnames(pwInf)),fixed=T)
  
  for(i in pwInf$pathway){
    q <- length(unique(sigPW$koID[sigPW$pw==i]))
    m <- length(unique(bgSet$ko[bgSet$pw==i]))
    n <- length(unique(bgSet$ko)) - m
    k <- length(unique(sigPW$koID))
    pwInf$numberInSig[pwInf$pathway==i] <- q
    pwInf$KOs[pwInf$pathway==i] <- paste(sigPW$koID[which(sigPW$pw==i)],sep=";",collapse=";")
    pwInf$numberTotal[pwInf$pathway==i] <- m
    pwInf$pval[pwInf$pathway==i] <- phyper(q-1,m,n,k,lower.tail=F) 
  }
  pwInf[,6] <- p.adjust(pwInf[,5],"fdr")
  colnames(pwInf)[6] <- "adjpval"
  
  pwEn <- pwInf[ pwInf$adjpval<=pval,]
  pwEn <- merge(pwEn,pwn,by.x="pathway",by.y="pw",all.x=T)
  return(pwEn)
}

#PCO analysis + plot
#####
JSD.pair <- function (x, y) { # from phyloseq
  u <- x/sum(x)
  v <- y/sum(y)
  m <- (u + v)/2
  if (all(u * v > 0)) {
    d <- (u * log(u/m) + v * log(v/m))/2
  }
  else {
    P1 <- u * log(u/m)
    P2 <- v * log(v/m)
    P1[is.nan(P1)] <- 0
    P2[is.nan(P2)] <- 0
    d <- (P1 + P2)/2
  }
  return(sum(d))
}
JSD <- function(mat){ # from phyloseq
  spn <- combn(colnames(mat), 2, simplify=FALSE)
  
  # initialize DistMat with NAs
  DistMat <- matrix(NA, ncol(mat), ncol(mat))
  # define the rows/cols of DistMat with the sample names (rownames)    
  rownames(DistMat) <- colnames(mat)
  colnames(DistMat) <- colnames(mat)
  
  ## Format coercion
  # Coerce to the vegan orientation, with species as columns
  OTU <- t(mat)
  # Coerce OTU to matrix for calculations.
  OTU <- as(OTU, "matrix")
  
  # optionally-parallel implementation with foreach
  distlist <- foreach( i = spn) %dopar% {
    A <- i[1]
    B <- i[2]
    return( JSD.pair(OTU[A, ], OTU[B, ]) )
  }
  # return(distlist)
  # This is in serial, but it is quick.
  distlist2distmat <- function(i, spn, DL){
    DistMat[ spn[[i]][2], spn[[i]][1] ] <<- DL[[i]]
  }
  junk <- sapply(1:length(spn), distlist2distmat, spn, distlist)
  
  return(as.dist(DistMat))
}


pcoaSuper <- function(mat,title,Mcol,Mchar,visFac,indiv,distMeth="bray",bina=F,silent=TRUE,ax=1:2){
  mBC <- vegdist(t(mat),distMeth,bina)
  mPCOA <- pcoa(mBC)
  mPCOAval <- 100*c(mPCOA$values[(1:3),2]/sum(mPCOA$values[mPCOA$values[,2]>0,2]))
  pdf(paste(gsub("-","",sub("20","",Sys.Date())),title,"_PCoA.pdf",sep=""),width=4.2,height=3.5,pointsize=8)
  par(mar=c(1.5,1.9,0.2,0.2))
  layout(matrix(c(1,2),1,2,byrow=T),c(5/6,1/6))
  annaSclass(mPCOA$vector[,ax], as.factor(visFac),col=Mcol[as.numeric(levels(as.factor((visFac))))], pch=Mchar[visFac],
             clabel=0, grid=0,mar=c(1.5,1.9,0.2,0.2),lwd=3,xlab=paste("PC ",ax[1]," (",sprintf("%.1f",mPCOAval[ax[1]]),"%)",sep=""),
             ylab=paste("PC ",ax[2]," (",sprintf("%.1f",mPCOAval[ax[2]]),"%)",sep=""))
  par(mar=c(0.2,0.2,0.2,0.2))
  plot(1,1,ann=F,type="n",axes=F)
  legend("topleft",legend=indiv,bty="n",pch=Mchar[as.numeric(levels(as.factor((visFac))))],
         col=Mcol[as.numeric(levels(as.factor((visFac))))],cex=1,xpd=T,y.intersp=1)
  dev.off()
  pdf(paste(gsub("-","",sub("20","",Sys.Date())),title,"_PCoA_3dplot.pdf",sep=""),width=4,height=4,pointsize=8)
  par(mar=c(2,3.2,0.4,1),mgp=c(2,0.6,0))
  scatterplot3d(mPCOA$vectors[,1:3],color=Mcol[visFac],pch=Mchar[visFac],main="",
                xlab=paste("PC 1 (",sprintf("%.1f",mPCOAval[1]),"%)",sep=""),
                ylab=paste("PC 2 (",sprintf("%.1f",mPCOAval[2]),"%)",sep=""),
                zlab=paste("PC 3 (",sprintf("%.1f",mPCOAval[3]),"%)",sep=""))
  dev.off() 
  
  jsdd <-JSD(mat)
  mPCOAJS <- pcoa(jsdd)
  mPCOAvalJS <- 100*c(mPCOAJS$values[(1:3),2]/sum(mPCOAJS$values[mPCOAJS$values[,2]>0,2]))
  pdf(paste(gsub("-","",sub("20","",Sys.Date())),title,"_PCoA_JS.pdf",sep=""),width=4.2,height=3.5,pointsize=8)
  par(mar=c(1.5,1.9,0.2,0.2))
  layout(matrix(c(1,2),1,2,byrow=T),c(5/6,1/6))
  annaSclass(mPCOAJS$vector[,ax], as.factor(visFac),col=Mcol[as.numeric(levels(as.factor((visFac))))], pch=Mchar[visFac],
             clabel=0, grid=0,mar=c(1.5,1.9,0.2,0.2),lwd=3,xlab=paste("PC ",ax[1]," (",sprintf("%.1f",mPCOAvalJS[ax[1]]),"%)",sep=""),
             ylab=paste("PC ",ax[2]," (",sprintf("%.1f",mPCOAvalJS[ax[2]]),"%)",sep=""))
  par(mar=c(0.2,0.2,0.2,0.2))
  plot(1,1,ann=F,type="n",axes=F)
  legend("topleft",legend=indiv,bty="n",pch=Mchar[as.numeric(levels(as.factor((visFac))))],
         col=Mcol[as.numeric(levels(as.factor((visFac))))],cex=1,xpd=T,y.intersp=1)
  dev.off()
  pdf(paste(gsub("-","",sub("20","",Sys.Date())),title,"_PCoA_JS_3dplot.pdf",sep=""),width=4,height=4,pointsize=8)
  par(mar=c(2,3.2,0.4,1),mgp=c(2,0.6,0))
  scatterplot3d(mPCOAJS$vectors[,1:3],color=Mcol[visFac],pch=Mchar[visFac],main="",
                xlab=paste("PC 1 (",sprintf("%.1f",mPCOAvalJS[1]),"%)",sep=""),
                ylab=paste("PC 2 (",sprintf("%.1f",mPCOAvalJS[2]),"%)",sep=""),
                zlab=paste("PC 3 (",sprintf("%.1f",mPCOAvalJS[3]),"%)",sep=""))
  dev.off() 
  
  #mScale <- cenLR(min(mat[mat>0])+t(mat))$x.clr
  mScale <- clr(min(mat[mat>0])+t(mat))
  mPCA <- prcomp(mScale)
  mPCAval <- 100 * summary(mPCA)$importance[2,1:max(ax)]
  pdf(paste(gsub("-","",sub("20","",Sys.Date())),title,"_PCA_3dplot.pdf",sep=""),width=4,height=4,pointsize=8)
  par(mar=c(2,3.2,0.4,1),mgp=c(2,0.6,0))
  scatterplot3d(mPCA$x[,1:3],color=Mcol[visFac],pch=Mchar[visFac],main="",
                xlab=paste("PC 1 (",sprintf("%.1f",mPCAval[1]),"%)",sep=""),
                ylab=paste("PC 2 (",sprintf("%.1f",mPCAval[2]),"%)",sep=""),
                zlab=paste("PC 3 (",sprintf("%.1f",mPCAval[3]),"%)",sep=""))
  dev.off() 
  pdf(paste(gsub("-","",sub("20","",Sys.Date())),title,"_PCA_biplot.pdf",sep=""),width=4.2,height=3.5,pointsize=8)
  par(mar=c(1.5,1.9,0.2,0.2))
  layout(matrix(c(1,2),1,2,byrow=T),c(5/6,1/6))
  annaPrcompBiplot(mPCA,fac=as.factor(visFac),col=Mcol[as.numeric(levels(as.factor((visFac))))],pch=Mchar[visFac],choices=ax)
  par(mar=c(0.2,0.2,0.2,0.2))
  plot(1,1,ann=F,type="n",axes=F)
  legend("topleft",legend=indiv,bty="n",pch=Mchar[as.numeric(levels(as.factor((visFac))))],
         col=Mcol[as.numeric(levels(as.factor((visFac))))],cex=1,xpd=T,y.intersp=1)
  dev.off() 
  if(!silent){
    par(mar=c(1.5,1.9,0.2,0.2))
    layout(matrix(c(1,2),1,2,byrow=T),c(5/6,1/6))
    annaSclass(mPCOA$vector[,ax], as.factor(visFac),col=Mcol[as.numeric(levels(as.factor((visFac))))], pch=Mchar[visFac],
               clabel=0, grid=0,mar=c(1.5,1.9,0.2,0.2),lwd=3,xlab=paste("PC ",ax[1]," (",sprintf("%.1f",mPCOAval[ax[1]]),"%)",sep=""),
               ylab=paste("PC ",ax[2]," (",sprintf("%.1f",mPCOAval[ax[2]]),"%)",sep=""))
    par(mar=c(0.2,0.2,0.2,0.2))
    plot(1,1,ann=F,type="n",axes=F)
    legend("topleft",legend=indiv,bty="n",pch=Mchar[as.numeric(levels(as.factor((visFac))))],
           col=Mcol[as.numeric(levels(as.factor((visFac))))],cex=1,xpd=T,y.intersp=1)
    Sys.sleep(0.1)
    par(mar=c(1.5,1.9,0.2,0.2))
    layout(matrix(c(1,2),1,2,byrow=T),c(5/6,1/6))
    annaSclass(mPCOAJS$vector[,ax], as.factor(visFac),col=Mcol[as.numeric(levels(as.factor((visFac))))], pch=Mchar[visFac],
               clabel=0, grid=0,mar=c(1.5,1.9,0.2,0.2),lwd=3,xlab=paste("PC ",ax[1]," (",sprintf("%.1f",mPCOAvalJS[ax[1]]),"%)",sep=""),
               ylab=paste("PC ",ax[2]," (",sprintf("%.1f",mPCOAvalJS[ax[2]]),"%)",sep=""))
    par(mar=c(0.2,0.2,0.2,0.2))
    plot(1,1,ann=F,type="n",axes=F)
    legend("topleft",legend=indiv,bty="n",pch=Mchar[as.numeric(levels(as.factor((visFac))))],
           col=Mcol[as.numeric(levels(as.factor((visFac))))],cex=1,xpd=T,y.intersp=1)
    Sys.sleep(0.1)
    par(mar=c(1.5,1.9,0.2,0.2),mfrow=c(1,1))
    layout(matrix(c(1,2),1,2,byrow=T),c(5/6,1/6))
    par(mar=c(0.2,0.2,0.2,0.2))
    annaPrcompBiplot(mPCA,fac=as.factor(visFac),col=Mcol[as.numeric(levels(as.factor((visFac))))],pch=Mchar[visFac],choices=ax)
    par(mar=c(0.2,0.2,0.2,0.2))
    plot(1,1,ann=F,type="n",axes=F)
    legend("topleft",legend=indiv,bty="n",pch=Mchar[as.numeric(levels(as.factor((visFac))))],
           col=Mcol[as.numeric(levels(as.factor((visFac))))],cex=1,xpd=T,y.intersp=1)
    par(mar=c(3,4,1,0.2),mfrow=c(1,1))
  }
}

pcoaAnyDist <- function(dist,title,Mcol,Mchar,visFac,indiv,silent=TRUE,ax=1:2){
  mPCOA <- pcoa(dist)
  mPCOAval <- 100*c(mPCOA$values[(1:3),2]/sum(mPCOA$values[mPCOA$values[,2]>0,2]))
  pdf(paste(gsub("-","",sub("20","",Sys.Date())),title,"_dist_PCoA.pdf",sep=""),width=4.2,height=3.5,pointsize=8)
  par(mar=c(1.5,1.9,0.2,0.2))
  layout(matrix(c(1,2),1,2,byrow=T),c(5/6,1/6))
  annaSclass(mPCOA$vector[,ax], as.factor(visFac),col=Mcol[as.numeric(levels(as.factor((visFac))))], pch=Mchar[visFac],
             clabel=0, grid=0,mar=c(1.5,1.9,0.2,0.2),lwd=3,xlab=paste("PC ",ax[1]," (",sprintf("%.1f",mPCOAval[ax[1]]),"%)",sep=""),
             ylab=paste("PC ",ax[2]," (",sprintf("%.1f",mPCOAval[ax[2]]),"%)",sep=""))
  par(mar=c(0.2,0.2,0.2,0.2))
  plot(1,1,ann=F,type="n",axes=F)
  legend("topleft",legend=indiv,bty="n",pch=Mchar[as.numeric(levels(as.factor((visFac))))],
         col=Mcol[as.numeric(levels(as.factor((visFac))))],cex=1,xpd=T,y.intersp=1)
  dev.off()
  pdf(paste(gsub("-","",sub("20","",Sys.Date())),title,"_dist_PCoA_3dplot.pdf",sep=""),width=4,height=4,pointsize=8)
  par(mar=c(2,3.2,0.4,1),mgp=c(2,0.6,0))
  scatterplot3d(mPCOA$vectors[,1:3],color=Mcol[visFac],pch=Mchar[visFac],main="",
                xlab=paste("PC 1 (",sprintf("%.1f",mPCOAval[1]),"%)",sep=""),
                ylab=paste("PC 2 (",sprintf("%.1f",mPCOAval[2]),"%)",sep=""),
                zlab=paste("PC 3 (",sprintf("%.1f",mPCOAval[3]),"%)",sep=""))
  dev.off() 
  if(!silent){
    par(mar=c(1.5,1.9,0.2,0.2))
    layout(matrix(c(1,2),1,2,byrow=T),c(5/6,1/6))
    annaSclass(mPCOA$vector[,ax], as.factor(visFac),col=Mcol[as.numeric(levels(as.factor((visFac))))], pch=Mchar[visFac],
               clabel=0, grid=0,mar=c(1.5,1.9,0.2,0.2),lwd=3,xlab=paste("PC ",ax[1]," (",sprintf("%.1f",mPCOAval[ax[1]]),"%)",sep=""),
               ylab=paste("PC ",ax[2]," (",sprintf("%.1f",mPCOAval[ax[2]]),"%)",sep=""))
    par(mar=c(0.2,0.2,0.2,0.2))
    plot(1,1,ann=F,type="n",axes=F)
    legend("topleft",legend=indiv,bty="n",pch=Mchar[as.numeric(levels(as.factor((visFac))))],
           col=Mcol[as.numeric(levels(as.factor((visFac))))],cex=1,xpd=T,y.intersp=1)
  }
}

# Wilcoxon test with and without multiple testing adjustment + plots
wilcoxSuper <- function(mat,x,grpName,sigMT=0.05,sig=0.01,outName,pch=Mchar[visFac],col=Mcol[visFac]
                        ,secy=matrix("",ncol=2,nrow=1),unit="[%]",silent=TRUE,font.secy=1,res=F){
  resMT <- annaMustWilcPlotMT(mat=mat,x=x,grpName=grpName,sigMT=sigMT,outName=outName,pch=pch,col=col,secy=secy,unit=unit,silent=silent,font.secy=font.secy,res=res)
  resWOMT <- annaMustWilcPlotWOMT(mat=mat,x=x,grpName=grpName,sig=sig,outName=outName,pch=pch,col=col,secy=secy,unit=unit,silent=silent,font.secy=font.secy,res=res)
  list(resMT,resWOMT)
}
#makes a group of 2 and the rest: (works for fam 1,3,4 and can be used for fam 2)
wilcShuf <- function(mat,ind1,ind2,famID){
  indi1 <- paste(famID,ind1,sep="-")
  indi2 <- paste(famID,ind2,sep="-")
  wilcRes <- vector(length=nrow(mat))
  names(wilcRes) <- rownames(mat)
  for(i in 1:nrow(mat)){
    wilcRes[i] <- wilcox.test(as.matrix(mat)[i,c(grep(indi1,colnames(mat)),grep(indi2,colnames(mat)))],
                              as.matrix(mat)[i,-c(grep(indi1,colnames(mat)),grep(indi2,colnames(mat)))],exact=F)$p.value
  }
  wilcRes[is.na(wilcRes)] <- 1
  print(c(indi1,indi2))
  print(c("lowest p-Value",min(wilcRes,na.rm=T)))
  print(c("below 0.05",length(which(wilcRes<0.05))))
  print(c("below 0.01",length(which(wilcRes<0.01))))
  print(c("below 0.05 adjusted",length(which(p.adjust(wilcRes,"fdr")<0.05))))
}
#makes a group of 3 and the rest: (can be used for fam 3 and 4)
wilcShuf3 <- function(mat,ind1,ind2,ind3,famID){
  indi1 <- paste(famID,ind1,sep="-")
  indi2 <- paste(famID,ind2,sep="-")
  indi3 <- paste(famID,ind3,sep="-")
  wilcRes <- vector(length=nrow(mat))
  names(wilcRes) <- rownames(mat)
  for(i in 1:nrow(mat)){
    wilcRes[i] <- wilcox.test(as.matrix(mat)[i,c(grep(indi1,colnames(mat)),grep(indi2,colnames(mat)),grep(indi3,colnames(mat)))],
                              as.matrix(mat)[i,-c(grep(indi1,colnames(mat)),grep(indi2,colnames(mat)),grep(indi3,colnames(mat)))],exact=F)$p.value
  }
  wilcRes[is.na(wilcRes)] <- 1
  print(c(indi1,indi2,indi3))
  print(c("lowest p-Value",min(wilcRes,na.rm=T)))
  print(c("below 0.05",length(which(wilcRes<0.05))))
  print(c("below 0.01",length(which(wilcRes<0.01))))
  print(c("below 0.05 adjusted",length(which(p.adjust(wilcRes,"fdr")<0.05))))
}


#DESeq2 wrapper:
annaMustDESeq <- function(countData,colData=miMeta[visFac,],designFor=as.formula("~ DIABETESTY1"),sigOnly=T,all=F){
  dds <- DESeqDataSetFromMatrix(countData=countData,colData=colData,design=designFor)
  dds <- DESeq(dds,quiet=T)
  resDS <- results(dds)
  resDS <- resDS[order(resDS$padj),]
  if(all) return (dds) else if(sigOnly) return(resDS[resDS$padj<= 0.05 & !is.na(resDS$padj),]) else return(resDS)
}
#DESeq2 wrapper that also plots and writes a table with stats:
#usage example: annaMustPlotDESeq(clusCountMR,"testIA2",clusAbundM,colData=data.frame("waist"=(physexByV$ave.WAISTCIRC>=70)[visFac]),designFor=as.formula("~waist"),secy=motuAnn2,font.secy=3,silent=T)
annaMustPlotDESeq <- function(countData,name,plotData,colData=miMeta[visFac,],designFor=as.formula("~ DIABETESTY1"),w2d=T,
                              pch=Mchar[visFac],col=Mcol[visFac],unit="[%]",
                              secy=matrix("",ncol=2,nrow=1),font.secy=3,silent=TRUE,sig=0.05,res=c("count","all","sig","none")){
  ddsRes <- annaMustDESeq(countData,colData=colData,designFor=designFor,sigOnly=F)
  if(w2d) write.table(ddsRes,paste(gsub("-","",sub("20","",Sys.Date())),"_",name,".tsv",sep=""),sep="\t")
  system(paste("mkdir ",gsub("-","",sub("20","",Sys.Date())),"_",name,"DESeq",sig,"_Boxplots",sep=""))
  if(unit=="[%]") perc <- 100 else perc <- 1
  setwd(paste(gsub("-","",sub("20","",Sys.Date())),"_",name,"DESeq",sig,"_Boxplots",sep=""))
  for(i in which(ddsRes$padj<sig & !is.na(ddsRes$padj))){
    plotName <- rownames(ddsRes)[i]
    annaMustPlotX(plotData[rownames(plotData)==plotName,]*perc,colData[,as.character(designFor)[2]],
                  name=gsub("/","",paste("boxplot","_",plotName,"_",name,sep="")),
                  ylab=paste(plotName, unit),xlab=as.character(designFor)[2],
                  pval=paste("adj p-val =",round(ddsRes$padj[i],3)),pch=pch,col=col,
                  secy=secy[plotName==secy[,1],2],font.secy=font.secy,silent=silent)
    Sys.sleep(0.1)
  }
  setwd("../")
  if(res[1] == "count"){
    return(length(ddsRes$padj[ddsRes$padj<sig & !is.na(ddsRes$padj)]))
  }else if(res[1] == "all"){
    return(ddsRes)
  }else if(res[1] == "sig"){
    if(nrow(ddsRes[ddsRes$padj<sig & !is.na(ddsRes$padj),])==0) warning("No significant findings.")
    return(ddsRes[ddsRes$padj<sig & !is.na(ddsRes$padj),])
  }
}

#DESeq2 to shuffle groups of 1 + x within families
annaMustShuf1DESeq <- function(countData,ind1,famID,levels=c(0.01,0.05,0.1),sep="-"){
  famData <- countData[,grep(famID,colnames(countData))]
  indi1 <- paste(famID,ind1,sep=sep)
  colData <- data.frame("indOI"=colnames(famData) %in% grep(indi1,colnames(famData),value=T))
  tab <- annaMustDESeq(famData,colData,designFor=as.formula("~ indOI"),sigOnly=F)
  return(c(tab$padj[1],sapply(levels, function(x) length(tab$padj[tab$padj<x & !is.na(tab$padj)]))))
}
#DESeq2 to shuffle groups of 2 + x within families
annaMustShuf2DESeq <- function(countData,ind1,ind2,famID,levels=c(0.01,0.05,0.1),sep="-"){
  famData <- countData[,grep(famID,colnames(countData))]
  indi1 <- paste(famID,ind1,sep=sep)
  indi2 <- paste(famID,ind2,sep=sep)
  colData <- data.frame("indOI"=colnames(famData) %in% 
                          c(grep(indi1,colnames(famData),value=T),grep(indi2,colnames(famData),value=T)))
  tab <- annaMustDESeq(famData,colData,designFor=as.formula("~ indOI"),sigOnly=F)
  return(c(tab$padj[1],sapply(levels, function(x) length(tab$padj[tab$padj<x & !is.na(tab$padj)]))))
}
#DESeq2 to shuffle groups of 3 + x within families
annaMustShuf3DESeq <- function(countData,ind1,ind2,ind3,famID,levels=c(0.01,0.05,0.1),sep="-"){
  famData <- countData[,grep(famID,colnames(countData))]
  indi1 <- paste(famID,ind1,sep=sep)
  indi2 <- paste(famID,ind2,sep=sep)
  indi3 <- paste(famID,ind3,sep=sep)
  colData <- data.frame("indOI"=colnames(famData) %in% 
                          c(grep(indi1,colnames(famData),value=T),grep(indi2,colnames(famData),value=T),
                            grep(indi3,colnames(famData),value=T)))
  tab <- annaMustDESeq(famData,colData,designFor=as.formula("~ indOI"),sigOnly=F)
  return(c(tab$padj[1],sapply(levels, function(x) length(tab$padj[tab$padj<x & !is.na(tab$padj)]))))
}


#limma/voom wrapper:
annaMustVoom <- function(countData,colData=miMeta[visFac,],designFor=as.formula("~ DIABETESTY1"),block=visFac,coef="DIABETESTY1Yes",
                         sig=0.05,sigOnly=T,all=F){
  y <- DGEList(counts = countData)
  y <- calcNormFactors(y)
  design <- model.matrix(designFor,data=colData)
  v <- voom(y, design)
  if(length(block)) {
    cor <- duplicateCorrelation(v, design, block = block)
    #cor$consensus
    fit <- lmFit(v, design, block = block, correlation = cor$consensus)
  } else {
    fit <- lmFit(v, design)
  }
  fit <- eBayes(fit)
  #summary(decideTests(fit))
  if(all) {
    resVoom <- fit
  } else if(sigOnly) {
    resVoom <- topTable(fit, coef = coef, n = Inf, sort = "p", p = sig)
  } else {
    resVoom <- topTable(fit, coef = coef, n = Inf, sort = "p", p = 1)
  }
  return (resVoom)
}
#limma/voom wrapper that also plots and writes a table with stats:
#usage example: annaMustPlotVoom(clusCountMR,"testIA2",clusAbundM,colData=data.frame("waist"=(physexByV$ave.WAISTCIRC>=70)[visFac]),designFor=as.formula("~waist"),secy=motuAnn2,font.secy=3,silent=T)
annaMustPlotVoom <- function(countData,name,plotData,colData=data.frame(miMetaF[visFac,]),designFor=as.formula("~ DIABETESTY1"),
                             block=visFac,coef="DIABETESTY1Yes",w2d=T,plot=T,
                             pch=Mchar[visFac],col=Mcol[visFac],unit="[%]",
                             secy=matrix("",ncol=2,nrow=1),font.secy=3,silent=TRUE,sig=0.05,res=c("count","all","sig","none")){
  voomRes <- annaMustVoom(countData,colData=colData,designFor=designFor,block=block,coef=coef,sigOnly=F)
  if(w2d) write.table(voomRes,paste(gsub("-","",sub("20","",Sys.Date())),"_",name,".tsv",sep=""),sep="\t")
  if(plot){
    system(paste("mkdir ",gsub("-","",sub("20","",Sys.Date())),"_",name,"Voom",sig,"_Boxplots",sep=""))
    if(unit=="[%]") perc <- 100 else perc <- 1
    setwd(paste(gsub("-","",sub("20","",Sys.Date())),"_",name,"Voom",sig,"_Boxplots",sep=""))
    for(i in which(voomRes$adj.P.Val<sig)){
      plotName <- rownames(voomRes)[i]
      annaMustPlotX(plotData[rownames(plotData)==plotName,]*perc,colData[,as.character(designFor)[2]],
                    name=gsub("/","",paste("boxplot","_",plotName,"_",name,sep="")),
                    ylab=paste(plotName, unit),xlab=as.character(designFor)[2],
                    pval=paste("adj p-val =",round(voomRes$adj.P.Val[i],3)),pch=pch,col=col,
                    secy=secy[plotName==secy[,1],2],font.secy=font.secy,silent=silent)
      Sys.sleep(0.1)
    }
    setwd("../")
  }
  if(res[1] == "count"){
    return(length(voomRes$adj.P.Val[voomRes$adj.P.Val<sig]))
  }else if(res[1] == "all"){
    return(voomRes)
  }else if(res[1] == "sig"){
    if(nrow(voomRes[voomRes$adj.P.Val<sig,])==0) warning("No significant findings.")
    return(voomRes[voomRes$adj.P.Val<sig,])
  }
}



#rarefaction curve
annaMustRarecurve <- function (x, step = 1, sample, xlab = "Sample Size", ylab = "Species", 
                               label = TRUE, col=Mcol[visFac], ...) {
  tot <- rowSums(x)
  S <- specnumber(x)
  nr <- nrow(x)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  par(mar=c(3.3,4.5,0.3,0.2),mgp=c(3.5,0.6,0))
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "", ylab = ylab, 
       type = "n", las=1, ...)
  mtext(xlab,1,2)
  if (!missing(sample)) {
    abline(v = sample)
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), 
                                           y = z, xout = sample, rule = 1)$y)
    abline(h = rare, lwd = 0.5)
  }
  for (ln in seq_len(length(out))) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col=col[ln],...)
  }
  if (label) {
    ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}
#distance measures
chaoDist <- function(mat,estimR){
  newdistChao <- vector()
  for(i in 1:(ncol(mat)-1)){
    for(j in (i+1):ncol(mat)){
      newdistChao <- append(newdistChao,((estimR[2,i]+estimR[2,j]-
                                            2*(length(which(mat[,i]>0&mat[,j]>0)) + estimR[3,i] + estimR[3,j]))
                                         /(estimR[2,i]+estimR[2,j])))
    }
  }
  return(newdistChao)
}
aceDist <- function(mat,estimR){
  newdistACE <- vector()
  for(i in 1:(ncol(mat)-1)){
    for(j in (i+1):ncol(mat)){
      newdistACE <- append(newdistACE,(estimR[4,i]+estimR[4,j]-
                                         2*(length(which(mat[,i]>0&mat[,j]>0)) + estimR[5,i] + estimR[5,j]))
                           /(estimR[4,i]+estimR[4,j]))
    }
  }
  return(newdistACE)
}

#confusion matrix:
annaMustPlotCM <- function(confMat){
  bwsca <- colorRampPalette(brewer.pal(9,"Greys"))
  dcol <- matrix(bwsca(256)[as.numeric(cut(c(confMat,0,sum(confMat)),breaks = 256))[1:length(c(confMat))]],nrow=dim(confMat)[1],ncol=dim(confMat)[2],dimnames=dimnames(confMat))
  xbound <- c(0:nrow(confMat))
  ybound <- c(ncol(confMat):0)
  par(mar=c(0.2,4.5,3,.2),mgp=c(5,0.1,0))
  plot(1,1,type="n",xlim=c(0,nrow(confMat)),ylim=c(0,nrow(confMat)),ann=F,axes=F,
       xaxs="i",yaxs="i")
  abline(h=ybound,lwd=0.5)
  abline(v=xbound,lwd=0.5)
  box(bty="o",lwd=0.5)
  for(i in 2:length(xbound)){
    for(j in 2:length(ybound)){
      polygon(x=c(xbound[i-1],rep(xbound[i],2),xbound[i-1]),y=rep(c(ybound[j-1],ybound[j]),each=2),
              col=dcol[j-1,i-1],lwd=0.1)
    }
  }
  text(rep((xbound[2:length(xbound)]+xbound[1:(length(xbound)-1)])/2,each=ncol(confMat)),
       rep((ybound[2:length(ybound)]+ybound[1:(length(ybound)-1)])/2,nrow(confMat)),label=c(confMat),cex=1)
  axis(3,at=(xbound[2:length(xbound)]+xbound[1:(length(xbound)-1)])/2,labels=colnames(confMat),las=1,tick=F,cex.axis=1)
  axis(2,at=(ybound[2:length(ybound)]+ybound[1:(length(ybound)-1)])/2,labels=rownames(confMat),las=2,tick=F,cex.axis=1)
  mtext("actual class",2,3.5,cex=1.2)
  mtext("predicted class",3,2,cex=1.2)
}

superRFC <- function(mat,classFac,title,ntree=20000,mtry=floor(sqrt(nrow(mat))),Mcol,Mchar,visFac,indiv,seed=NULL,silent=TRUE,ax=1:2,
                     axPerc=F){
  if(!is.null(seed)) set.seed(seed)
  classRF <- randomForest(t(mat),as.factor(classFac),importance=T,proximity=T,ntree=ntree,mtry=mtry)
  pdf(paste(gsub("-","",sub("20","",Sys.Date())),title,"_ConfusionMatrix.pdf",sep=""),width=3.5,height=3.0,pointsize=8)
  annaMustPlotCM(classRF$confusion[,-ncol(classRF$confusion)])
  dev.off()
  classMDS <- pcoa(1-classRF$proximity)
  classMDSval <- 100*c(classMDS$values[(1:3),2])
  if(axPerc){
    xlab <- paste("PC ",ax[1]," (",sprintf("%.1f",classMDSval[ax[1]]),"%)",sep="") 
    ylab <- paste("PC ",ax[2]," (",sprintf("%.1f",classMDSval[ax[2]]),"%)",sep="")
  }else{
    xlab <- paste("PC",ax[1])
    ylab <- paste("PC",ax[2])
  } 
  pdf(paste(gsub("-","",sub("20","",Sys.Date())),title,"_PCoAProximity.pdf",sep=""),width=4.2,height=3.5,pointsize=8)
  par(mar=c(1.5,1.9,0.2,0.2))
  layout(matrix(c(1,2),1,2,byrow=T),c(5/6,1/6))
  annaSclass(classMDS$vector[,ax], as.factor(visFac),col=Mcol[as.numeric(levels(as.factor((visFac))))], pch=Mchar[visFac],
             clabel=0, grid=0,mar=c(1.5,1.9,0.2,0.2),lwd=3,xlab=xlab,ylab=ylab)
  par(mar=c(0.2,0.2,0.2,0.2))
  plot(1,1,ann=F,type="n",axes=F)
  legend("topleft",legend=indiv,bty="n",pch=Mchar[as.numeric(levels(as.factor((visFac))))],
         col=Mcol[as.numeric(levels(as.factor((visFac))))],cex=1,xpd=T,y.intersp=1)
  dev.off()
  if(!silent){
    annaMustPlotCM(classRF$confusion[,-ncol(classRF$confusion)])
    Sys.sleep(0.1)
    par(mar=c(1.5,1.9,0.2,0.2))
    layout(matrix(c(1,2),1,2,byrow=T),c(5/6,1/6))
    annaSclass(classMDS$vector[,ax], as.factor(visFac),col=Mcol[as.numeric(levels(as.factor((visFac))))], pch=Mchar[visFac],
               clabel=0, grid=0,mar=c(1.5,1.9,0.2,0.2),lwd=3,xlab=xlab,ylab=ylab)
    par(mar=c(0.2,0.2,0.2,0.2))
    plot(1,1,ann=F,type="n",axes=F)
    legend("topleft",legend=indiv,bty="n",pch=Mchar[as.numeric(levels(as.factor((visFac))))],
           col=Mcol[as.numeric(levels(as.factor((visFac))))],cex=1,xpd=T,y.intersp=1)
    par(mar=c(3,4,1,0.2),mfrow=c(1,1))
  }
  classRF
}

superRFU <- function(mat,title,ntree=20000,mtry=floor(sqrt(nrow(mat))),Mcol,Mchar,visFac,indiv,seed=NULL,silent=TRUE,ax=1:2,
                     axPerc=F){
  if(!is.null(seed)) set.seed(seed)
  classRF <- randomForest(t(mat),importance=T,proximity=T,ntree=ntree,mtry=mtry)
  classMDS <- pcoa(1-classRF$proximity)
  classMDSval <- 100*c(classMDS$values[(1:3),2])
  if(axPerc){
    xlab <- paste("PC ",ax[1]," (",sprintf("%.1f",classMDSval[ax[1]]),"%)",sep="") 
    ylab <- paste("PC ",ax[2]," (",sprintf("%.1f",classMDSval[ax[2]]),"%)",sep="")
  }else{
    xlab <- paste("PC",ax[1])
    ylab <- paste("PC",ax[2])
  } 
  pdf(paste(gsub("-","",sub("20","",Sys.Date())),title,"_PCoAUnsupProximity.pdf",sep=""),width=4.2,height=3.5,pointsize=8)
  par(mar=c(1.5,1.9,0.2,0.2))
  layout(matrix(c(1,2),1,2,byrow=T),c(5/6,1/6))
  annaSclass(classMDS$vector[,ax], as.factor(visFac),col=Mcol[as.numeric(levels(as.factor((visFac))))], pch=Mchar[visFac],
             clabel=0, grid=0,mar=c(1.5,1.9,0.2,0.2),lwd=3,xlab=xlab,ylab=ylab)
  par(mar=c(0.2,0.2,0.2,0.2))
  plot(1,1,ann=F,type="n",axes=F)
  legend("topleft",legend=indiv,bty="n",pch=Mchar[as.numeric(levels(as.factor((visFac))))],
         col=Mcol[as.numeric(levels(as.factor((visFac))))],cex=1,xpd=T,y.intersp=1)
  dev.off()
  if(!silent){
    par(mar=c(1.5,1.9,0.2,0.2))
    layout(matrix(c(1,2),1,2,byrow=T),c(5/6,1/6))
    annaSclass(classMDS$vector[,ax], as.factor(visFac),col=Mcol[as.numeric(levels(as.factor((visFac))))], pch=Mchar[visFac],
               clabel=0, grid=0,mar=c(1.5,1.9,0.2,0.2),lwd=3,xlab=xlab,ylab=ylab)
    par(mar=c(0.2,0.2,0.2,0.2))
    plot(1,1,ann=F,type="n",axes=F)
    legend("topleft",legend=indiv,bty="n",pch=Mchar[as.numeric(levels(as.factor((visFac))))],
           col=Mcol[as.numeric(levels(as.factor((visFac))))],cex=1,xpd=T,y.intersp=1)
    par(mar=c(3,4,1,0.2),mfrow=c(1,1))
  }
  classRF
}

superRFR <- function(mat,regressVec,title,ntree=20000,mtry=max(floor(nrow(mat)/3)),Mcol,Mchar,visFac,seed=NULL,
                     silent=T,axlab="",legend=Mleg){
  if(!is.null(seed)) set.seed(seed)
  classRF <- randomForest(t(mat),regressVec,importance=T,ntree=ntree,mtry=mtry)
  y <- predict(classRF)
  xlab <- paste("actual",axlab)
  ylab <- paste("predicted",axlab)
  annaMustPlotI(regressVec,y,name=title,xlab=xlab,ylab=ylab,col=Mcol[visFac],pch=Mchar[visFac],
                type="pdf",legend=legend,legpch=Mchar,legcol=Mcol,main="",fmain=1,cmain=0.8,silent=silent,
                secy="",font.secy=1,font.leg=1,log="",cex=1)
  classRF
}

annaMustHeinzNW <- function(file,fdr=0.05,graph=koGraph,outname="Heinz"){
  data <- read.delim(file,stringsAsFactors=F,row.names=1)
  subnet <- subNetwork(rownames(data)[!is.na(data$pvalue)],graph)
  subnet <- rmSelfLoops(subnet)
  pval <- data$pvalue[!is.na(data$pvalue)]
  names(pval) <- rownames(data)[!is.na(data$pvalue)]
  fb <- fitBumModel(pval, plot = FALSE)
  scores <- scoreNodes(subnet, fb, fdr = fdr)
  writeHeinzEdges(network = subnet, file = paste(gsub(".tsv","",file),outname,sub("0.","",fdr),"Edges",sep="."),use.score = FALSE)
  writeHeinzNodes(network = subnet, file = paste(gsub(".tsv","",file),outname,sub("0.","",fdr),"Nodes",sep="."),node.scores = scores)
  system(paste("scp ",paste(gsub(".tsv","",file),outname,sub("0.","",fdr),"*es","txt",sep=".")," aheintzbuschart@gaia-cluster:/work/users/aheintzbuschart/testHeinz/",sep=""))
  return(list(file,paste(gsub(".tsv","",file),outname,sub("0.","",fdr),"Output.txt",sep=".")))
}

annaMustHeinzPlotNW <- function(file,outname="module.pdf",graph=koGraph,w2d=F){
  data <- read.delim(file[[1]],stringsAsFactors=F,row.names=1)
  subnet <- subNetwork(rownames(data)[!is.na(data$pvalue)],graph)
  subnet <- rmSelfLoops(subnet)
  fn <- unlist(strsplit(file[[2]],split="/"))
  fn <- fn[length(fn)]
  system(paste("scp aheintzbuschart@gaia-cluster:/work/users/aheintzbuschart/testHeinz/",fn," .",sep=""))
  module <- readHeinzGraph(node.file = fn, network = subnet)
  logFC <- data$log2FoldChange[!is.na(data$pvalue)]
  names(logFC) <- rownames(data)[!is.na(data$pvalue)]
  if(w2d){
    pdf(outname,width=3.5,height=3.5,pointsize=8)
    plotModule(module, scores = scores, diff.expr = logFC)
    dev.off()
  }
  plotModule(module, scores = scores, diff.expr = logFC)
  return(nodes(module))
}

annaMustPlotNW <- function(file,outname="module.pdf",fdr=0.05,graph=koGraph,w2d=F){
  data <- read.delim(file,stringsAsFactors=F,row.names=1)
  subnet <- subNetwork(rownames(data)[!is.na(data$pvalue)],graph)
  subnet <- rmSelfLoops(subnet)
  pval <- data$pvalue[!is.na(data$pvalue)]
  names(pval) <- rownames(data)[!is.na(data$pvalue)]
  fb <- fitBumModel(pval, plot = FALSE)
  scores <- scoreNodes(subnet, fb, fdr = fdr)
  module <- runFastHeinz(subnet, scores)
  logFC <- data$log2FoldChange[!is.na(data$pvalue)]
  names(logFC) <- rownames(data)[!is.na(data$pvalue)]
  if(w2d){
    pdf(outname,width=3.5,height=3.5,pointsize=8)
    plotModule(module, scores = scores, diff.expr = logFC)
    dev.off()
  }
  plotModule(module, scores = scores, diff.expr = logFC)
  return(nodes(module))
}

annaMustPrepModule <- function(nodes,diffRes,graph=koGraph,labCol=hmcol){
  modNW <- igraph.from.graphNEL(subNetwork(nodes,graph),name=T,weight=T)
  modLab <- merge(nodes,diffRes,by.x=1,by.y=0,all.x=T)$log2FoldChange 
  modLabC1 <- labCol[cut(c(modLab,-max(abs(modLab)),max(abs(modLab))),length(labCol))]
  modLabC2 <- labCol[length(labCol):1][cut(c(modLab,-max(abs(modLab)),max(abs(modLab))),length(labCol))]
  layo <- layout.fruchterman.reingold(modNW)
  return(list(modNW,modLabC1,modLabC2,layo))
}

annaMustPlotIndiModule <- function(nodes,nw,labCol,layout,indiv,outname,ome="T",silent=T,compc=compCol,clusterInf=cI,edge.width=2,edge.color="grey"){
  if(!all(ome %in% c("G","T","P"))){
    stop("invalid omic level, should be G,T or P")
  }else{
    if("G" %in% ome){
      indivExpr <- paste("/Users/anna.buschart/Documents/LCSB/Results/Must/KOnetwork/perSampleNodes2Pop/",
                         indiv,".DNAclusterPerNode.tsv",sep="")
      if(!file.exists(indivExpr)){
        stop("invalid individual, should be MX.X-VX")
      }else{
        nwPropG <- read.delim(indivExpr,stringsAsFactors=F)
        nwPropGmod <- merge(nodes,nwPropG,by.x=1,by.y=0,all.x=T)
        rownames(nwPropGmod) <- nwPropGmod[,1]
        nwPropGmod <- nwPropGmod[,-1]
        nwPropGmod[is.na(nwPropGmod)] <- 0
        nwPropGmod$absent <- ifelse(rowSums(nwPropGmod)==0,1,0)
        modCol <- c(sapply(colnames(nwPropGmod)[-ncol(nwPropGmod)],
                           function(x) compc[as.numeric(clusterInf$uniqueEss[clusterInf$cluster==x&clusterInf$sample==indiv])+1]),"grey")
        modAlph <- vector("numeric",length(modCol))
        first <- T
        for(i in 1:length(modAlph)){
          if(first) {
            first <- F
            modAlph[i] <- 1
          } else {
            if(cC==modCol[i]){
              modAlph[i] <- modAlph[i-1]*0.7
            }else{
              modAlph[i] <- 1
            }
          }
          cC <- modCol[i]
        }
        modCol <- list(alpha(modCol,modAlph))
        
        pdf(paste(outname,indiv,"metaG.pdf",sep="_"),width=7,height=7,pointsize=8)
        plot(nw,vertex.shape="pie",vertex.pie=as.list(as.data.frame(t(nwPropGmod))),
             vertex.size=4*log10(0.01+rowSums(nwPropGmod[,-1])), vertex.label.dist=.6,vertex.label.cex=1,
             vertex.label.family="sans",vertex.pie.color=modCol,vertex.label.color=labCol,vertex.label.font=2,layout=layout,
             edge.width=edge.width,edge.color=edge.color)
        dev.off()
        if(!silent) plot(nw,vertex.shape="pie",vertex.pie=as.list(as.data.frame(t(nwPropGmod))),
                         vertex.size=4*log10(0.01+rowSums(nwPropGmod[,-1])), vertex.label.dist=.5,vertex.label.cex=.6,
                         vertex.label.family="sans",vertex.pie.color=modCol,vertex.label.color=labCol,vertex.label.font=2,
                         layout=layout,edge.width=edge.width,edge.color=edge.color)
      }
    }
    if("T" %in% ome){
      indivExpr <- paste("/Users/anna.buschart/Documents/LCSB/Results/Must/KOnetwork/perSampleNodes2Pop/",
                         indiv,".RNAclusterPerNode.tsv",sep="")
      if(!file.exists(indivExpr)){
        stop("invalid individual, should be MX.X-VX")
      }else{
        nwPropT <- read.delim(indivExpr,stringsAsFactors=F)
        nwPropTmod <- merge(nodes,nwPropT,by.x=1,by.y=0,all.x=T)
        rownames(nwPropTmod) <- nwPropTmod[,1]
        nwPropTmod <- nwPropTmod[,-1]
        nwPropTmod[is.na(nwPropTmod)] <- 0
        nwPropTmod$absent <- ifelse(rowSums(nwPropTmod)==0,1,0)
        modCol <- c(sapply(colnames(nwPropTmod)[-ncol(nwPropTmod)],
                           function(x) compc[as.numeric(clusterInf$uniqueEss[clusterInf$cluster==x&clusterInf$sample==indiv])+1]),"grey")
        modAlph <- vector("numeric",length(modCol))
        first <- T
        for(i in 1:length(modAlph)){
          if(first) {
            first <- F
            modAlph[i] <- 1
          } else {
            if(cC==modCol[i]){
              modAlph[i] <- modAlph[i-1]*0.7
            }else{
              modAlph[i] <- 1
            }
          }
          cC <- modCol[i]
        }
        modCol <- list(alpha(modCol,modAlph))
        
        pdf(paste(outname,indiv,"metaT.pdf",sep="_"),width=7,height=7,pointsize=8)
        plot(nw,vertex.shape="pie",vertex.pie=as.list(as.data.frame(t(nwPropTmod))),
             vertex.size=4*log10(0.01+rowSums(nwPropTmod[,-1])), vertex.label.dist=.6,vertex.label.cex=1,
             vertex.label.family="sans",vertex.pie.color=modCol,vertex.label.color=labCol,vertex.label.font=2,layout=layout,
             edge.width=edge.width,edge.color=edge.color)
        dev.off()
        if(!silent) plot(nw,vertex.shape="pie",vertex.pie=as.list(as.data.frame(t(nwPropTmod))),
                         vertex.size=4*log10(0.01+rowSums(nwPropTmod[,-1])), vertex.label.dist=.5,vertex.label.cex=.6,
                         vertex.label.family="sans",vertex.pie.color=modCol,vertex.label.color=labCol,vertex.label.font=2,
                         layout=layout,edge.width=edge.width,edge.color=edge.color)
      }
    }
    if("P" %in% ome){
      indivExpr <- paste("/Users/anna.buschart/Documents/LCSB/Results/Must/KOnetwork/perSampleNodes2Pop/",
                         indiv,".PROTclusterPerNode.tsv",sep="")
      if(!file.exists(indivExpr)){
        stop("invalid individual, should be MX.X-VX")
      }else{
        nwPropP <- read.delim(indivExpr,stringsAsFactors=F)
        nwPropPmod <- merge(nodes,nwPropP,by.x=1,by.y=0,all.x=T)
        rownames(nwPropPmod) <- nwPropPmod[,1]
        nwPropPmod <- nwPropPmod[,-1]
        nwPropPmod[is.na(nwPropPmod)] <- 0
        nwPropPmod$absent <- ifelse(rowSums(nwPropPmod)==0,1,0)
        modCol <- c(sapply(colnames(nwPropPmod)[-ncol(nwPropPmod)],
                           function(x) compc[as.numeric(clusterInf$uniqueEss[clusterInf$cluster==x&clusterInf$sample==indiv])+1]),"grey")
        modAlph <- vector("numeric",length(modCol))
        first <- T
        for(i in 1:length(modAlph)){
          if(first) {
            first <- F
            modAlph[i] <- 1
          } else {
            if(cC==modCol[i]){
              modAlph[i] <- modAlph[i-1]*0.7
            }else{
              modAlph[i] <- 1
            }
          }
          cC <- modCol[i]
        }
        modCol <- list(alpha(modCol,modAlph))
        
        pdf(paste(outname,indiv,"metaP.pdf",sep="_"),width=7,height=7,pointsize=8)
        plot(nw,vertex.shape="pie",vertex.pie=as.list(as.data.frame(t(nwPropPmod))),
             vertex.size=4*log10(0.01+rowSums(nwPropPmod[,-1])), vertex.label.dist=.5,vertex.label.cex=.6,
             vertex.label.family="sans",vertex.pie.color=modCol,vertex.label.color=labCol,vertex.label.font=2,layout=layout,
             edge.width=edge.width,edge.color=edge.color)
        dev.off()
        if(!silent) plot(nw,vertex.shape="pie",vertex.pie=as.list(as.data.frame(t(nwPropPmod))),
                         vertex.size=4*log10(0.01+rowSums(nwPropPmod[,-1])), vertex.label.dist=.5,vertex.label.cex=.6,
                         vertex.label.family="sans",vertex.pie.color=modCol,vertex.label.color=labCol,vertex.label.font=2,
                         layout=layout,edge.width=edge.width,edge.color=edge.color)
      }
    }
  }
}

plotModSuper <- function(nodes,diffRes,indiv,outname,groups=1,ome="T",silent=T,compc=compCol,clusterInf=cI,graph=koGraph,labCol=hmcol){
  if(any(groups>2)){
    stop("only 2 groups possible")
  }else{
    val <- annaMustPrepModule(nodes,diffRes,graph=koGraph,labCol=hmcol)
    if(length(groups)>1){
      for(i in unique(groups)){
        indi <- indiv[groups==i]
        for(ind in indi){
          annaMustPlotIndiModule(nodes,val[[1]],val[[1+i]],val[[4]],ind,outname,ome)
        }
      }
    }else{
      for(ind in indiv){
        annaMustPlotIndiModule(nodes,val[[1]],val[[1+groups]],val[[4]],ind,outname,ome)
      }
    }
  }
  return(list(val[[1]],val[[4]]))
}

diffPathHM <- function(diffTab,outName,sigFC=1,sigDiff=0.05,bgset=rownames(koAbundR),sigEnrich=0.05,visFac=visFac,
                       plotMat=koAbundN,colsep=NULL,ordervec=c(1:ncol(plotMat)),returnSig=T){
  diffTab <- data.frame(diffTab,stringsAsFactors=F)
  diffTab$KO <- rownames(diffTab)
  diffTabSig <- diffTab[abs(diffTab$log2FoldChange)>=sigFC & !is.na(diffTab$padj) & diffTab$padj<=sigDiff,]
  
  heatmap.2a(as.matrix(plotMat[rownames(plotMat) %in% diffTabSig$KO,ordervec])[order(diffTabSig$log2FoldChange[
    order(diffTabSig$KO)]),],trace="none",Colv="none",dendrogram="none",Rowv="none",
    col=hmcol,scale="row",sepcolor="black",colsep=colsep,
    ColSideColor=rbind(Mcol[visFac],c("white","black")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="Yes")])[,ordervec])
  pdf(paste(gsub("-","",sub("20","",Sys.Date())),outName,"_heatmapFC.pdf",sep=""),width=5,height=5,pointsize=8)
  heatmap.2a(as.matrix(plotMat[rownames(plotMat) %in% diffTabSig$KO,ordervec])[order(diffTabSig$log2FoldChange[
    order(diffTabSig$KO)]),],trace="none",Colv="none",dendrogram="none",Rowv="none",
    col=hmcol,scale="row",sepcolor="black",colsep=colsep,
    ColSideColor=rbind(Mcol[visFac],c("white","black")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="Yes")])[,ordervec])
  dev.off()
  
  enrichSig <- annaEnrichPath(diffTabSig,bgset,"KO",sigEnrich)
  pathCol <- diffTab$log2FoldChange
  names(pathCol) <- diffTab$KO
  print(enrichSig$pathway)
  for(i in enrichSig$pathway[!enrichSig$pathway %in% c("ko00121","ko01120","ko01054","ko01200")]){ 
    pv <- pathview(pathCol,pathway.id=i,species="ko",out.suffix=outName,low = list(gene = "#2C7BB6", cpd = "green"), 
                   mid =list(gene = "#FFFFBF", cpd = "yellow"), high = list(gene = "#D7191C", cpd ="pink"))
  }
  if(returnSig) return(diffTabSig)
} 

diffPathHMNW <- function(diffTab,outName,sigFC=1,sigDiff=0.05,sigEnrich=0.05,visFac=visFac,
                         plotMat=koAbundNNW,colsep=NULL,ordervec=c(1:ncol(plotMat)),returnSig=T){
  diffTab <- data.frame(diffTab,stringsAsFactors=F)
  diffTab$KO <- rownames(diffTab)
  diffTabSig <- diffTab[abs(diffTab$log2FoldChange)>=sigFC & !is.na(diffTab$padj) & diffTab$padj<=sigDiff,]
  
  heatmap.2a(as.matrix(plotMat[rownames(plotMat) %in% diffTabSig$KO,ordervec])[order(diffTabSig$log2FoldChange[
    order(diffTabSig$KO)]),],trace="none",Colv="none",dendrogram="none",Rowv="none",
    col=hmcol,scale="row",sepcolor="black",colsep=colsep,margins=c(5,10),
    ColSideColor=rbind(Mcol[visFac],c("white","black")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="Yes")])[,ordervec])
  pdf(paste(gsub("-","",sub("20","",Sys.Date())),outName,"_heatmapFC.pdf",sep=""),width=5,height=5,pointsize=8)
  heatmap.2a(as.matrix(plotMat[rownames(plotMat) %in% diffTabSig$KO,ordervec])[order(diffTabSig$log2FoldChange[
    order(diffTabSig$KO)]),],trace="none",Colv="none",dendrogram="none",Rowv="none",
    col=hmcol,scale="row",sepcolor="black",colsep=colsep,margins=c(5,10),
    ColSideColor=rbind(Mcol[visFac],c("white","black")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="Yes")])[,ordervec])
  dev.off()
  
  sigKO <- unlist(sapply(diffTabSig$KO,function(x) unlist(strsplit(x,split="-"))))
  sigPW <- merge(sigKO,pw,by.x=1,by.y="ko")
  sigPW <- unique(sigPW)
  colnames(sigPW)[1] <- "koID"
  bgSet <- merge(unlist(sapply(nwNodes$V1,function(x) unlist(strsplit(x,split="-")))),pw,by.x=1,by.y="ko")
  colnames(bgSet)[1] <- "ko" 
  bgSet <- unique(bgSet)
  pwInf <- data.frame("pathway"=unique(sigPW$pw),"numberInSig","KOs","numberTotal","pval",stringsAsFactors=F)
  colnames(pwInf) <- gsub(".","",gsub("X.","",colnames(pwInf)),fixed=T)
  
  for(i in pwInf$pathway){
    q <- length(unique(sigPW$koID[sigPW$pw==i]))
    m <- length(unique(bgSet$ko[bgSet$pw==i]))
    n <- length(unique(bgSet$ko)) - m
    k <- length(unique(sigPW$koID))
    pwInf$numberInSig[pwInf$pathway==i] <- q
    pwInf$KOs[pwInf$pathway==i] <- paste(sigPW$koID[which(sigPW$pw==i)],sep=";",collapse=";")
    pwInf$numberTotal[pwInf$pathway==i] <- m
    pwInf$pval[pwInf$pathway==i] <- phyper(q-1,m,n,k,lower.tail=F) 
  }
  pwInf[,6] <- p.adjust(pwInf[,5],"fdr")
  colnames(pwInf)[6] <- "adjpval"
  
  pwEn <- pwInf[ pwInf$adjpval<=sigEnrich,]
  pwEn <- merge(pwEn,pwn,by.x="pathway",by.y="pw",all.x=T)
  
  allKO <- unlist(sapply(diffTab$KO,function(x) unlist(strsplit(x,split="-"))))
  pathCol <- sapply(allKO,function(x) diffTab$log2FoldChange[grep(x,diffTab$KO)])
  names(pathCol) <- allKO
  print(pwEn$pathway)
  for(i in pwEn$pathway[!pwEn$pathway %in% c("ko00121","ko01120","ko01054","ko01200")]){ 
    pv <- pathview(pathCol,pathway.id=i,species="ko",out.suffix=outName,low = list(gene = "#2C7BB6", cpd = "green"), 
                   mid =list(gene = "#FFFFBF", cpd = "yellow"), high = list(gene = "#D7191C", cpd ="pink"))
  }
  if(returnSig) return(diffTabSig)
}

nodes2Names <- function(nodes){
  kon[kon$koID %in% unlist(sapply(nodes,function(x) unlist(strsplit(x,split="-")))),]
}

betaSor.multi.adj <- function (x, missing, index.family = "sorensen") {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.numeric(x)) 
    stop("The data in x is not numeric.", call. = TRUE)
  xvals <- unique(as.vector(x))
  if (any(!is.element(xvals, c(0, 1)))) 
    stop("The table contains values other than 0 and 1: data should be presence/absence.",call. = TRUE)
  x <- x[,order(colSums(x),decreasing=T)]
  addCol <- (rowSums(x)+missing)-ncol(x)
  addCol[addCol<0] <- 0
  while(max(addCol)>0){
    x <- cbind(x,rep(0,nrow(x)))
    addCol <- addCol -1 
  }
  for(i in 1:nrow(x)){
    if(missing[i]>0){
      x[i,x[i,]==0][1:missing[i]] <- 1
    }
  }
  shared <- x %*% t(x)
  not.shared <- abs(sweep(shared, 2, diag(shared)))
  sumSi <- sum(diag(shared))
  St <- sum(colSums(x) > 0)
  a <- sumSi - St
  max.not.shared <- pmax(not.shared, t(not.shared))
  min.not.shared <- pmin(not.shared, t(not.shared))
  x <- list(a = a, max.not.shared = max.not.shared, min.not.shared = min.not.shared)
  maxbibj <- sum(x$max.not.shared[lower.tri(x$max.not.shared)])
  minbibj <- sum(x$min.not.shared[lower.tri(x$min.not.shared)])
  beta.sor <- (minbibj + maxbibj)/(minbibj + maxbibj + (2 * x$a))
  return(beta.sor)
}

anna5tab2venn <- function(data,noCount="unknown",col=brewer.pal(5,"Set1")){
  plot(1,type="n",ann=F,axes=F)
  draw.quintuple.venn(area1 =length(which(data[,1] != noCount)),
                      area2 =length(which(data[,2] != noCount)),
                      area3 =length(which(data[,3] != noCount)),
                      area4 =length(which(data[,4] != noCount)),
                      area5 =length(which(data[,5] != noCount)),
                      n12 =length(which(data[,1] != noCount & data[,2] != noCount)), 
                      n13=length(which(data[,1] != noCount & data[,3] != noCount)), 
                      n14=length(which(data[,1] != noCount & data[,4] != noCount)), 
                      n15=length(which(data[,1] != noCount & data[,5] != noCount)),
                      n23=length(which(data[,2] != noCount & data[,3] != noCount)), 
                      n24=length(which(data[,2] != noCount & data[,4] != noCount)), 
                      n25=length(which(data[,2] != noCount & data[,5] != noCount)), 
                      n34=length(which(data[,3] != noCount & data[,4] != noCount)), 
                      n35=length(which(data[,3] != noCount & data[,5] != noCount)), 
                      n45=length(which(data[,4] != noCount & data[,5] != noCount)), 
                      n123=length(which(data[,1] != noCount & data[,2] != noCount& data[,3] != noCount)), 
                      n124=length(which(data[,1] != noCount & data[,2] != noCount& data[,4] != noCount)), 
                      n125=length(which(data[,1] != noCount & data[,2] != noCount& data[,5] != noCount)), 
                      n134=length(which(data[,1] != noCount & data[,3] != noCount& data[,4] != noCount)), 
                      n135=length(which(data[,1] != noCount & data[,3] != noCount& data[,5] != noCount)),
                      n145=length(which(data[,1] != noCount & data[,4] != noCount& data[,5] != noCount)), 
                      n234=length(which(data[,2] != noCount & data[,3] != noCount& data[,4] != noCount)), 
                      n235=length(which(data[,2] != noCount & data[,3] != noCount& data[,5] != noCount)), 
                      n245=length(which(data[,2] != noCount & data[,4] != noCount& data[,5] != noCount)), 
                      n345=length(which(data[,3] != noCount & data[,4] != noCount& data[,5] != noCount)), 
                      n1234=length(which(data[,1] != noCount & data[,2] != noCount& data[,3] != noCount& data[,4] != noCount )), 
                      n1235=length(which(data[,1] != noCount & data[,2] != noCount& data[,3] != noCount& data[,5] != noCount )), 
                      n1245=length(which(data[,1] != noCount & data[,2] != noCount& data[,4] != noCount& data[,5] != noCount )), 
                      n1345=length(which(data[,1] != noCount & data[,3] != noCount& data[,4] != noCount& data[,5] != noCount )), 
                      n2345=length(which(data[,2] != noCount & data[,3] != noCount& data[,4] != noCount& data[,5] != noCount )), 
                      n12345=length(which(data[,1] != noCount & data[,2] != noCount& data[,3] != noCount& data[,4] != noCount & data[,5] != noCount )),
                      category=colnames(data),col=col,fil=col,alpha=rep(0.3,5),cex=0.8,fontfamily=rep("sans",31),
                      cat.fontfamily=rep("sans",5),cat.just=list(c(0.5,4),c(0,.5),c(0,-3),c(1,-2),c(1,.5))
  )
}

anna4tab2venn <- function(data,noCount="unknown",col=brewer.pal(4,"Set1")){
  plot(1,type="n",ann=F,axes=F)
  draw.quad.venn(area1 =length(which(data[,1] != noCount)),
                 area2 =length(which(data[,2] != noCount)),
                 area3 =length(which(data[,3] != noCount)),
                 area4 =length(which(data[,4] != noCount)),
                 n12 =length(which(data[,1] != noCount & data[,2] != noCount)), 
                 n13=length(which(data[,1] != noCount & data[,3] != noCount)), 
                 n14=length(which(data[,1] != noCount & data[,4] != noCount)), 
                 n23=length(which(data[,2] != noCount & data[,3] != noCount)), 
                 n24=length(which(data[,2] != noCount & data[,4] != noCount)), 
                 n34=length(which(data[,3] != noCount & data[,4] != noCount)), 
                 n123=length(which(data[,1] != noCount & data[,2] != noCount& data[,3] != noCount)), 
                 n124=length(which(data[,1] != noCount & data[,2] != noCount& data[,4] != noCount)), 
                 n134=length(which(data[,1] != noCount & data[,3] != noCount& data[,4] != noCount)), 
                 n234=length(which(data[,2] != noCount & data[,3] != noCount& data[,4] != noCount)), 
                 n1234=length(which(data[,1] != noCount & data[,2] != noCount& data[,3] != noCount& data[,4] != noCount )), 
                 category=colnames(data),col=col,fil=col[c(1,4,2,3)],alpha=rep(0.3,4),cex=0.8,fontfamily=rep("sans",15),
                 cat.fontfamily=rep("sans",4))
}

anna3tab2venn <- function(data,noCount="unknown",col=brewer.pal(3,"Set1")){
  plot(1,type="n",ann=F,axes=F)
  draw.triple.venn(area1 =length(which(data[,1] != noCount)),
                   area2 =length(which(data[,2] != noCount)),
                   area3 =length(which(data[,3] != noCount)),
                   n12 =length(which(data[,1] != noCount & data[,2] != noCount)), 
                   n13=length(which(data[,1] != noCount & data[,3] != noCount)), 
                   n23=length(which(data[,2] != noCount & data[,3] != noCount)), 
                   n123=length(which(data[,1] != noCount & data[,2] != noCount& data[,3] != noCount)), 
                   category=colnames(data),col=col,fil=col,alpha=rep(0.3,3),cex=0.8,fontfamily=rep("sans",7),
                   cat.fontfamily=rep("sans",3)
  )
}


diffHM <- function(diffTab,outName,sigFC=1,sigDiff=0.05,visFac=visFac,plotMat,colsep=NULL,ordervec=c(1:ncol(plotMat))){
  diffTab <- data.frame(diffTab,stringsAsFactors=F)
  diffTab$ID <- rownames(diffTab)
  diffTabSig <- diffTab[abs(diffTab$log2FoldChange)>=sigFC & !is.na(diffTab$padj) & diffTab$padj<=sigDiff,]
  
  heatmap.2a(as.matrix(plotMat[rownames(plotMat) %in% diffTabSig$ID,ordervec])[order(diffTabSig$log2FoldChange[
    order(diffTabSig$ID)]),],trace="none",Colv="none",dendrogram="none",Rowv="none",margin=c(4,8),
    col=hmcol,scale="row",sepcolor="black",colsep=colsep,
    ColSideColor=rbind(Mcol[visFac],c("white","black")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="Yes")])[,ordervec])
  pdf(paste(gsub("-","",sub("20","",Sys.Date())),outName,"_heatmapFC.pdf",sep=""),width=5,height=5,pointsize=8)
  heatmap.2a(as.matrix(plotMat[rownames(plotMat) %in% diffTabSig$ID,ordervec])[order(diffTabSig$log2FoldChange[
    order(diffTabSig$ID)]),],trace="none",Colv="none",dendrogram="none",Rowv="none",
    col=hmcol,scale="row",sepcolor="black",colsep=colsep,margin=c(4,8),
    ColSideColor=rbind(Mcol[visFac],c("white","black")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="Yes")])[,ordervec])
  dev.off()
}


# kruskal test with and without multiple testing adjustment + plots
kruskalSuper <- function(mat,x,grpName,sigMT=0.05,sig=0.01,outName,pch=Mchar[visFac],col=Mcol[visFac]
                         ,secy=matrix("",ncol=2,nrow=1),unit="[%]",silent=TRUE,font.secy=1){
  annaMustKrusPlotMT(mat=mat,x=x,grpName=grpName,sigMT=sigMT,outName=outName,pch=pch,col=col,secy=secy,unit=unit,silent=silent,font.secy=font.secy)
  annaMustKrusPlotWOMT(mat=mat,x=x,grpName=grpName,sig=sig,outName=outName,pch=pch,col=col,secy=secy,unit=unit,silent=silent,font.secy=font.secy)
}

# kruskal test and boxplot of significant observations (with multiple testing adjustment)
annaMustKrusPlotMT <- function(mat,x,grpName,sigMT,outName,pch=Mchar[visFac],col=Mcol[visFac],unit="[%]",
                               secy=matrix("",ncol=2,nrow=1),font.secy=1,silent=TRUE){
  krusVec <- vector(length=nrow(mat))
  names(krusVec) <- rownames(mat)
  for(i in 1:nrow(mat)) krusVec[i] <- kruskal.test(as.matrix(mat)[i,],x)$p.value
  krusVec[is.na(krusVec)] <- 1
  write.table(data.frame("name"=rownames(mat),"pval"=krusVec,"adj.pval"=p.adjust(krusVec,"fdr"),stringsAsFactors=F),
              paste(gsub("-","",sub("20","",Sys.Date())),"-",outName,"MTlt",sigMT,".tsv",sep=""),sep="\t",row.names=F,quote=F)
  if(length(krusVec[p.adjust(krusVec,"fdr")<=sigMT])>0){
    system(paste("mkdir ",gsub("-","",sub("20","",Sys.Date())),"_",outName,"MTlt",sigMT,"_Boxplots",sep=""))
    if(unit=="[%]") perc <- 100 else perc <- 1
    setwd(paste(gsub("-","",sub("20","",Sys.Date())),"_",outName,"MTlt",sigMT,"_Boxplots",sep=""))
    for(i in which(p.adjust(krusVec,"fdr")<=sigMT)) {
      pdf(paste(gsub("-","",sub("20","",Sys.Date())),"-",gsub("/","",paste("boxplot",rownames(mat)[i],grpName,sep="")),".pdf",sep=""),width=3.5,height=2.5,pointsize=8)
      par(mar=c(4.3,4.5,0.3,0.2),mgp=c(3.5,0.6,0))
      miny <- 0
      maxy <- 0
      y <- as.matrix(mat)[i,]*perc
      if(min(y)>=0) miny <- min(y)*0.9 else miny <- min(y)*1.1
      if(max(y)>=0) maxy <- max(y)*1.1 else maxy <- max(y)*0.9
      boxplot(y ~ x, ylab=paste(rownames(mat)[i], unit),bty="l",las=1,outline=F,ylim=c(miny,maxy))
      points(jitter(as.numeric(as.factor(x)),0.5),y,pch=pch,col=col,cex=0.8)
      mtext(paste("adj p-val =",round(p.adjust(krusVec[i],"fdr"),3)),1,3)
      mtext(grpName,1,2)
      if(length(secy)>2) mtext(strtrim(secy[rownames(mat)[i]==secy[,1],2],40),2,2.5,font=font.secy)    
      dev.off()
    }
    setwd("../")
  }else warning("No significant observations!")
}
# Kruskal test and boxplot of significant observations (without multiple testing adjustment)
annaMustKrusPlotWOMT <- function(mat,x,grpName,sig,outName,pch=Mchar[visFac],col=Mcol[visFac],unit="[%]",
                                 secy=matrix("",ncol=2,nrow=1),font.secy=1,silent=TRUE){
  krusVec <- vector(length=nrow(mat))
  names(krusVec) <- rownames(mat)
  for(i in 1:nrow(mat)) krusVec[i] <- kruskal.test(as.matrix(mat)[i,],x)$p.value
  krusVec[is.na(krusVec)] <- 1
  if(length(krusVec[krusVec<=sig])>0){
    system(paste("mkdir ",gsub("-","",sub("20","",Sys.Date())),"_",outName,"pvallt",sig,"_Boxplots",sep=""))
    if(unit=="[%]") perc <- 100 else perc <- 1
    setwd(paste(gsub("-","",sub("20","",Sys.Date())),"_",outName,"pvallt",sig,"_Boxplots",sep=""))
    for(i in which(krusVec<=sig)) {
      pdf(paste(gsub("-","",sub("20","",Sys.Date())),"-",gsub("/","",paste("boxplot",rownames(mat)[i],grpName,sep="")),".pdf",sep=""),width=3.5,height=2.5,pointsize=8)
      par(mar=c(4.3,4.5,0.3,0.2),mgp=c(3.5,0.6,0))
      miny <- 0
      maxy <- 0
      y <- as.matrix(mat)[i,]*perc
      if(min(y)>=0) miny <- min(y)*0.9 else miny <- min(y)*1.1
      if(max(y)>=0) maxy <- max(y)*1.1 else maxy <- max(y)*0.9
      boxplot(y ~ x, ylab=paste(rownames(mat)[i], unit),bty="l",las=1,outline=F,ylim=c(miny,maxy))
      points(jitter(as.numeric(as.factor(x)),0.5),y,pch=pch,col=col,cex=0.8)
      mtext(paste("p-val =",round(krusVec[i],3)),1,3)
      mtext(grpName,1,2)
      if(length(secy)>2) mtext(strtrim(secy[rownames(mat)[i]==secy[,1],2],40),2,2.5,font=font.secy)
      dev.off()
    }
    setwd("../")
  }else warning("No significant observations!")
}


bxpR <- function (z, notch = FALSE, width = NULL, varwidth = FALSE, outline = TRUE, 
                  notch.frac = 0.5, log = "", border = par("fg"), pars = NULL, 
                  frame.plot = axes, horizontal = FALSE, add = FALSE, at = NULL, 
                  show.names = NULL, bgcol=par("bg"), ...) 
{
  pars <- c(list(...), pars)
  pars <- pars[unique(names(pars))]
  bplt <- function(x, wid, stats, out, conf, notch, xlog, i) {
    ok <- TRUE
    if (!any(is.na(stats))) {
      xP <- if (xlog) 
        function(x, w) x * exp(w)
      else function(x, w) x + w
      wid <- wid/2
      if (notch) {
        ok <- stats[2L] <= conf[1L] && conf[2L] <= stats[4L]
        xx <- xP(x, wid * c(-1, 1, 1, notch.frac, 1, 
                            1, -1, -1, -notch.frac, -1))
        yy <- c(stats[c(2, 2)], conf[1L], stats[3L], 
                conf[2L], stats[c(4, 4)], conf[2L], stats[3L], 
                conf[1L])
      }
      else {
        xx <- xP(x, wid * c(0, 1, 1, 0))
        yy <- stats[c(2, 2, 4, 4)]
      }
      if (!notch) 
        notch.frac <- 1
      wntch <- notch.frac * wid
      xypolygon(xx, yy, lty = "blank", col = boxfill[i])
      xysegments(xP(x, 0), stats[3L], xP(x, +wntch), 
                 stats[3L], lty = medlty[i], lwd = medlwd[i], 
                 col = medcol[i], lend = 1)
      xypoints(x, stats[3L], pch = medpch[i], cex = medcex[i], 
               col = medcol[i], bg = medbg[i])
      xysegments(rep.int(x, 2), stats[c(1, 5)], rep.int(x, 
                                                        2), stats[c(2, 4)], lty = whisklty[i], lwd = whisklwd[i], 
                 col = whiskcol[i])
      xysegments(rep.int(xP(x, 0), 2), 
                 stats[c(1, 5)], rep.int(xP(x, +wid * staplewex[i]), 
                                         2), stats[c(1, 5)], lty = staplelty[i], lwd = staplelwd[i], 
                 col = staplecol[i])
      xypolygon(xx, yy, lty = boxlty[i], lwd = boxlwd[i], 
                border = boxcol[i])
      if ((nout <- length(out))) {
        xysegments(rep(x - wid * outwex, nout), out, 
                   rep(x + wid * outwex, nout), out, lty = outlty[i], 
                   lwd = outlwd[i], col = outcol[i])
        xypoints(rep.int(x, nout), out, pch = outpch[i], 
                 lwd = outlwd[i], cex = outcex[i], col = outcol[i], 
                 bg = outbg[i])
      }
      if (any(inf <- !is.finite(out))) {
        warning(sprintf(ngettext(length(unique(out[inf])), 
                                 "Outlier (%s) in boxplot %d is not drawn", 
                                 "Outliers (%s) in boxplot %d are not drawn"), 
                        paste(unique(out[inf]), collapse = ", "), i), 
                domain = NA)
      }
    }
    return(ok)
  }
  if (!is.list(z) || 0L == (n <- length(z$n))) 
    stop("invalid first argument")
  if (is.null(at)) 
    at <- 1L:n
  else if (length(at) != n) 
    stop("'at' must have same length as 'z$n', i.e. ", n)
  if (is.null(z$out)) 
    z$out <- numeric()
  if (is.null(z$group) || !outline) 
    z$group <- integer()
  if (is.null(pars$ylim)) 
    ylim <- range(z$stats[is.finite(z$stats)], if (outline) z$out[is.finite(z$out)], 
                  if (notch) z$conf[is.finite(z$conf)])
  else {
    ylim <- pars$ylim
    pars$ylim <- NULL
  }
  if (length(border) == 0L) 
    border <- par("fg")
  dev.hold()
  on.exit(dev.flush())
  if (!add) {
    if (is.null(pars$xlim)) 
      xlim <- range(at, finite = TRUE) + c(-0.5, 0.5)
    else {
      xlim <- pars$xlim
      pars$xlim <- NULL
    }
    plot.new()
    if (horizontal) 
      plot.window(ylim = xlim, xlim = ylim, log = log, 
                  xaxs = pars$yaxs)
    else plot.window(xlim = xlim, ylim = ylim, log = log, 
                     yaxs = pars$yaxs)
  }
  xlog <- (par("ylog") && horizontal) || (par("xlog") && !horizontal)
  pcycle <- function(p, def1, def2 = NULL) rep(if (length(p)) p else if (length(def1)) def1 else def2, 
                                               length.out = n)
  p <- function(sym) pars[[sym, exact = TRUE]]
  boxlty <- pcycle(pars$boxlty, p("lty"), par("lty"))
  boxlwd <- pcycle(pars$boxlwd, p("lwd"), par("lwd"))
  boxcol <- pcycle(pars$boxcol, border)
  boxfill <- pcycle(pars$boxfill, bgcol)
  boxwex <- pcycle(pars$boxwex, 0.8 * {
    if (n <= 1) 
      1
    else stats::quantile(diff(sort(if (xlog) 
      log(at)
      else at)), 0.1)
  })
  medlty <- pcycle(pars$medlty, p("lty"), par("lty"))
  medlwd <- pcycle(pars$medlwd, 3 * p("lwd"), 3 * par("lwd"))
  medpch <- pcycle(pars$medpch, NA_integer_)
  medcex <- pcycle(pars$medcex, p("cex"), par("cex"))
  medcol <- pcycle(pars$medcol, border)
  medbg <- pcycle(pars$medbg, p("bg"), par("bg"))
  whisklty <- pcycle(pars$whisklty, p("lty"), "dashed")
  whisklwd <- pcycle(pars$whisklwd, p("lwd"), par("lwd"))
  whiskcol <- pcycle(pars$whiskcol, border)
  staplelty <- pcycle(pars$staplelty, p("lty"), par("lty"))
  staplelwd <- pcycle(pars$staplelwd, p("lwd"), par("lwd"))
  staplecol <- pcycle(pars$staplecol, border)
  staplewex <- pcycle(pars$staplewex, 0.5)
  outlty <- pcycle(pars$outlty, "blank")
  outlwd <- pcycle(pars$outlwd, p("lwd"), par("lwd"))
  outpch <- pcycle(pars$outpch, p("pch"), par("pch"))
  outcex <- pcycle(pars$outcex, p("cex"), par("cex"))
  outcol <- pcycle(pars$outcol, border)
  outbg <- pcycle(pars$outbg, p("bg"), par("bg"))
  outwex <- pcycle(pars$outwex, 0.5)
  width <- if (!is.null(width)) {
    if (length(width) != n | any(is.na(width)) | any(width <= 
                                                     0)) 
      stop("invalid boxplot widths")
    boxwex * width/max(width)
  }
  else if (varwidth) 
    boxwex * sqrt(z$n/max(z$n))
  else if (n == 1) 
    0.5 * boxwex
  else rep.int(boxwex, n)
  if (horizontal) {
    xypoints <- function(x, y, ...) points(y, x, ...)
    xypolygon <- function(x, y, ...) polygon(y, x, ...)
    xysegments <- function(x0, y0, x1, y1, ...) segments(y0, 
                                                         x0, y1, x1, ...)
  }
  else {
    xypoints <- points
    xypolygon <- polygon
    xysegments <- segments
  }
  ok <- TRUE
  for (i in 1L:n) ok <- ok & bplt(at[i], wid = width[i], stats = z$stats[, 
                                                                         i], out = z$out[z$group == i], conf = z$conf[, i], notch = notch, 
                                  xlog = xlog, i = i)
  if (!ok) 
    warning("some notches went outside hinges ('box'): maybe set notch=FALSE")
  axes <- is.null(pars$axes)
  if (!axes) {
    axes <- pars$axes
    pars$axes <- NULL
  }
  if (axes) {
    ax.pars <- pars[names(pars) %in% c("xaxt", "yaxt", "xaxp", 
                                       "yaxp", "las", "cex.axis", "col.axis", "format")]
    if (is.null(show.names)) 
      show.names <- n > 1
    if (show.names) 
      do.call("axis", c(list(side = 1 + horizontal, at = at, 
                             labels = z$names), ax.pars))
    do.call("Axis", c(list(x = z$stats, side = 2 - horizontal), 
                      ax.pars))
  }
  do.call("title", pars[names(pars) %in% c("main", "cex.main", 
                                           "col.main", "sub", "cex.sub", "col.sub", "xlab", "ylab", 
                                           "cex.lab", "col.lab")])
  if (frame.plot) 
    box()
  invisible(at)
}

bxpL <- function (z, notch = FALSE, width = NULL, varwidth = FALSE, outline = TRUE, 
                  notch.frac = 0.5, log = "", border = par("fg"), pars = NULL, 
                  frame.plot = axes, horizontal = FALSE, add = FALSE, at = NULL, 
                  show.names = NULL, bgcol=par("bg"),...) 
{
  pars <- c(list(...), pars)
  pars <- pars[unique(names(pars))]
  bplt <- function(x, wid, stats, out, conf, notch, xlog, i) {
    ok <- TRUE
    if (!any(is.na(stats))) {
      xP <- if (xlog) 
        function(x, w) x * exp(w)
      else function(x, w) x + w
      wid <- wid/2
      if (notch) {
        ok <- stats[2L] <= conf[1L] && conf[2L] <= stats[4L]
        xx <- xP(x, wid * c(-1, 1, 1, notch.frac, 1, 
                            1, -1, -1, -notch.frac, -1))
        yy <- c(stats[c(2, 2)], conf[1L], stats[3L], 
                conf[2L], stats[c(4, 4)], conf[2L], stats[3L], 
                conf[1L])
      }
      else {
        xx <- xP(x, wid * c(-1, 0, 0, -1))
        yy <- stats[c(2, 2, 4, 4)]
      }
      if (!notch) 
        notch.frac <- 1
      wntch <- notch.frac * wid
      xypolygon(xx, yy, lty = "blank", col = boxfill[i])
      xysegments(xP(x, -wntch), stats[3L], xP(x, 0), 
                 stats[3L], lty = medlty[i], lwd = medlwd[i], 
                 col = medcol[i], lend = 1)
      xypoints(x, stats[3L], pch = medpch[i], cex = medcex[i], 
               col = medcol[i], bg = medbg[i])
      xysegments(rep.int(x, 2), stats[c(1, 5)], rep.int(x, 
                                                        2), stats[c(2, 4)], lty = whisklty[i], lwd = whisklwd[i], 
                 col = whiskcol[i])
      xysegments(rep.int(xP(x, -wid * staplewex[i]), 2), 
                 stats[c(1, 5)], rep.int(xP(x, 0), 
                                         2), stats[c(1, 5)], lty = staplelty[i], lwd = staplelwd[i], 
                 col = staplecol[i])
      xypolygon(xx, yy, lty = boxlty[i], lwd = boxlwd[i], 
                border = boxcol[i])
      if ((nout <- length(out))) {
        xysegments(rep(x - wid * outwex, nout), out, 
                   rep(x + wid * outwex, nout), out, lty = outlty[i], 
                   lwd = outlwd[i], col = outcol[i])
        xypoints(rep.int(x, nout), out, pch = outpch[i], 
                 lwd = outlwd[i], cex = outcex[i], col = outcol[i], 
                 bg = outbg[i])
      }
      if (any(inf <- !is.finite(out))) {
        warning(sprintf(ngettext(length(unique(out[inf])), 
                                 "Outlier (%s) in boxplot %d is not drawn", 
                                 "Outliers (%s) in boxplot %d are not drawn"), 
                        paste(unique(out[inf]), collapse = ", "), i), 
                domain = NA)
      }
    }
    return(ok)
  }
  if (!is.list(z) || 0L == (n <- length(z$n))) 
    stop("invalid first argument")
  if (is.null(at)) 
    at <- 1L:n
  else if (length(at) != n) 
    stop("'at' must have same length as 'z$n', i.e. ", n)
  if (is.null(z$out)) 
    z$out <- numeric()
  if (is.null(z$group) || !outline) 
    z$group <- integer()
  if (is.null(pars$ylim)) 
    ylim <- range(z$stats[is.finite(z$stats)], if (outline) z$out[is.finite(z$out)], 
                  if (notch) z$conf[is.finite(z$conf)])
  else {
    ylim <- pars$ylim
    pars$ylim <- NULL
  }
  if (length(border) == 0L) 
    border <- par("fg")
  dev.hold()
  on.exit(dev.flush())
  if (!add) {
    if (is.null(pars$xlim)) 
      xlim <- range(at, finite = TRUE) + c(-0.5, 0.5)
    else {
      xlim <- pars$xlim
      pars$xlim <- NULL
    }
    plot.new()
    if (horizontal) 
      plot.window(ylim = xlim, xlim = ylim, log = log, 
                  xaxs = pars$yaxs)
    else plot.window(xlim = xlim, ylim = ylim, log = log, 
                     yaxs = pars$yaxs)
  }
  xlog <- (par("ylog") && horizontal) || (par("xlog") && !horizontal)
  pcycle <- function(p, def1, def2 = NULL) rep(if (length(p)) p else if (length(def1)) def1 else def2, 
                                               length.out = n)
  p <- function(sym) pars[[sym, exact = TRUE]]
  boxlty <- pcycle(pars$boxlty, p("lty"), par("lty"))
  boxlwd <- pcycle(pars$boxlwd, p("lwd"), par("lwd"))
  boxcol <- pcycle(pars$boxcol, border)
  boxfill <- pcycle(pars$boxfill, bgcol)
  boxwex <- pcycle(pars$boxwex, 0.8 * {
    if (n <= 1) 
      1
    else stats::quantile(diff(sort(if (xlog) 
      log(at)
      else at)), 0.1)
  })
  medlty <- pcycle(pars$medlty, p("lty"), par("lty"))
  medlwd <- pcycle(pars$medlwd, 3 * p("lwd"), 3 * par("lwd"))
  medpch <- pcycle(pars$medpch, NA_integer_)
  medcex <- pcycle(pars$medcex, p("cex"), par("cex"))
  medcol <- pcycle(pars$medcol, border)
  medbg <- pcycle(pars$medbg, p("bg"), par("bg"))
  whisklty <- pcycle(pars$whisklty, p("lty"), "dashed")
  whisklwd <- pcycle(pars$whisklwd, p("lwd"), par("lwd"))
  whiskcol <- pcycle(pars$whiskcol, border)
  staplelty <- pcycle(pars$staplelty, p("lty"), par("lty"))
  staplelwd <- pcycle(pars$staplelwd, p("lwd"), par("lwd"))
  staplecol <- pcycle(pars$staplecol, border)
  staplewex <- pcycle(pars$staplewex, 0.5)
  outlty <- pcycle(pars$outlty, "blank")
  outlwd <- pcycle(pars$outlwd, p("lwd"), par("lwd"))
  outpch <- pcycle(pars$outpch, p("pch"), par("pch"))
  outcex <- pcycle(pars$outcex, p("cex"), par("cex"))
  outcol <- pcycle(pars$outcol, border)
  outbg <- pcycle(pars$outbg, p("bg"), par("bg"))
  outwex <- pcycle(pars$outwex, 0.5)
  width <- if (!is.null(width)) {
    if (length(width) != n | any(is.na(width)) | any(width <= 
                                                     0)) 
      stop("invalid boxplot widths")
    boxwex * width/max(width)
  }
  else if (varwidth) 
    boxwex * sqrt(z$n/max(z$n))
  else if (n == 1) 
    0.5 * boxwex
  else rep.int(boxwex, n)
  if (horizontal) {
    xypoints <- function(x, y, ...) points(y, x, ...)
    xypolygon <- function(x, y, ...) polygon(y, x, ...)
    xysegments <- function(x0, y0, x1, y1, ...) segments(y0, 
                                                         x0, y1, x1, ...)
  }
  else {
    xypoints <- points
    xypolygon <- polygon
    xysegments <- segments
  }
  ok <- TRUE
  for (i in 1L:n) ok <- ok & bplt(at[i], wid = width[i], stats = z$stats[, 
                                                                         i], out = z$out[z$group == i], conf = z$conf[, i], notch = notch, 
                                  xlog = xlog, i = i)
  if (!ok) 
    warning("some notches went outside hinges ('box'): maybe set notch=FALSE")
  axes <- is.null(pars$axes)
  if (!axes) {
    axes <- pars$axes
    pars$axes <- NULL
  }
  if (axes) {
    ax.pars <- pars[names(pars) %in% c("xaxt", "yaxt", "xaxp", 
                                       "yaxp", "las", "cex.axis", "col.axis", "format")]
    if (is.null(show.names)) 
      show.names <- n > 1
    if (show.names) 
      do.call("axis", c(list(side = 1 + horizontal, at = at, 
                             labels = z$names), ax.pars))
    do.call("Axis", c(list(x = z$stats, side = 2 - horizontal), 
                      ax.pars))
  }
  do.call("title", pars[names(pars) %in% c("main", "cex.main", 
                                           "col.main", "sub", "cex.sub", "col.sub", "xlab", "ylab", 
                                           "cex.lab", "col.lab")])
  if (frame.plot) 
    box()
  invisible(at)
}

interleave <- function(v1,v2){
  ord1 <- 2*(1:length(v1))-1
  ord2 <- 2*(1:length(v2))
  c(v1,v2)[order(c(ord1,ord2))]
}


deadwoodPlotNet <- function(net,metatable,ret=F,twowoods=F){
  coo <- data.frame("vertex"=V(net)$name,"xl"=-100,"xu"=100,"yl"=-100,"yu"=100,stringsAsFactors = F)
  dimw <- length(grep("_",coo$vertex))
  coo$xl[grep("_",coo$vertex)] <- jitter(6*cos(1:dimw*2*pi/dimw)[unlist(sapply(grep("_",coo$vertex,value=T), function(x) metatable$tworder1[metatable$treespecies_depth==x]))],1)
  coo$yl[grep("_",coo$vertex)] <- jitter(6*sin(1:dimw*2*pi/dimw)[unlist(sapply(grep("_",coo$vertex,value=T), function(x) metatable$tworder1[metatable$treespecies_depth==x]))],1)
  coo$xu[grep("_",coo$vertex)] <- coo$xl[grep("_",coo$vertex)]
  coo$yu[grep("_",coo$vertex)] <- coo$yl[grep("_",coo$vertex)]
  colay <-  layout_with_fr(net,minx=coo$xl,maxx=coo$xu,miny=coo$yl,maxy=coo$yu)
  plot(net,vertex.size=V(net)$size,vertex.label=NA,vertex.color=alpha(V(net)$color,0.5),vertex.frame.color=V(net)$frame.color,layout=colay,edge.color=alpha("grey50",0.15))
  
  if(twowoods){
    co2o <- data.frame("vertex"=V(net)$name,"xl"=-100,"xu"=100,"yl"=-100,"yu"=100,stringsAsFactors = F)
    co2o$xl[grep("_",co2o$vertex)] <- jitter(6*cos(1:dimw*2*pi/dimw)[unlist(sapply(grep("_",coo$vertex,value=T), function(x) metatable$tworder2[metatable$treespecies_depth==x]))],1)
    co2o$yl[grep("_",co2o$vertex)] <- jitter(6*sin(1:dimw*2*pi/dimw)[unlist(sapply(grep("_",coo$vertex,value=T), function(x) metatable$tworder2[metatable$treespecies_depth==x]))],1)
    co2o$xu[grep("_",co2o$vertex)] <- co2o$xl[grep("_",co2o$vertex)]
    co2o$yu[grep("_",co2o$vertex)] <- co2o$yl[grep("_",co2o$vertex)]
    co2lay <-  layout_with_fr(net,minx=co2o$xl,maxx=co2o$xu,miny=co2o$yl,maxy=co2o$yu)
    plot(net,vertex.size=V(net)$size,vertex.label=NA,vertex.color=alpha(V(net)$color,0.5),vertex.frame.color=V(net)$frame.color,layout=co2lay,edge.color=alpha("grey50",0.15))
    if(ret) return(list(colay,co2lay))
  }
  if(ret) return(colay)
}

deadwoodMakeIgraph <- function(binet,metatable,taxtable){
  require(igraph)
  require(bipartite)
  binettmp <- web2edges(binet,weight.column = F,return = T,is.one.mode = T)
  binettmp[,1] <- unlist(sapply(binettmp[,1],function(x) rownames(binet)[as.numeric(x)]))
  binettmp[,2] <- unlist(sapply(binettmp[,2],function(x) colnames(binet)[as.numeric(x)]))
  binet2 <- graph_from_edgelist(binettmp,directed = F)
  V(binet2)$type <- bipartite.mapping(binet2)$type
  V(binet2)$color <- unlist(sapply(V(binet2)$name,function(x){
    if(x %in% metatable$treespecies_depth) metatable$woodcol[metatable$treespecies_depth==x] else taxtable$phylcol[taxtable$OTU==x]
  }))
  V(binet2)$shape <- ifelse(V(binet2)$type, "circle", "square")
  V(binet2)$size <- ifelse(V(binet2)$type, 2, 3)
  V(binet2)$frame.color <- unlist(sapply(V(binet2)$name,function(x){
    if(x %in% metatable$treespecies_depth) metatable$woodframe[metatable$treespecies_depth==x] else taxtable$phylcol[taxtable$OTU==x]}
  ))
  return(binet2)
}