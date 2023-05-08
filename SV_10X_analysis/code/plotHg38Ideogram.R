
### modify SNPchip function "plotIdiogram"
plotIdiogram.hg38 <- function (chromosome, cytoband, seqinfo, cytoband.ycoords, xlim,
                               ylim = c(0, 2), new = TRUE, label.cytoband = TRUE, label.y = NULL,
                               srt, cex.axis = 1, outer = FALSE, taper = 0.15, verbose = FALSE,
                               unit = c("bp", "Mb"), is.lattice = FALSE, ...)
{
  def.par <- par(no.readonly = TRUE, mar = c(4.1, 0.1, 3.1,
                                             2.1))
  on.exit(def.par)
  if (is.lattice) {
    segments <- lsegments
    polygon <- lpolygon
  }
  
  cytoband <- cytoband[cytoband[, "chrom"] == chromosome, ]
  unit <- match.arg(unit)
  if (unit == "Mb") {
    cytoband$start <- cytoband$start/1e+06
    cytoband$end <- cytoband$end/1e+06
  }
  if (missing(cytoband.ycoords)) {
    cytoband.ycoords <- ylim
  }
  rownames(cytoband) <- as.character(cytoband[, "name"])
  sl <- seqlengths(seqinfo)[chromosome]
  if (missing(xlim)){
    xlim <- c(0, sl)}
  if (unit == "Mb"){
    xlim <- xlim/1e+06}
  cytoband_p <- cytoband[grep("^p", rownames(cytoband), value = TRUE),
                         ]
  cytoband_q <- cytoband[grep("^q", rownames(cytoband), value = TRUE),
                         ]
  p.bands <- nrow(cytoband_p)
  cut.left <- c()
  cut.right <- c()
  for (i in seq_len(nrow(cytoband))) {
    if (i == 1) {
      cut.left[i] <- TRUE
      cut.right[i] <- FALSE
    }
    else if (i == p.bands) {
      cut.left[i] <- FALSE
      cut.right[i] <- TRUE
    }
    else if (i == (p.bands + 1)) {
      cut.left[i] <- TRUE
      cut.right[i] <- FALSE
    }
    else if (i == nrow(cytoband)) {
      cut.left[i] <- FALSE
      cut.right[i] <- TRUE
    }
    else {
      cut.left[i] <- FALSE
      cut.right[i] <- FALSE
    }
  }
  for (i in seq_len(nrow(cytoband))) {
    if (as.character(cytoband[i, "gieStain"]) == "stalk") {
      cut.right[i - 1] <- TRUE
      cut.left[i] <- NA
      cut.right[i] <- NA
      cut.left[i + 1] <- TRUE
    }
  }
  include <- cytoband[, "end"] > xlim[1] & cytoband[, "start"] <
    xlim[2]
  cytoband <- cytoband[include, ]
  N <- nrow(cytoband)
  cut.left <- cut.left[include]
  cut.right <- cut.right[include]
  if (new) {
    xx <- c(0, cytoband[nrow(cytoband), "end"])
    yy <- cytoband.ycoords
    plot(xx, yy, xlim = xlim, type = "n", xlab = "", ylab = "",
         axes = FALSE, yaxs = "i", ylim = ylim, ...)
  }
  top <- cytoband.ycoords[2]
  bot <- cytoband.ycoords[1]
  h <- top - bot
  p <- taper
  for (i in seq_len(nrow(cytoband))) {
    start <- cytoband[i, "start"]
    last <- cytoband[i, "end"]
    delta = (last - start)/4
    getStain <- function(stain) {
      switch(stain, gneg = "grey100", gpos25 = "grey90",
             gpos50 = "grey70", gpos75 = "grey40", gpos100 = "grey0",
             gvar = "grey100", stalk = "brown3", acen = "brown4",
             "white")
    }
    color <- getStain(as.character(cytoband[i, "gieStain"]))
    if (is.na(cut.left[i]) & is.na(cut.right[i])) {
      delta <- (last - start)/3
      segments(start + delta, cytoband.ycoords[1], start +
                 delta, cytoband.ycoords[2])
      segments(last - delta, cytoband.ycoords[1], last -
                 delta, cytoband.ycoords[2])
    }
    else if (cut.left[i] & cut.right[i]) {
      yy <- c(bot + p * h, bot, bot, bot + p * h, top -
                p * h, top, top, top - p * h)
      polygon(c(start, start + delta, last - delta, last,
                last, last - delta, start + delta, start), yy,
              col = color)
    }
    else if (cut.left[i]) {
      yy <- c(bot + p * h, bot, bot, top, top, top - p *
                h)
      polygon(c(start, start + delta, last, last, start +
                  delta, start), yy, col = color)
    }
    else if (cut.right[i]) {
      yy <- c(bot, bot, bot + p * h, top - p * h, top,
              top)
      polygon(c(start, last - delta, last, last, last -
                  delta, start), yy, col = color)
    }
    else {
      polygon(c(start, last, last, start), c(bot, bot,
                                             top, top), col = color)
    }
  }
  my.x <- (cytoband[, "start"] + cytoband[, "end"])/2
  if (label.cytoband & !is.lattice) {
    if (is.null(label.y)) {
      axis(1, at = my.x, labels = rownames(cytoband), outer = outer,
           cex.axis = cex.axis, line = 1, las = 3, tick = FALSE)
      axis(1, at = cytoband$start, outer = outer, cex.axis = cex.axis,
           line = 1, las = 3, labels = FALSE)
    }
    else {
      if (!is.numeric(label.y)) {
        warning("label.y must be numeric -- using default y coordinates for cytoband labels")
        label.y <- bot - p * h
      }
      if (missing(srt))
        srt <- 90
      text(x = my.x, y = rep(label.y, length(my.x)), labels = rownames(cytoband),
           srt = srt)
    }
  }
  return()
}
