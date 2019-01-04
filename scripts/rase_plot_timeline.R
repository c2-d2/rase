#!/usr/bin/env Rscript

#
# Author: Karel Brinda <kbrinda@hsph.harvard.edu>
#
# License: MIT
#


library(optparse)


###############
# CLI PARSING #
###############

kIsRStudio <- Sys.getenv("RSTUDIO") == "1"

if (kIsRStudio) {
  src.file <- "../tests/spneumoniae_sparc.k18.predict.tsv"
} else {
  parser <-
    OptionParser(usage = "%prog [options] timeline.tsv plot.pdf")
  arguments <- parse_args(parser, positional_arguments = 2)

  opt <- arguments$options

  src.file <- arguments$args[1]
  out.file <- arguments$args[2]

  kWidth <- 4
  kHeight <- 10

  pdf(out.file,
      width = kWidth,
      height = kHeight)
}


#################
# CONFIGURATION #
#################

set.seed(42)

kRLUnitRatio <- 60

# first and second panels
kFirstMinutes <- 15
kLastHours <- 2

# remove endpoints more than ... far away
kEndpointFilter <- 3

# flag letters
kFlagCol <- "black"
kFlagSize <- 2.0

kLWD <- 2

kIndicSize <- 0.75
kIndicPos <- -1.2

kYLabDist <- 2.3

par(cex = 0.6)
par(oma = c(4, 4, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))



#############
# FUNCTIONS #
#############

LoadTimelineData <- function(src.file) {
  df <- read.delim(src.file, header = TRUE)
  df$datetime <-
    as.POSIXct(strptime(df$datetime, "%Y-%m-%d %H:%M:%S"))
  first.datetime <- df$datetime[1]

  df$time.mins <-
    difftime(df$datetime, first.datetime, units = "mins")

  # is it really sorted? should be..
  stopifnot(sort(df$datetime) == df$datetime)

  # remove too distant end points
  while (diff(tail(df, 2)$time) >= kEndpointFilter) {
    df <- head(df,-1)
  }

  df$inv.time.mins <- max(df$time.mins) - df$time.mins

  df
}

DfToAnts <- function(dataframe) {
  cols <- colnames(dataframe)
  antcols <- cols[grepl("_cat", cols)]
  ants <- gsub("_cat", "", antcols)
  ants
}

DfToFlags <- function(dataframe) {
  df.1 <- dataframe[grep("S:PG1", dataframe$flags), ]
  df.1$pch <- rep(4, nrow(df.1))

  df.2 <- dataframe[grep("S:PG2", dataframe$flags), ]
  df.2$pch <- rep(1, nrow(df.2))

  df.3 <- dataframe[grep("S:taxid", dataframe$flags), ]
  df.3$pch <- rep(0, nrow(df.3))

  df <- rbind(df.1, df.2, df.3)
  df
}

#
# plots vertical ablines
#
TimeAblines <- function(x) {
  abline(
    v = as.numeric(unlist(x)),
    lty = 2,
    col = "grey",
    lwd = 2
  )
}

#
# plots horizontal ablines
#
ThresholdAbline <- function(y) {
  abline(h = c(y),
         lty = 1,
         col = "grey")
}


#
# set margins for a subfigure
#
margin <- function(i) {
  if (i == 1) {
    par(mar = c(0, 0, 0, 0))
  }
  else{
    par(mar = c(0, 8, 0, 0))
  }
}


#
# plot res color boxes
#
RedBox <- function(df2, threshold) {
  mx <- max(df2$time.mins) + 15
  rect(-mx,-0.1,
       mx,
       threshold,
       col = rgb(1.0, 0, 0, alpha = 0.1),
       border = "NA")
}

GreenBox <- function(df2, threshold) {
  mx <- max(df2$time.mins) + 15
  rect(-mx,
       threshold,
       mx,
       1.1,
       col = rgb(0, 1.0, 0, alpha = 0.1),
       border = "NA")
}


#
# plots curve for nb of reads
#

PlotReads <- function(i) {
  margin(i)
  reads.ylim = c(0, max(df$read.count) / 1000)
  if (i == 1) {
    par(bty = "[")
    plot(
      df1$time.mins,
      df1$read.count / 1000,
      xlim = l.xlim,
      ylim = reads.ylim,
      type = 'l',
      las = 1,
      xaxt = 'n',
      ylab = NA,
      xlab = NA,
      lwd = kLWD,
      xaxs = "i"
    )
    points(
      df1.flag$time.mins,
      df1.flag$read.count / 1000,
      col = kFlagCol,
      pch = df1.flag$pch,
      cex = kFlagSize
    )
    mtext(
      "#reads (thousands)",
      side = 2,
      line = kYLabDist,
      cex.lab = 1,
      cex = 0.7,
      las = 3
    )
  }
  else {
    par(bty = "]")
    plot(
      df2$time.mins / kRLUnitRatio,
      df2$read.count / 1000,
      xlim = r.xlim,
      ylim = reads.ylim,
      type = 'l',
      las = 1,
      xaxt = 'n',
      yaxt = 'n',
      ylab = NA,
      xlab = NA,
      lwd = kLWD,
      xaxs = "i"
    )
    points(
      df2.flag$time.mins / kRLUnitRatio,
      df2.flag$read.count / 1000,
      col = kFlagCol,
      pch = df2.flag$pch,
      cex = kFlagSize
    )
  }

  TimeAblines(kVerticalAblines[i])

  if (i == 1) {
    legend("topleft",
           c("Predicted PG stabilized", "Alternative PG  stabilized", "Predicted isolate stabilized"),
           bg="white",
           pch = c(4, 1, 0))
  }
}


#
# plots curve for PG score
#
#df.flag=DfToFlags(df)
PlotPG <- function(i) {
  last_pg_predicted = tail(df, n = 1)["PG_score"] >= 0.6
  margin(i)
  if (i == 1) {
    par(bty = "[")
    plot(
      df1$time.mins,
      df1$PG_score,
      type = 'l',
      xlim = l.xlim,
      ylim = c(0, 1),
      las = 1,
      ylab = NA,
      xlab = NA,
      xaxt = 'n',
      lwd = kLWD,
      xaxs = "i"
    )

    mtext(
      "PG score",
      side = 2,
      line = kYLabDist,
      cex.lab = 1,
      cex = 0.7,
      las = 3
    )

    mtext(
      "fail",
      side = 2,
      line = kIndicPos,
      cex = kIndicSize,
      at = 0.3
    )

    mtext(
      "pass",
      side = 2,
      line = kIndicPos,
      cex = kIndicSize,
      at = 0.85
    )
  }
  else {
    par(bty = "]")
    plot(
      df2$time.mins / kRLUnitRatio,
      df2$PG_score,
      type = 'l',
      xlim = r.xlim,
      ylim = c(0, 1),
      xlab = NA,
      ylab = NA,
      las = 1,
      xaxt = 'n',
      yaxt = 'n',
      lwd = kLWD,
      xaxs = "i"
    )
  }
  ThresholdAbline(0.6)
  if (last_pg_predicted) {
    GreenBox(df2, 0.6)
  } else{
    RedBox(df2, 0.6)
  }
  TimeAblines(kVerticalAblines[i])
}


#
# plots curve for an antibiotic
#
PlotAntibiotic <- function(ant, i, is.last) {
  antcol <- paste(ant, "_susc_score", sep = "")
  print(paste(ant, antcol))

  last_is_resistant <- tail(df, n = 1)[antcol] <= 0.6

  par(bty = "l")
  margin(i)
  if (i == 1) {
    plot(
      df1$time.mins,
      df1[, antcol],
      type = 'l',
      xlim = l.xlim,
      ylim = c(0.0, 1),
      ylab = NA,
      xlab = NA,
      yaxt = "s",
      las = 1,
      xaxt = 'n',
      bty = "[",
      lwd = kLWD,
      xaxs = "i"
    )

    mtext(
      paste(ant, "susc score"),
      side = 2,
      line = kYLabDist,
      cex.lab = 1,
      cex = 0.7,
      las = 3
    )

    mtext(
      "non-susc",
      side = 2,
      line = kIndicPos,
      cex = kIndicSize,
      at = 0.6 / 2
    )
    mtext(
      "susc",
      side = 2,
      line = kIndicPos,
      cex = kIndicSize,
      at = 0.8
    )
  }
  else{
    plot(
      df2$time.mins / kRLUnitRatio,
      df2[, antcol],
      type = 'l',
      xlim = r.xlim,
      ylim = c(0.0, 1),
      yaxt = "n",
      xlab = NA,
      ylab = NA,
      las = 1,
      xaxt = 'n',
      bty = "]",
      lwd = kLWD,
      xaxs = "i"
    )
  }

  if (last_is_resistant) {
    RedBox(df2, 0.6)
  } else{
    GreenBox(df2, 0.6)
  }

  ThresholdAbline(0.6)
  TimeAblines(kVerticalAblines[i])


  # last row => plot labels
  if (is.last) {
    if (i == 1) {
      axis(1, lwd = 0.5)

      mtext("minutes",
            side = 1,
            line = 2)

    }
    else {
      axis(1, lwd = 0.5)

      mtext("hours",
            side = 1,
            line = 2)

    }
  }
}


############
# PLOTTING #
############

df <- LoadTimelineData(src.file)

df1 <- df[df$time.min <= kFirstMinutes,]
df2 <- df[df$inv.time.min <= kRLUnitRatio * kLastHours,]

df1.flag <- DfToFlags(df1)
df2.flag <- DfToFlags(df2)

last.min <- max(df$time.mins)

l.xlim1 <- -0.10 * kFirstMinutes # left padding
l.xlim2 <- kFirstMinutes
l.xlim <- c(l.xlim1, l.xlim2)

r.xlim1 <-
  max(l.xlim[2] / kRLUnitRatio, last.min / kRLUnitRatio - kLastHours)
r.xlim2 <- last.min / kRLUnitRatio
r.xlim2 <- r.xlim2 + (r.xlim2 - r.xlim1) * 0.05 # right padding
r.xlim <- c(r.xlim1, r.xlim2)

kVerticalAblines <-
  list(A = c(1.0, 5.0), B = c(last.min / kRLUnitRatio))

ants <- DfToAnts(df)

par(mfrow = c(length(ants) + 2, 2), tcl = -0.5)


################
# READS
################

for (i in c(1, 2)) {
  PlotReads(i)
}

################
# PG SCORE
################

for (i in c(1, 2)) {
  PlotPG(i)
}

################
# RESISTANCE
################

last.ant <- tail(ants, 1)
for (ant in ants) {
  is.last <- ant == last.ant
  for (i in c(1, 2)) {
    PlotAntibiotic(ant, i, is.last)
  }
}

if (!kIsRStudio) {
  dev.off()
}