#!/usr/bin/env Rscript

#' Plot RASE prediction timeline.
#'
#' Author: Karel Brinda <kbrinda@hsph.harvard.edu>
#'
#' License: MIT
#'

suppressMessages(suppressWarnings(library(optparse)))


# CLI-parsing -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

kIsRStudio <- Sys.getenv("RSTUDIO") == "1"

pgs.thres.pass <- 0.5
ssc.thres.shiconf <- 0.6
ssc.thres.sr <- 0.5
ssc.thres.rhiconf <- 0.4

if (kIsRStudio) {
    src.file <- "pipeline/tests/predict.tsv"
} else {
    option_list <- list(
        make_option(
            c("--pgs-thres-pass"),
            dest = "pgs.thres.pass",
            default = pgs.thres.pass,
            help = "phylogroup score threshold [default %default]",
            metavar = "FLOAT"
        ),
        make_option(
            c("--ssc-thres-shiconf"),
            dest = "ssc.thres.shiconf",
            default = ssc.thres.shiconf,
            help = "phylogroup score threshold [default %default]",
            metavar = "FLOAT"
        ),
        make_option(
            c("--ssc-thres-sr"),
            dest = "ssc.thres.sr",
            default = ssc.thres.sr,
            help = "phylogroup score threshold [default %default]",
            metavar = "FLOAT"
        ),
        make_option(
            c("--ssc-thres-rhiconf"),
            dest = "ssc.thres.rhiconf",
            default = ssc.thres.rhiconf,
            help = "phylogroup score threshold [default %default]",
            metavar = "FLOAT"
        )
    )

    parser <-
        OptionParser(usage = "%prog [options] prediction.tsv timeline.pdf", option_list =
                         option_list)
    arguments <- parse_args(parser, positional_arguments = 2)

    opt <- arguments$options

	pgs.thres.pass <- opt$pgs.thres.pass
	ssc.thres.shiconf <- opt$ssc.thres.shiconf
	ssc.thres.sr <- opt$ssc.thres.sr
	ssc.thres.rhiconf <- opt$ssc.thres.rhiconf

    src.file <- arguments$args[1]
    out.file <- arguments$args[2]

    kWidth <- 4
    kHeight <- 10

    pdf(out.file, width = kWidth, height = kHeight)
}



# CONFIGURATION -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

set.seed(42)

kRLUnitRatio <- 60

# first and second panels
kFirstMinutes <- 15
kLastHours <- 2

# remove endpoints more than ... far away
kEndpointFilter <- 3

# flag letters
kFlagCol <- "black"
kFlagSize <- 2

kLWD <- 2

kIndicSize <- 0.75
kIndicPos <- -1.2

kYLabDist <- 2.3

par(cex = 0.6)
par(oma = c(4, 4, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))



# FUNCTIONS-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Load output of RASE prediction
#'
#' @param src.file TSV table with the output
#'
#' @return Data frame with loaded prediction data
#'
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
    if (length(df$time) > 5) {
        while (diff(tail(df, 2)$time) >= kEndpointFilter) {
            df <- head(df,-1)
        }
    }

    df$inv.time.mins <- max(df$time.mins) - df$time.mins

    df
}


#' Extract antibiotics from a data frame
#'
#' @param df Input dataframe
#'
#' @return List of antibiotics
#'
DfToAnts <- function(df) {
    cols <- colnames(df)
    antcols <- cols[grepl("sus", cols)]
    ants <- gsub("_sus", "", antcols)
    ants
}


#' Create a data frame with flag points
#'
#' @param df Input dataframe
#'
#' @return Dataframe with flag points to plot
#'
DfToFlags <- function(dataframe) {
    df.1 <- dataframe[grep("S:pg", dataframe$flags),]
    df.1$pch <- rep(4, nrow(df.1))

    #df.2 <- dataframe[grep("S:alt_pg", dataframe$flags), ]
    #df.2$pch <- rep(1, nrow(df.2))

    df.3 <- dataframe[grep("S:bm", dataframe$flags),]
    df.3$pch <- rep(0, nrow(df.3))

    df <- rbind(df.1, df.3)
    df
}


#' Plot vertical ablines (displayed snapshots)
#'
#' @param x Snapshots (time)
#'
TimeAblines <- function(x) {
    abline(
        v = as.numeric(unlist(x)),
        lty = 2,
        col = "grey",
        lwd = 2
    )
}


#' Plots horizontal ablines
#'
#' @param y SSC thresholds
#'
ThresholdAbline <- function(y) {
    abline(h = c(y),
           lty = 1,
           col = "grey")
}


#' Set margins for a subfigure
#'
#' @param i Position: 1=left, 2=right
#'
margin <- function(i) {
    if (i == 1) {
        par(mar = c(0, 0, 0, 0))
    } else {
        par(mar = c(0, 8, 0, 0))
    }
}


#' Plot resistance color boxes (if predicted as resistance)
#'
#' @param df2 Dataframe
#' @param threshold Threshold of resistance
#'
RedBox <- function(df2, threshold) {
    mx <- max(df2$time.mins) + 15
    rect(-mx,
         -0.1,
         mx,
         threshold,
         col = rgb(1, 0, 0, alpha = 0.1),
         border = "NA")
}

GreenBox <- function(df2, threshold) {
    mx <- max(df2$time.mins) + 15
    rect(-mx,
         threshold,
         mx,
         1.1,
         col = rgb(0, 1, 0, alpha = 0.1),
         border = "NA")
}


#' Plot nb of reads
#'
#' @param i Position: 1=left, 2=right
#'
PlotReads <- function(df1, i) {
    margin(i)
    reads.ylim <- c(0, max(df$reads) / 1000)
    if (i == 1) {
        par(bty = "[")
        plot(
            df1$time.mins,
            df1$reads / 1000,
            xlim = l.xlim,
            ylim = reads.ylim,
            type = "l",
            las = 1,
            xaxt = "n",
            ylab = NA,
            xlab = NA,
            lwd = kLWD,
            xaxs = "i"
        )
        points(
            df1.flag$time.mins,
            df1.flag$reads / 1000,
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
    } else {
        par(bty = "]")
        plot(
            df2$time.mins / kRLUnitRatio,
            df2$reads / 1000,
            xlim = r.xlim,
            ylim = reads.ylim,
            type = "l",
            las = 1,
            xaxt = "n",
            yaxt = "n",
            ylab = NA,
            xlab = NA,
            lwd = kLWD,
            xaxs = "i"
        )
        points(
            df2.flag$time.mins / kRLUnitRatio,
            df2.flag$reads / 1000,
            col = kFlagCol,
            pch = df2.flag$pch,
            cex = kFlagSize
        )
    }

    TimeAblines(kVerticalAblines[i])

    if (i == 1) {
        legend(
            "topleft",
            c(
                "Predicted PG stabilized",
                #"Alternative PG  stabilized",
                "Predicted isolate stabilized"
            ),
            bg = "white",
            pch = c(4,
                    #1,
                    0)
        )
    }
}


#' Plot nb of reads
#'
#' @param i Position: 1=left, 2=right
#'
PlotProportion <- function(df1, i) {
    margin(i)
    prop.ylim <- c(0, 0.10)
    if (i == 1) {
        par(bty = "[")
        plot(
            df1$time.mins,
            df1$matched_prop,
            xlim = l.xlim,
            ylim = prop.ylim,
            type = "l",
            las = 1,
            xaxt = "n",
            ylab = NA,
            xlab = NA,
            lwd = kLWD,
            xaxs = "i"
        )
        mtext(
            "proportion of k-mers",
            side = 2,
            line = kYLabDist,
            cex.lab = 1,
            cex = 0.7,
            las = 3
        )
    } else {
        par(bty = "]")
        plot(
            df2$time.mins / kRLUnitRatio,
            df2$matched_prop,
            xlim = r.xlim,
            ylim = prop.ylim,
            type = "l",
            las = 1,
            xaxt = "n",
            yaxt = "n",
            ylab = NA,
            xlab = NA,
            lwd = kLWD,
            xaxs = "i"
        )
    }

    TimeAblines(kVerticalAblines[i])
}


#' Plot PGS
#'
#' @param i
#'
#' @return
#' @export
#'
#' @examples
PlotPGS <- function(i) {
    last_pg_predicted <- tail(df, n = 1)["pgs"] >= pgs.thres.pass
    margin(i)
    if (i == 1) {
        par(bty = "[")
        plot(
            df1$time.mins,
            df1$pgs,
            type = "l",
            xlim = l.xlim,
            ylim = c(0, 1),
            las = 1,
            ylab = NA,
            xlab = NA,
            xaxt = "n",
            lwd = kLWD,
            xaxs = "i"
        )

        mtext(
            "PGS",
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
    } else {
        par(bty = "]")
        plot(
            df2$time.mins / kRLUnitRatio,
            df2$pgs,
            type = "l",
            xlim = r.xlim,
            ylim = c(0, 1),
            xlab = NA,
            ylab = NA,
            las = 1,
            xaxt = "n",
            yaxt = "n",
            lwd = kLWD,
            xaxs = "i"
        )
    }
    ThresholdAbline(pgs.thres.pass)
    if (last_pg_predicted) {
        GreenBox(df2, pgs.thres.pass)
    } else {
        RedBox(df2, pgs.thres.pass)
    }
    TimeAblines(kVerticalAblines[i])
}


#' Plot SSC
#'
#' @param ant
#' @param i
#' @param is.last
#'
#' @return
#' @export
#'
#' @examples
PlotAntibiotic <- function(ant, i, is.last) {
    antcol <- paste0(ant, "_sus")
    last_is_resistant <- tail(df, n = 1)[antcol] <= ssc.thres.sr
    par(bty = "l")
    margin(i)
    if (i == 1) {
        plot(
            df1$time.mins,
            df1[, antcol],
            type = "l",
            xlim = l.xlim,
            ylim = c(0, 1),
            ylab = NA,
            xlab = NA,
            yaxt = "s",
            las = 1,
            xaxt = "n",
            bty = "[",
            lwd = kLWD,
            xaxs = "i"
        )

        mtext(
            paste(toupper(ant), "SSC"),
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
    } else {
        plot(
            df2$time.mins / kRLUnitRatio,
            df2[, antcol],
            type = "l",
            xlim = r.xlim,
            ylim = c(0, 1),
            yaxt = "n",
            xlab = NA,
            ylab = NA,
            las = 1,
            xaxt = "n",
            bty = "]",
            lwd = kLWD,
            xaxs = "i"
        )
    }

    if (last_is_resistant) {
        RedBox(df2, ssc.thres.sr)
    } else {
        GreenBox(df2, ssc.thres.sr)
    }

    ThresholdAbline(ssc.thres.rhiconf)
    ThresholdAbline(ssc.thres.sr)
    ThresholdAbline(ssc.thres.shiconf)
    TimeAblines(kVerticalAblines[i])


    # last row => plot labels
    if (is.last) {
        if (i == 1) {
            axis(1, lwd = 0.5)

            mtext("minutes", side = 1, line = 2)

        } else {
            axis(1, lwd = 0.5)

            mtext("hours", side = 1, line = 2)

        }
    }
}


# PLOTTING -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

df <- LoadTimelineData(src.file)

df1 <- df[df$time.min <= kFirstMinutes,]
df2 <- df[df$inv.time.min <= kRLUnitRatio * kLastHours,]

df1.flag <- DfToFlags(df1)
df2.flag <- DfToFlags(df2)

last.min <- max(df$time.mins)

l.xlim1 <- -0.1 * kFirstMinutes  # left padding
l.xlim2 <- kFirstMinutes
l.xlim <- c(l.xlim1, l.xlim2)

r.xlim1 <-
    max(l.xlim[2] / kRLUnitRatio, last.min / kRLUnitRatio - kLastHours)
r.xlim2 <- last.min / kRLUnitRatio
r.xlim2 <- r.xlim2 + (r.xlim2 - r.xlim1) * 0.05  # right padding
r.xlim <- c(r.xlim1, r.xlim2)

kVerticalAblines <-
    list(A = c(1, 5), B = c(last.min / kRLUnitRatio))

ants <- DfToAnts(df)

par(mfrow = c(length(ants) + 3, 2), tcl = -0.5)

# 1) reads
PlotReads(df1, 1)
PlotReads(df2, 2)

# 2) proportion of k-mers
PlotProportion(df1, 1)
PlotProportion(df2, 2)

# 3) pgs
PlotPGS(1)
PlotPGS(2)

# 4) ssc
last.ant <- tail(ants, 1)
for (ant in ants) {
    is.last <- ant == last.ant
    PlotAntibiotic(ant, 1, is.last)
    PlotAntibiotic(ant, 2, is.last)
}

if (!kIsRStudio) {
    dev.off()
}
