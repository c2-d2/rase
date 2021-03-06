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

ls.thres.pass <- 0.5
ss.thres.shiconf <- 0.6
ss.thres.sr <- 0.5
ss.thres.rhiconf <- 0.4

if (kIsRStudio) {
    src.file <- "tests/predict.tsv"
} else {
    option_list <- list(
        make_option(
            c("--ls-thres-pass"),
            dest = "ls.thres.pass",
            default = ls.thres.pass,
            help = "lineage score threshold for high-confidence calls [default %default]",
            metavar = "FLOAT"
        ),
        make_option(
            c("--ss-thres-shiconf"),
            dest = "ss.thres.shiconf",
            default = ss.thres.shiconf,
            help = "susceptibility score threshold for high-confidence susceptible [default %default]",
            metavar = "FLOAT"
        ),
        make_option(
            c("--ss-thres-sr"),
            dest = "ss.thres.sr",
            default = ss.thres.sr,
            help = "susceptibility score threshold for susceptible/non-susceptible [default %default]",
            metavar = "FLOAT"
        ),
        make_option(
            c("--ss-thres-rhiconf"),
            dest = "ss.thres.rhiconf",
            default = ss.thres.rhiconf,
            help = "susceptibility score threshold for high-confidence non-susceptible [default %default]",
            metavar = "FLOAT"
        )
    )

    parser <-
        OptionParser(usage = "%prog [options] prediction.tsv timeline.pdf", option_list =
                option_list)
    arguments <- parse_args(parser, positional_arguments = 2)

    opt <- arguments$options

    ls.thres.pass <- opt$ls.thres.pass
    ss.thres.shiconf <- opt$ss.thres.shiconf
    ss.thres.sr <- opt$ss.thres.sr
    ss.thres.rhiconf <- opt$ss.thres.rhiconf

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
    antcols <- cols[grepl("_ss", cols)]
    ants <- gsub("_ss", "", antcols)
    ants
}


#' Create a data frame with flag points
#'
#' @param df Input dataframe
#'
#' @return Dataframe with flag points to plot
#'
DfToFlags <- function(dataframe) {
    df.1 <- dataframe[grep("S:lineage", dataframe$flags),]
    df.1$pch <- rep(4, nrow(df.1))

    #df.2 <- dataframe[grep("S:alt_lineage", dataframe$flags), ]
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
#' @param y SS thresholds
#'
ThresholdAbline <- function(y, lty = 1, col = "grey", ...) {
    abline(h = c(y),
        lty = lty,
        col = col,
        ...)
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
PlotReads <- function(df,
    df.flag,
    i,
    ylim,
    xlab = NA,
    ylab = NA) {
    margin(i)
    if (i == 1) {
        par(bty = "[")
        plot(
            df$time.mins,
            df$reads / 1000,
            xlim = l.xlim,
            ylim = ylim,
            type = "l",
            las = 1,
            xaxt = "n",
            ylab = ylab,
            xlab = xlab,
            lwd = kLWD,
            xaxs = "i"
        )
        mtext(
            "#reads (thousands)",
            side = 2,
            line = kYLabDist,
            cex.lab = 1,
            cex = 0.7,
            las = 3
        )
        points(
            df.flag$time.mins,
            df.flag$reads / 1000,
            col = kFlagCol,
            pch = df.flag$pch,
            cex = kFlagSize
        )
    } else {
        par(bty = "]")
        plot(
            df$time.mins / kRLUnitRatio,
            df$reads / 1000,
            xlim = r.xlim,
            ylim = ylim,
            type = "l",
            las = 1,
            xaxt = "n",
            yaxt = "n",
            ylab = ylab,
            xlab = xlab,
            lwd = kLWD,
            xaxs = "i"
        )
        points(
            df.flag$time.mins / kRLUnitRatio,
            df.flag$reads / 1000,
            col = kFlagCol,
            pch = df.flag$pch,
            cex = kFlagSize
        )
    }

    TimeAblines(kVerticalAblines[i])

    if (i == 1) {
        legend(
            "topleft",
            c(
                "Predicted lineage stabilized",
                #"Alternative lineage  stabilized",
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
PlotKS <-
    function(df,
        i,
        ylim,
        type = "l",
        xlab = NA,
        ylab = NA,
        xaxt = "n",
        las = 1,
        lwd = kLWD,
        xaxs = "i",
        ...)
    {
        margin(i)
        if (i == 1) {
            par(bty = "[")
            plot(
                df$time.mins,
                df$ks,
                xlim = l.xlim,
                ylim = ylim,
                type = type,
                las = las,
                xaxt = xaxt,
                xlab = xlab,
                ylab = ylab,
                lwd = lwd,
                xaxs = xaxs,
                ...
            )
            mtext(
                "KS",
                side = 2,
                line = kYLabDist,
                cex.lab = 1,
                cex = 0.7,
                las = 3
            )
        } else {
            par(bty = "]")
            plot(
                df$time.mins / kRLUnitRatio,
                df$ks,
                xlim = r.xlim,
                ylim = ylim,
                type = type,
                las = las,
                xaxt = xaxt,
                yaxt = "n",
                xlab = xlab,
                ylab = ylab,
                lwd = lwd,
                xaxs = xaxs,
                ...
            )
        }

        TimeAblines(kVerticalAblines[i])
    }


#' Plot LS
#'
#' @param i
#'
#' @return
#' @export
#'
#' @examples
PlotLS <- function(df, i) {
    last_lineage_predicted <- tail(df, n = 1)["ls"] >= ls.thres.pass
    margin(i)
    if (i == 1) {
        par(bty = "[")
        plot(
            df$time.mins,
            df$ls,
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
            "LS",
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
            df$time.mins / kRLUnitRatio,
            df$ls,
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
    ThresholdAbline(ls.thres.pass)
    if (last_lineage_predicted) {
        GreenBox(df2, ls.thres.pass)
    } else {
        RedBox(df2, ls.thres.pass)
    }
    TimeAblines(kVerticalAblines[i])
}


#' Plot SS
#'
#' @param ant
#' @param i
#' @param is.last
#'
#' @return
#' @export
#'
#' @examples
PlotAntibiotic <- function(df, ant, i, is.last) {
    antcol <- paste0(ant, "_ss")
    print(antcol)
    last_is_resistant <- tail(df, n = 1)[antcol] <= ss.thres.sr
    par(bty = "l")
    margin(i)
    if (i == 1) {
        plot(
            df$time.mins,
            df[, antcol],
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
            paste(toupper(ant), "SS"),
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
            at = 0.22
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
            df$time.mins / kRLUnitRatio,
            df[, antcol],
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
        RedBox(df2, ss.thres.sr)
    } else {
        GreenBox(df2, ss.thres.sr)
    }

    ThresholdAbline(ss.thres.rhiconf, lty = 3)
    ThresholdAbline(ss.thres.shiconf, lty = 3)
    ThresholdAbline(ss.thres.sr)
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
reads.ylim <- c(0, max(df2$reads) / 1000)
PlotReads(df1, df1.flag, 1, ylim=reads.ylim)
PlotReads(df2, df2.flag, 2, ylim=reads.ylim)

# 2) ks
y.lim <- c(0, max(pretty(c(
    df1$ks, df2$ks
))))
PlotKS(df1, 1, ylim = y.lim)
PlotKS(df2, 2, ylim = y.lim)

# 3) ls
PlotLS(df1, 1)
PlotLS(df2, 2)

# 4) ss
last.ant <- tail(ants, 1)
for (ant in ants) {
    is.last <- ant == last.ant
    PlotAntibiotic(df1, ant, 1, is.last)
    PlotAntibiotic(df2, ant, 2, is.last)
}

if (!kIsRStudio) {
    dev.off()
}
