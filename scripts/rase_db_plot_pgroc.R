#!/usr/bin/env Rscript

#
# Author:  Karel Brinda <kbrinda@hsph.harvard.edu>
#
# License: MIT
#

library(RColorBrewer)
library(data.table)

set.seed(42)

getAnts <- function(df1) {
  c = colnames(df1)
  m = grepl(".*_Rr" , c, perl = T)
  c2 = c[m]
  c3 = sub("_Rr", "", c2)
  c3
}

tableCases <- function(df1, ant, norm = NA) {
  # Load table with phylogroup cases. Normalize to norm isolates per phylogroup.
  print("")
  print(paste("-------", ant, "-------"))
  print("")
  df2 = subset(df1, select = c("pg", "count", paste0(ant, "_Rr"), paste0(ant, "_Ss")))
  print(df2)
  colnames(df2)[3] = "R"
  colnames(df2)[4] = "S"

  if (!is.na(norm)) {
    df2$R = norm * df2$R / df2$count
    df2$S = norm * df2$S / df2$count
    df2$count = df2$S + df2$R
  }

  df2$prop = df2$R / df2$count
  print(df2)

  df2
}

tableCumul <- function(df1) {
  df2 = df1[order(-df1[, "prop"]), ]
  df3 = subset(df2, select = c("pg", "R", "S"))
  df4 = rbind(c(0, 0), df3)
  df5 = cumsum(df4)
  df5$pg = c(0, df2$pg)
  rownames(df5) <- NULL
  print(df5)
  df5
}

tableRoc <- function(df1) {
  setnames(df1, old = "R", new = "TP")
  setnames(df1, old = "S", new = "FP")
  df1$FN = max(df1$TP) - df1$TP
  df1$TN = max(df1$FP) - df1$FP
  print(df1)
  df1
}

maxResIsolates <- function(df, ants, norm) {
  res = c()
  for (ant in ants) {
    dfc = tableCases(df, ant, norm)
    res = c(res, sum(dfc$R))
  }
  max(res)
}

#df <- read.delim(fn)
all.ants <- scan("antibiotics.txt",
                 what = "",
                 sep = "\n",
                 quiet = T)

plotToc = function(df,
                   all.ants,
                   main = NA,
                   norm = NA,
                   y.max = NA,
                   xlab = "autocomplete",
                   ylab = "autocomplete") {
  ants <- getAnts(df)
  l = length(all.ants)
  palette(brewer.pal(l, "Paired"))

  if (identical(xlab, "autocomplete")) {
    if (is.na(norm)) {
      xlab = "no. of isolates"
    }
    else {
      xlab = "no. of phylogroups"
    }
  }

  if (identical(ylab, "autocomplete")) {
    if (is.na(norm)) {
      ylab = "cumulative no. of resistant isolates"
    } else {
      ylab = "cumulative no. of fully resistant phylogroups"
    }
  }


  if (is.na(norm)) {
    x.max = sum(df$count)
  }
  else {
    x.max = norm * length(df[, 1])
  }

  if (is.na(y.max)) {
    y.max = 1.2 * maxResIsolates(df, ants, norm)
  }

  plot(
    NA,
    xlim = c(0, x.max),
    ylim = c(0, y.max),
    xlab = xlab,
    ylab = ylab,
    main = main
    #axes = F
  )

  i = 1
  mask = c()
  for (a in all.ants) {
    if (a %in% ants) {
      mask = c(mask, T)

      dfc = tableCases(df, a, norm)
      dfl = tableCumul(dfc)
      dfr = tableRoc(dfl)

      no.resistant.isolates = sum(dfc$R)

      points(dfr$TP + dfr$FP,
             dfr$TP,
             #pch = NA,
             pch = i,
             col = i)
      lines(dfr$TP + dfr$FP,
            dfr$TP,
            col = i)

      m = max(dfr$TP)
      lines(c(m, x.max), c(m, m), lt = 6)

    }
    else{
      mask = c(mask, F)

    }
    i = i + 1
  }

  lines(c(0, y.max), c(0, y.max), lt = 2)
  lines(c(0, x.max), c(0, 0), lt = 1)

  legend(
    "topleft",
    legend = c(all.ants[mask], "fully res.", "fully sus.", "best"),
    col = c(seq(l)[mask], "black", "black", "black"),
    pch = c(seq(l)[mask], NA, NA, NA),
    lt = c(rep(1, l)[mask], 2, 3, 6),
    bg = "white"
  )
}

plotRoc = function(df,
                   all.ants,
                   norm = NA,
                   xlab = "FPR (1-specificity)",
                   ylab = "TPR (sensitivity)",
                   xlim = c(0, 1),
                   ylim = c(0, 1),
                   ...) {
  ants <- getAnts(df)
  l = length(all.ants)
  palette(brewer.pal(l, "Paired"))

  plot(
    NA,
    xlim = xlim,
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    ...
  )

  lines(c(0, 1), c(0, 1), lt = 2)
  lines(c(0, 0), c(0, 1), lt = 6)
  lines(c(0, 1), c(1, 1), lt = 6)

  i = 1
  mask = c()
  for (a in all.ants) {
    if (a %in% ants) {
      mask = c(mask, T)

      dfc = tableCases(df, a, norm)
      dfl = tableCumul(dfc)
      dfr = tableRoc(dfl)

      points(
        dfr$FP / (dfr$FP + dfr$TN),
        dfr$TP / (dfr$TP + dfr$FN),
        pch = i,
        col = i,
        type = "b"
      )
    }
    else{
      mask = c(mask, F)

    }
    i = i + 1
  }

  legend(
    "bottomright",
    legend = c(all.ants[mask], "random", "best"),
    col = c(seq(l)[mask], "black", "black"),
    pch = c(seq(l)[mask], NA, NA),
    lt = c(rep(1, l)[mask], 2, 6)
  )
}

integrateRoc <- function(fpr, tpr) {
  x = fpr
  y = tpr
  delta.x = diff(x)
  avg.y = (head(y, n = -1) + tail(y, n = -1)) / 2
  sum(delta.x * avg.y)
}

ants.to.aucs <- function(df, all.ants = NA, norm = NA)
{
  #area under roc curve for all antibiotics
  ants <- getAnts(df)
  if (identical(all.ants,NA)) {
    all.ants = ants
  }
  aucs = c()
  i = 1
  cols = c()
  for (a in all.ants) {
    if (a %in% ants) {
      cols = c(cols, i)
      dfc = tableCases(df, a, norm)
      dfl = tableCumul(dfc)
      dfr = tableRoc(dfl)

      auc = integrateRoc(dfr$FP / (dfr$FP + dfr$TN),
                         dfr$TP / (dfr$TP + dfr$FN))
      aucs = c(aucs, auc)
    }
    i = i + 1
  }
  dfa = data.frame(v1 = ants, v2 = aucs, v3 = cols)
  dfa
}

plotAuc <- function(df,
                    all.ants = NA,
                    norm = NA,
                    ylim = c(0.5, 1),
                    annotate = F,
                    ...)
{
  ants <- getAnts(df)
  aucs = ants.to.aucs(df, all.ants, norm)
  x = barplot(
    aucs$v2,
    names.arg = aucs$v1,
    col = aucs$v3,
    ylim = ylim,
    xpd = F,
    ...
  )
  abline(h = seq(0.5, 1.0, 0.1), lt = 2)
  offset = 1.4
  if (annotate) {
    text(
      length(aucs[, 1]) + offset,
      seq(.55, .95, 0.10),
      c("fail", "poor", "fair", "good", "excellent"),
      pos=2
    )
  }
}

###########
# DEPREC
###########

plot.yonlike <-
  function(fn, xlab = "no. of phylogroups", ylab = "cumulative fraction of total resistant isolates", ...) {
    df <- read.delim(fn)

    df2 = transpose(df)
    l = length(colnames(df2))
    colnames(df2) = df2[1, ] # the first row will be the header

    df2 = df2[-1, ]          # removing the first row.

    cols.act = colSums(is.na(df2)) == 0 # active columns, to be plotted


    palette(brewer.pal(l, "Paired"))
    matplot(
      df2,
      xlab = xlab,
      ylab = ylab,
      type = "b",
      col = seq(l),
      pch = seq(l),
      lty = 1,
      axes = F,
      ...
    )

    axis(
      side = 1,
      at = 1:nrow(df2),
      labels = 0:(nrow(df2) - 1)
    )
    axis(2)
    axis(
      side = 3,
      at = 1:nrow(df2),
      labels = F,
      col.ticks = NA
    )
    axis(4, labels = F, col.ticks = NA)

    legend(
      "bottomright",
      legend = colnames(df2)[cols.act],
      col = seq(l)[cols.act],
      pch = seq(l)[cols.act]
    )
  }
