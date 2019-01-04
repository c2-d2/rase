#!/usr/bin/env Rscript

#
# Author: Karel Brinda <kbrinda@hsph.harvard.edu>
#
# License: MIT
#

library(optparse)


#################
# CONFIGURATION #
#################

set.seed(42)

kWidth = 8
kHeight = 7

kSelected <- 70 # strains to display

kResGridHeight <-
  0.045 # heigth of a resistance cell, fraction of the height of the graph above
kGridShift <- 2

# c(
#   "#1b9e77",
#   "#d95f02",
#   "#7570b3",
#   "#e7298a",
#   "#66a61e",
#   "#e6ab02",
#   "#a6761d",
#   "#666666"
# ),


kPalette <- rep(c("#d95f02",
                  "#7570b3",
                  "#222222",
                  "#e7298a"),
                times = 2)

kRStudio <- Sys.getenv("RSTUDIO") == "1"

# R / S / NA
kCatColors <- c("#ff0000", "#0000aa", '#00aa00')


#############
# FUNCTIONS #
#############

CatToColor <- function(cat) {
  catu = toupper(cat)
  ifelse(catu == "R", kCatColors[1], (ifelse(catu == "S", kCatColors[2], kCatColors[3])))
}

CatToNumber <- function(cat) {
  catu = toupper(cat)
  ifelse(catu == "R", 1, (ifelse(catu == "S", 2, 3)))
}

CatToCatName <- function(cat) {
  catu = toupper(cat)
  ifelse(catu == "S", "susceptible", (ifelse(
    catu ==
      "R", "non-susceptible", "unknown"
  )))
}

DfToAnts <- function(dataframe) {
  cols = colnames(dataframe)
  antcols = cols[grepl("_cat", cols)]
  ants = gsub("_cat", "", antcols)
  ants = ants[ants != "CHL"]
  ants = ants[ants != "chl"]
  ants
}


###############
# CLI PARSING #
###############

if (kRStudio) {
  src.file <- "../tests/test.timestamp.tsv"
  res.file <- "../tests/res_cat.tsv"
  timestamp <- 1505967676
  sample.desc <- "Sample desc."
} else {
  option.list <- list(make_option(c("-d", "--desc"),
                                  default = F,
                                  help = "sample description"))
  parser <-
    OptionParser(usage = "%prog [options] res_cats.tsv snapshot.tsv plot.pdf", option_list = option.list)

  arguments <- parse_args(parser, positional_arguments = 3)
  opt <- arguments$options

  sample.desc <- ""
  if (opt$desc) {
    sample.desc <- opt$desc
  }

  res.file <- arguments$args[1]
  src.file <- arguments$args[2]
  out.file <- arguments$args[3]

  timestamp <- as.integer(rev(strsplit(src.file, '[\\./]')[[1]])[2])

  pdf(out.file,
      width = kWidth,
      height = kHeight)

}



############
# PLOTTING #
############

palette(kPalette)

dfres <- read.delim(res.file, header = T)
dfsnap_with_unassigned <- read.delim(src.file, header = TRUE)
dfsnap <-
  dfsnap_with_unassigned[dfsnap_with_unassigned$taxid != "_unassigned_", ]

stopifnot(length(dfsnap[, 1]) == length(dfres[, 1])) # are the lengths the same?
stopifnot(data.frame(lapply(dfsnap[order(dfsnap$taxid), ][['taxid']], as.character)) ==
            data.frame(lapply(dfres[order(dfres$taxid), ][['taxid']], as.character))) # are the taxids the same?


df <- merge(dfsnap, dfres, by = "taxid")

sel <- df[with(df, order(-h1_norm)),][1:kSelected, ]

first.phylogroup <- sel$phylogroup[[1]]
first.serotype <- sel$Serotype.From.Reads[[1]]
first.cat <- sel$pnc[[2]]
first.res <- CatToCatName(first.cat)

second.phylogroup <- unique(sel$phylogroup)[[2]]

phylogroups.masked =
  3 * as.integer(sel$phylogroup > -1) -
  2 * as.integer(sel$phylogroup == first.phylogroup) -
  1 * as.integer(sel$phylogroup == second.phylogroup)

labs <- paste(sel$taxid, sel$h1)  ## create labels

phylogroups <- max(sel$phylogroup)
vals.to.plot <- sel$h1_norm
maxval <- max(vals.to.plot)
maxtick <- (ceiling(100 * maxval)) / 100.0
#print(maxtick)

par(mar = c(2.5, 1.5 , 1.5 , 0.0))

antibiotics <- DfToAnts(df)
nants = length(antibiotics)
AbsResGridHeight = kResGridHeight * maxval

x <- barplot(
  height = vals.to.plot,
  #col = CatToColor(sel$pnc_cat),
  col = phylogroups.masked,
  border = NA,
  space = 0,
  axes = F,
  ylim = c(-(nants + kGridShift) * AbsResGridHeight, maxtick),
  xlim = c(-2, kSelected),
  beside = T,
  mgp = c(4.8, 0, 0)
)

i = 0
for (ant in rev(antibiotics)) {
  clm <- paste(ant, "_cat", sep = "")
  y <- AbsResGridHeight * (nants - i + kGridShift)
  xx <- barplot(
    height = 0 * vals.to.plot - y,
    col = CatToColor(sel[[clm]]),
    border = NA,
    space = 0.0,
    axes = F,
    add = T
  )

  text(
    x = -0.3,
    y = -y + 0.5 * AbsResGridHeight,
    col = "#000000",
    labels = toupper(ant),
    pos = 2,
    offset = 0.15,
    srt = 00,
    cex = 0.90
  )

  i = i + 1
}

y <- kGridShift * AbsResGridHeight
xx <- barplot(
  height = 0 * vals.to.plot - y,
  col = "#ffffff",
  border = NA,
  space = 0.0,
  axes = F,
  add = T
)

#ablines(c(0,0.01), col="white")

# grid(lwd = 1,
#      ny = maxval / 0.01,
#      nx = NA)



axis(
  side = 2,
  pos = -0.2,
  cex.axis = 0.8,
  at = seq(0, maxtick, 0.01),
  mgp = c(50, 0.5, 0)
)


line1 <- paste0("Test sample: ", sample.desc)
line2 <- paste0(
  "Results after 5 mins: Penicillin ",
  first.res,
  ", lineage ",
  first.phylogroup,
  ", serotype ",
  first.serotype
)


mtext(
  "Relative similarity to sample",
  side = 2,
  line = -0.3,
  cex.lab = 2,
  las = 3,
  cex = 1.6
)


# mtext(
#   "(arbitrary units)",
#   side = 2,
#   line = -0.0,
#   cex.lab = 1,
#   cex = 1.15,
#   las = 3
# )


## taxid
text(
  #y = tmp,
  y = 0,
  x = x,
  labels = sel$taxid,
  pos = 1,
  cex = 0.05,
  offset = 0.35,
  col = "#eeeeee"
  #srt = 40
)

#tmp <- sel$h1_norm
# # serotype
# text(
#   y = tmp,
#   x = x + 0.6,
#   labels = sel$Serotype.From.Reads,
#   pos = 3,
#   cex = 0.75,
#   offset = 0.35,
#   col = "#000000",
#   srt = 40
# )



# # phylogroup
# text(
#   x = xx,
#   y = 0,
#   col = "#ffffff",
#   labels = as.numeric(sel$phylogroup),
#   pos = 1,
#   offset = -0.35,
#   srt = 00,
#   cex = 0.45
# )

abline(v = seq(0:kSelected),
       col = "white")

text(
  y = -kGridShift * 0.5 * AbsResGridHeight,
  x = 35,
  labels = "Isolate",
  cex = 1.75,
  col = "#000000"
)


bps.total <- sum(dfsnap_with_unassigned[c("ln")])
bps.voting <- sum(dfsnap_with_unassigned[c("h1")])
reads <- sum(dfsnap_with_unassigned[c("count")])
subtitle <- paste(
  "Reads: ",
  format(as.integer(reads), big.mark = ","),
  "       Bps: ",
  format(as.integer(bps.total), big.mark = ","),
  "       Useful bps: ",
  round(100 * as.double(bps.voting) / as.double(bps.total)),
  "%",
  sep = ""
)


# bps. etc.
mtext(
  side = 1,
  text = subtitle,
  line = 0.8,
  cex = 0.95
)


# # legend - resistance
#par(lend = 1)
legend(
  x = "topright",
  title = "Susceptibility",
  legend = c(CatToCatName("R"), CatToCatName("S")),
  cex = 1.5,
  fill = kCatColors,
  y.intersp = 0.8,
  #lty = 1,
  #lwd = 10,
  yjust = 3,
  border = NA,
  box.col = NA,
  inset = c(0, 0.25),
  title.adj = 0
)


# legend - seq. phylogroups

legend(
  x = "topright",
  title = "Phylogroup",
  legend = c(first.phylogroup, second.phylogroup, "Others"),
  cex = 1.5,
  fill = c(1, 2, 3),
  y.intersp = 0.8,
  box.col = NA,
  border = NA,
  inset = c(0.142, 0.00),
  title.adj = 0
)




# timestamp
mtext(
  side = 1,
  text = paste(as.POSIXct(timestamp, origin = "1970-01-01"), "  "),
  line = 0.8,
  adj = 1,
  cex = 0.95
)


if (!kRStudio) {
  dev.off()
}
