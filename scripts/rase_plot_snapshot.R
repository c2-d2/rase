#!/usr/bin/env Rscript

#' Plot RASE prediction snapshot.
#'
#' Author: Karel Brinda <kbrinda@hsph.harvard.edu>
#'
#' License: MIT
#'

suppressMessages(suppressWarnings(library(optparse)))

# CONFIGURATION -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

set.seed(42)

kWidth <- 8
kHeight <- 7
kSelected <- 70  # strains to display
kResGridHeight <-
    0.045  # heigth of a resistance cell, fraction of the height of the graph above
kGridShift <- 2
kPalette <-
    rep(c("#d95f02", "#7570b3", "#222222", "#e7298a"), times = 2)
kRStudio <- Sys.getenv("RSTUDIO") == "1"
kCatColors <- c("#ff0000", "#0000aa", "#00aa00")  # R / S / NA


# FUNCTIONS -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

CatToColor <- function(cat) {
    catu <- toupper(cat)
    ifelse(catu == "R", kCatColors[1], (ifelse(catu == "S", kCatColors[2], kCatColors[3])))
}

CatToNumber <- function(cat) {
    catu <- toupper(cat)
    ifelse(catu == "R", 1, (ifelse(catu == "S", 2, 3)))
}

CatToCatName <- function(cat) {
    catu <- toupper(cat)
    ifelse(catu == "S", "susceptible", (ifelse(
        catu == "R", "non-susceptible", "unknown"
    )))
}

DfToAnts <- function(dataframe) {
    cols <- colnames(dataframe)
    antcols <- cols[grepl("_cat", cols)]
    ants <- gsub("_cat", "", antcols)
    ants <- ants[ants != "CHL"]
    ants <- ants[ants != "chl"]
    ants
}

LastLine <- function(ct, bps.total, kmers.matched) {
    paste(
        "Reads: ",
        format(as.integer(ct), big.mark = ","),
        "       Bps: ",
        format(as.integer(bps.total), big.mark = ","),
        "       Matched k-mers: ",
        round(100 * as.double(kmers.matched) / as.double(bps.total)),
        "%",
        sep = ""
    )
}


# CLI-parsing -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if (kRStudio) {
    src.file <- "tests/snapshot.tsv"
    met.file <- "tests/metadata.tsv"
    timestamp <- 1505967676
    sample.desc <- "Sample desc."
} else {
    option.list <-
        list(make_option(c("-d", "--desc"), default = F, help = "sample description"))
    parser <-
        OptionParser(usage = "%prog [options] metadata.tsv snapshot.tsv plot.pdf", option_list = option.list)

    arguments <- parse_args(parser, positional_arguments = 3)
    opt <- arguments$options

    sample.desc <- ""
    if (opt$desc) {
        sample.desc <- opt$desc
    }

    met.file <- arguments$args[1]
    src.file <- arguments$args[2]
    out.file <- arguments$args[3]

    timestamp <-
        as.integer(rev(strsplit(src.file, "[\\./]")[[1]])[2])

    pdf(out.file, width = kWidth, height = kHeight)

}


# PLOTTING-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

palette(kPalette)

df.met <- read.delim(met.file, header = T, stringsAsFactors = F)
# support for legacy formats
colnames(df.met)[colnames(df.met) == "phylogroup"] <- "lineage"
colnames(df.met)[colnames(df.met) == "pg"] <- "lineage"

df.snap.all <- read.delim(src.file, header = T, stringsAsFactors = F)
# support for legacy formats
colnames(df.snap.all)[colnames(df.snap.all) == "phylogroup"] <-
    "lineage"
colnames(df.snap.all)[colnames(df.snap.all) == "pg"] <-
    "lineage"

df.snap <- df.snap.all[df.snap.all$taxid != "_unassigned_",]

stopifnot(length(df.snap[, 1]) == length(df.met[, 1]))  # are the lengths the same?
taxids.snap = data.frame(lapply(df.snap[order(df.snap$taxid),][["taxid"]], as.character))
taxids.met = data.frame(lapply(df.met[order(df.met$taxid),][["taxid"]], as.character))
stopifnot(taxids.snap == taxids.met)  # are the taxids the same?

df <-
    merge(df.snap,
          df.met,
          by = c("taxid", "lineage"))

sel <- df[with(df, order(-weight)),][1:kSelected,]

first.lineage <- sel$lineage[[1]]
first.cat <- sel$pnc[[2]]
first.resistant <- CatToCatName(first.cat)

second.lineage <- unique(sel$lineage)[[2]]

lineages.masked <-
    3 * as.integer(sel$lineage > -1) - 2 * as.integer(sel$lineage == first.lineage) - 1 * as.integer(sel$lineage == second.lineage)

labs <- paste(sel$taxid, sel$weight)  ## create labels

vals.to.plot <- sel$weight_norm
maxval <- max(vals.to.plot)
maxtick <- (ceiling(100 * maxval)) / 100

par(mar = c(2.5, 1.5, 1.5, 0))

antibiotics <- DfToAnts(df)
nants <- length(antibiotics)
AbsResGridHeight <- kResGridHeight * maxval

x <-
    barplot(
        height = vals.to.plot,
        col = lineages.masked,
        border = NA,
        space = 0,
        axes = F,
        ylim = c(-(nants + kGridShift) *
                     AbsResGridHeight, maxtick),
        xlim = c(-2, kSelected),
        beside = T,
        mgp = c(4.8, 0, 0)
    )

i <- 0
for (ant in rev(antibiotics)) {
    clm <- paste(ant, "_cat", sep = "")
    y <- AbsResGridHeight * (nants - i + kGridShift)
    xx <-
        barplot(
            height = 0 * vals.to.plot - y,
            col = CatToColor(sel[[clm]]),
            border = NA,
            space = 0,
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
        srt = 0,
        cex = 0.9
    )

    i <- i + 1
}

y <- kGridShift * AbsResGridHeight
xx <-
    barplot(
        height = 0 * vals.to.plot - y,
        col = "#ffffff",
        border = NA,
        space = 0,
        axes = F,
        add = T
    )

axis(
    side = 2,
    pos = -0.2,
    cex.axis = 0.8,
    at = seq(0, maxtick, 0.01),
    mgp = c(50, 0.5, 0)
)

mtext(
    "Weight (normalized)",
    side = 2,
    line = -0.3,
    cex.lab = 2,
    las = 3,
    cex = 1.6
)



## taxid
text(
    y = 0,
    x = x,
    labels = sel$taxid,
    pos = 1,
    cex = 0.05,
    offset = 0.35,
    col = "#eeeeee"
)


abline(v = seq(0:kSelected), col = "white")

text(
    y = -kGridShift * 0.5 * AbsResGridHeight,
    x = 35,
    labels = "Isolate",
    cex = 1.75,
    col = "#000000"
)

bps.total <- sum(df.snap.all[c("length")])
kmers.matched <- sum(df.snap.all[c("weight")])
reads <- sum(df.snap.all[c("count")])
subtitle <- LastLine(reads, bps.total, kmers.matched)

# bps. etc.
mtext(
    side = 1,
    text = subtitle,
    line = 0.8,
    cex = 0.95
)


# legend - resistance, susc.
legend(
    x = "topright",
    title = "Susceptibility",
    legend = c(CatToCatName("R"), CatToCatName("S")),
    cex = 1.5,
    fill = kCatColors,
    y.intersp = 0.8,
    yjust = 3,
    border = NA,
    box.col = NA,
    inset = c(0, 0.25),
    title.adj = 0
)

# legend - lineages
legend(
    x = "topright",
    title = "Lineage",
    legend = c(first.lineage, second.lineage, "Others"),
    cex = 1.5,
    fill = c(1, 2, 3),
    y.intersp = 0.8,
    box.col = NA,
    border = NA,
    inset = c(0.142, 0),
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
