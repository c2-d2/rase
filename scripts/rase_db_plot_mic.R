#!/usr/bin/env Rscript --vanilla

#
# Author:  Karel Brinda <kbrinda@hsph.harvard.edu>
#
# License: MIT
#

suppressMessages(suppressWarnings(library(optparse)))


# CONFIGURATION -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

set.seed(42)

kPch <- 3

kWidth <- 8
kHeight <- 7

kHeadRstudio <- 150
#kHeadRstudio <- 200

kRStudio <- Sys.getenv("RSTUDIO") == "1"

# R / S / r / s / NA
kCatColors <-
    c("#ff0000", "#0000aa", "#ff9999", "#9999ff", '#00aa00')


# FUNCTIONS-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

CatToColor <- function(cat) {
    catu <- toupper(cat)
    ifelse(cat == catu,
        ifelse(catu == "R", kCatColors[1], (
            ifelse(catu == "S", kCatColors[2], kCatColors[5])
        )),
        ifelse(catu == "R", kCatColors[3], (ifelse(
            catu == "S", kCatColors[4], 5
        ))))
}

CatToNumber <- function(cat) {
    catu <- toupper(cat)
    ifelse(cat == catu,
        ifelse(catu == "R", 1, (ifelse(catu == "S", 2, 5))),
        ifelse(catu == "R", 3, (ifelse(catu == "s", 4, 5))))
}

CatToCatName <- function(cat) {
    catu <- toupper(cat)
    ifelse(catu == "S", "susceptible", (ifelse(
        catu ==
            "R", "non-susceptible", "unknown"
    )))
}

DfToAnts <- function(dataframe) {
    cols <- colnames(dataframe)
    antcols <- cols[grepl("_cat", cols)]
    ants <- gsub("_cat", "", antcols)
    ants
}

InfToBig <- function(vec) {
    ifelse(vec == Inf,
        9999,
        ifelse(vec == 0, 0.00000000000000000001, vec))
}

GetPhylogroupStarts <- function(df) {
    starts <- c(1)
    for (i in seq(nrow(df) - 1)) {
        if (df[i,]$phylogroup != df[i + 1,]$phylogroup) {
            starts <- c(starts, i + 1)
        }
    }
    starts <- c(starts, nrow(df) + 1)
    starts
}

GetSubPhylogroupStarts <- function(df) {
    if ("subsubphylogroup" %in% colnames(df)) {
        starts <- c(1)
        for (i in seq(nrow(df) - 1)) {
            if (df[i,]$subphylogroup != df[i + 1,]$subphylogroup) {
                starts <- c(starts, i + 1)
            }
        }
        starts <- c(starts, nrow(df) + 1)
        starts
    } else{
        c()
    }
}


# CLI-parsing -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

zoomed <- T

if (kRStudio) {
    res.file <- "../spneumoniae_abc/resistance/res_cat.tsv"
    kHead <- kHeadRstudio
    breakpoint <- 2
    ant <- "pen"
    df.len <- nrow(read.delim(res.file, header = T))
} else {
    option.list <- list(
        make_option(
            c("-b", "--breakpoint"),
            default = NA,
            dest = 'breakpoint',
            type = 'double',
            help = "breakpoint to plot"
        ),
        make_option(
            c("-z", "--zoom"),
            action = "store_true",
            default = F,
            dest = 'zoomed',
            help = "zoomed version"
        )
    )
    parser <-
        OptionParser(usage = "%prog [options] res_cats.tsv antibiotic plot.pdf", option_list = option.list)

    arguments <- parse_args(parser, positional_arguments = 3)
    opt <- arguments$options

    # sample.desc <- ""
    # if (opt$desc) {
    #   sample.desc <- opt$desc
    # }
    breakpoint = opt$breakpoint
    zoomed = opt$zoomed
    #print(opt)

    res.file <- arguments$args[1]
    ant <- arguments$args[2]
    out.file <- arguments$args[3]

    df.len <- nrow(read.delim(res.file, header = T))

    if (zoomed) {
        pdf(out.file,
            width = df.len / 10,
            height = kHeight)
    } else {
        pdf(out.file,
            width = kWidth,
            height = kHeight)

    }

    kHead <- df.len
}

print(paste("File", res.file))
print(paste("Antibiotic", ant))
print(paste("Breakpoint", breakpoint))


# PLOTTING -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

df_big <- read.delim(res.file, header = T)

par(mar = c(4.4, 4.2, 1.5 , 1.0))

if ("order" %in% colnames(df_big)) {
    the.order <- order(df_big$order)
} else if ("subsubphylogroup" %in% colnames(df_big)) {
    the.order <- order(df_big$phylogroup,
        df_big$subphylogroup,
        df_big$subsubphylogroup)
} else if ("subsubphylogroup" %in% colnames(df_big)) {
    the.order <- order(df_big$phylogroup, df_big$subphylogroup)
} else {
    the.order <- order(df_big$phylogroup)
}

df <- head(df_big[the.order, ], kHead)

ant_int_col <- paste(ant, "_int", sep = "")
ant_cat_col <- paste(ant, "_cat", sep = "")
ant_intervals <- df[[ant_int_col]]
l_int <- as.numeric(gsub("(.*)-(.*)", "\\1", ant_intervals))
r_int <- as.numeric(gsub("(.*)-(.*)", "\\2", ant_intervals))
df$l <- InfToBig(l_int)
df$r <- InfToBig(r_int)
l_int_df <- data.frame(l_int)
r_int_df <- data.frame(r_int)
rmax <- max(r_int_df[r_int_df != Inf, ])
rmin <- min(r_int_df[r_int_df != 0, ])

k.ylim <- c(log2(rmin) - 1, log2(rmax) + 3)

plot(
    1,
    type = "n",
    xaxs = "i",
    xlab = "Isolate",
    ylab = bquote(log[2] ~ .(paste(
        toupper(ant), "MIC", sep = "_"
    ))),
    xlim = c(1, nrow(df)),
    ylim = k.ylim
)

abline(h = log2(breakpoint))

segments(seq(nrow(df)),
    log2(df$l),
    seq(nrow(df)),
    log2(df$r),
    col = CatToColor(df[, ant_cat_col, ]),
    lwd = 1)

points(log2(df$l), col = CatToColor(df[, ant_cat_col, ]), pch = kPch)
points(log2(df$r), col = CatToColor(df[, ant_cat_col, ]), pch = kPch)

pg_starts <- GetPhylogroupStarts(df)
spg_starts <- GetSubPhylogroupStarts(df)
if (zoomed) {
    for (i in seq(lengths(spg_starts))) {
        middle <- floor((spg_starts[i] + spg_starts[i + 1]) / 2.0)
        text(
            x = middle ,
            y = log2(rmin) - 0.6,
            cex = 0.6,
            label = paste("SPG", df[middle,]$subphylogroup)
        )
    }
}

for (i in seq(lengths(pg_starts))) {
    middle <- floor((pg_starts[i] + pg_starts[i + 1]) / 2.0)
    pg = df[middle,]$phylogroup
    if (zoomed) {
        pg.label = paste("PG", pg)
    } else{
        pg.label = pg
    }
    text(x = middle ,
        y = log2(rmin) - 1,
        label = pg.label)
}

if (zoomed) {
    abline(
        v = spg_starts - 0.5,
        col = "black",
        lty = 2,
        lwd = 1
    )

    abline(
        v = pg_starts - 0.5,
        col = "black",
        lty = 1,
        lwd = 2
    )
} else {
    abline(
        v = pg_starts - 0.5,
        col = "black",
        lty = 2,
        lwd = 2
    )
}

print(breakpoint)

legend(
    "topright",
    c(
        "non-susc.",
        "susceptible",
        "inferred non-susc.",
        "inferred susc."
    ),
    col = kCatColors,
    pch = kPch,
    lwd = NA,
    bg = "white"
)

if (zoomed) {
    taxidPos = pmin(pmax(head(log2(df$r), kHead) + 0.1, 0 * head(df$r, kHead) +
            k.ylim[1] + 1), k.ylim[2] - 1)
    text(
        labels = head(df$taxid, kHead),
        x = seq(kHead),
        pos = 4,
        offset = 0.12,
        y = taxidPos,
        cex = 0.6,
        srt = 90
    )
}

if (!kRStudio) {
    dev.off()
}
