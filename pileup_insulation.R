#!/usr/bin/Rscript
suppressPackageStartupMessages({
library(tidyverse, quietly = T)
library(regioneR, quietly = T)
library(plyranges, quietly = T)
library(optparse, quietly = T)
})
# command line options for multiple alignment pipeline
option_list = list(
  make_option(c("-i", "--insulation"),
              type="character",
              default=NULL,
              help="score.bedgraph file output from hicFindTads",
              metavar = "character"),
  make_option(c("-f", "--features"),
              type="character",
              default = NULL,
              help="Feature file (3 column BED) to aggregate insulation score over",
              metavar = "character"),
  make_option(c("-l", "--feature_label"),
              type="character",
              default = "features",
              help="Feature file to aggregate insulation score over",
              metavar = "character"),
  make_option(c("-s", "--score_species"),
              type="character",
              default = "hg38",
              help="Label for insulation score species (hg38)",
              metavar = "character"),
  make_option(c("-o", "--outfile"),
              type="character",
              default = "output.pdf",
              help = "Output file (.pdf)",
              metavar = "character"))

# parse command line arguments
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# require args
if (is.null(opt$insulation) | is.null(opt$features) | is.null(opt$outfile)) {
  print_help(opt_parser)
  stop("Error, must supply all command line arguments.", call. = FALSE)
}

# reead breaks of syn into df
read_feature <- function(bos) {
  df <- read_delim(bos, delim="\t",
                   col_names = c("seqnames", "start", "end"))
  return(as_granges(df))
}

# read insulation
read_score <- function(ins) {
  df <- read_delim(ins, delim='\t', col_names = c("seqnames", "start", "end", "score"))
  return(as_granges(df))
}

features <- read_feature(opt$features)
ins <- read_score(opt$insulation)

# aggregate insulation score over extended feature ranges
features %>% plyranges::mutate(start=start-3e5, end=end+3e5) %>%
  plyranges::filter(start>0, end>0) %>%
  tile_ranges(width=width(.)/500) %>%
  plyranges::group_by(partition) %>%
  plyranges::mutate(bin=1:n()) %>%
  ungroup() %>%
  join_overlap_left(ins) %>%
  plyranges::group_by(bin) %>%
  plyranges::summarise(medianScore=median(score, na.rm = T)) %>%
  as_tibble() -> bos_bin_summ

# ensure 500 rows
bos_bin_summ <- head(bos_bin_summ, 500)

# plot features
p <- ggplot(bos_bin_summ, aes(x=bin, y=medianScore)) +
  geom_point() +
  geom_smooth() +
  labs(title=paste("Median Insulation Score Near",opt$feature_label,"Features"),
       y="median score",
       x=paste(opt$feature_label, "position (Kb)")) +
  scale_x_continuous(breaks=c(0, 125, 250, 375, 500), labels=c("-300","-150","0","150","300")) +
  theme_classic()

ggsave(plot=p, opt$outfile, width = 8, height = 5)
