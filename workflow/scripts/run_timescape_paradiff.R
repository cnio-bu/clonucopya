#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
library("timescape")
library("htmlwidgets")

# create parser object
parser <- ArgumentParser()

parser$add_argument("--prev", action="store")
parser$add_argument("--mutations", action="store")
parser$add_argument("--locations", action="store")
parser$add_argument("--tree", action="store")
parser$add_argument("--out_file", action="store")

args <- parser$parse_args()

#results_clonal_prev <- read.csv(file = '/home/cleon/PycharmProjects/MSc-Project/out/timescape/results_clonal_prev.csv')
results_clonal_prev_ts <- read.csv(file = args$prev)

#results_mutations <- read.csv(file = '/home/cleon/PycharmProjects/MSc-Project/out/timescape/results_mutations.csv')
results_mutations_ts <- read.csv(file = args$mutations)

#results_tree <- read.csv(file = '/home/cleon/PycharmProjects/MSc-Project/out/timescape/results_tree.csv')
results_tree_ts <- read.csv(file = args$tree)

#TODO: make it for multi-sample
results_perturbations <- data.frame(pert_name=c("16T320"), prev_tp=c("18T59-I"))

res_results_ts = timescape(clonal_prev = results_clonal_prev_ts, tree_edges = results_tree_ts, perturbations = results_perturbations, mutations = results_mutations_ts)

res_results_ts = timescape(clonal_prev = results_clonal_prev_ts,
                       tree_edges = results_tree_ts,
                       perturbations = results_perturbations,
                       mutations = results_mutations_ts)

options(browser="/usr/bin/google-chrome")

##browserURL(res_results)

res_results_ts

# Export
#results_img_ref <- "/home/cleon/PycharmProjects/MSc-Project/out/timescape/933124/res_results_ts.html"
results_ts_output_html <- args$out_file

saveWidget(res_results_ts,
            results_ts_output_html,
            selfcontained = FALSE,
            libdir = NULL,
            background = "white",
            title = class(res_results_ts)[[1]],
            knitrOptions = list())
