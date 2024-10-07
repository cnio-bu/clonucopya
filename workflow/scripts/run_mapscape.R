#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
library("mapscape")
library("htmlwidgets")

# create parser object
parser <- ArgumentParser()

parser$add_argument("--prev", action="store")
parser$add_argument("--mutations", action="store")
parser$add_argument("--locations", action="store")
parser$add_argument("--tree", action="store")
parser$add_argument("--image", action="store")
parser$add_argument("--out_file", action="store")

args <- parser$parse_args()

#results_clonal_prev <- read.csv(file = '/home/cleon/PycharmProjects/MSc-Project/out/mapscape/results_clonal_prev.csv')
results_clonal_prev <- read.csv(file = args$prev)

#results_mutations <- read.csv(file = '/home/cleon/PycharmProjects/MSc-Project/out/mapscape/results_mutations.csv')
results_mutations <- read.csv(file = args$mutations)

#results_sample_locations <- read.csv(file = '/home/cleon/PycharmProjects/MSc-Project/out/mapscape/results_sample_locations.csv')
results_sample_locations <- read.csv(file = args$locations)

#results_tree <- read.csv(file = '/home/cleon/PycharmProjects/MSc-Project/out/mapscape/results_tree.csv')
results_tree <- read.csv(file = args$tree)

#results_img_ref <- "/home/cleon/PycharmProjects/MSc-Project/out/mapscape/A21_anatomical_image.png"
results_img_ref <- args$image

#TODO: make it for multi-sample Paradiff
#results_sample_ids <- c("16T320","18T59-I")
#TODO: make it for multi-sample AML
results_sample_ids <- c("tumor","relapse")

#res_results = mapscape(clonal_prev = results_clonal_prev, tree_edges = results_tree, sample_locations = results_sample_locations, mutations = results_mutations, img_ref = results_img_ref, sample_ids = results_sample_ids)
res_results = mapscape(clonal_prev = results_clonal_prev,
                       tree_edges = results_tree,
                       sample_locations = results_sample_locations,
                       mutations = results_mutations,
                       show_low_prev_gtypes = FALSE,
                       img_ref = results_img_ref,
                       sample_ids = results_sample_ids,
                       width = 960, height = 960)

options(browser="/usr/bin/google-chrome")

##browserURL(res_results)

res_results

# Export
#results_img_ref <- "/home/cleon/PycharmProjects/MSc-Project/out/mapscape/res_results.html"
results_output_html <- args$out_file

saveWidget(res_results,
            results_output_html,
            selfcontained = FALSE,
            libdir = NULL,
            background = "white",
            title = class(res_results)[[1]],
            knitrOptions = list())
