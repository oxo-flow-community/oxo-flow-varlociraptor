#!/usr/bin/env Rscript
# Annotate COSMIC signature exposures per VAF threshold.
#
# Port of the upstream workflow/scripts/annotate_mutational_signatures.R
# (MIT, snakemake-workflows dna-seq-varlociraptor v6.10.0), adapted to argv:
#   annotate_mutational_signatures.R COSMIC.tsv CONTEXT.tsv
#     --group GROUP --vafs 5,10,...,100 --outputs OUT1.tsv,OUT2.tsv,... --log LOG
#
# The snakemake script API is replaced by positional inputs + flags (the
# output list and vaf thresholds are comma-joined). The siglasso fitting
# loop with the running prior update is preserved verbatim; stdout and
# messages are sunk to the log file, matching the upstream script's capture
# of the R session output. Thresholds are given in percent (divided by 100
# before filtering, as upstream).

args <- commandArgs(trailingOnly = TRUE)

get_flag <- function(flag) {
    idx <- match(flag, args)
    if (is.na(idx)) {
        stop(paste("missing required argument", flag))
    }
    args[idx + 1]
}

input_cosmic <- args[1]
input_context <- args[2]
group <- get_flag("--group")
min_vafs <- as.integer(strsplit(get_flag("--vafs"), ",")[[1]])
outputs <- strsplit(get_flag("--outputs"), ",")[[1]]
log_file <- get_flag("--log")

sink(log_file, type = "output")
sink(log_file, type = "message", append = TRUE)

library(siglasso)
library(tibble)
library(dplyr)
library(readr)
library(purrr)
library(stringr)

# Load COSMIC signatures
cosmic_signatures <- read_tsv(input_cosmic)
cosmic_signatures <-  as.matrix(cosmic_signatures
                                %>% mutate(Type = gsub("\\[|\\]", "", Type))
                                %>% column_to_rownames(var = "Type")
                    )

# Add a Prefix to sample names for correct handling of numerical group names
sample_substitutions <- (
    read_tsv(input_context)
    %>% mutate(Sample = paste0("X_", Sample))
    )
if (nrow(sample_substitutions) == 0) {
    for (output_file in outputs) {
        write_tsv(tibble(), output_file)
    }
} else {
    # Replace dots by - in data.frame header caused by siglasso
    replace_dots <- function(df) {
        names(df) <- str_replace_all(names(df), "\\.", "-")
        return(df)
    }

    prior <- rep(1, ncol(cosmic_signatures))
    for (i in seq_along(outputs)) {
        output_file <- outputs[[i]]
        min_vaf <- min_vafs[[i]] / 100
        print(min_vaf)
        filtered_substitions <- (
            sample_substitutions
            %>% filter(AF >= min_vaf)
            %>% mutate(AF = NULL)
        )
        # Skip entries with equal or less than one substitution
        if (nrow(filtered_substitions) <= 1) {
            write_tsv(tibble(), output_file)
        } else {
            spectrum <- context2spec(filtered_substitions, plot=FALSE)
            sample_signatures <- (
                as.data.frame(siglasso(spectrum, cosmic_signatures, prior=prior, plot=FALSE))
                %>% rownames_to_column(var="Signature")
                %>% replace_dots()
                %>% filter(!!sym(paste0("X_", group)) > 0)
                %>% add_column(Frequency = min_vaf)
            )
            write_tsv(sample_signatures, output_file, col_names=FALSE)
            prior <- colnames(cosmic_signatures) %>% map_dbl(~ if_else(any(sample_signatures$Signature == .x), 0.1, 1))
        }
    }
}
