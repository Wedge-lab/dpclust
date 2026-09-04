#!/usr/bin/env Rscript
# DPClust pipeline entry point
# All argument parsing, validation, logging, and execution are handled
# by DPClust::dpclust_cli(). See ?dpclust_cli for available options.
suppressMessages(library(DPClust))
DPClust::dpclust_cli()
