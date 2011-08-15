#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
parser <- OptionParser(usage = "%prog infname outfname", add_help_option=TRUE)

parser_out <- parse_args(parser, positional_arguments = TRUE)
args <- parser_out$args
opts <- parser_out$options

if (length(args) != 2) {
    print_help(parser)
    quit(save="no", status=1)
}

infname <- args[1]
outfname <- args[2]

if (!file.exists(infname))
    stop("Couldn't find input file ", infname)
if (file.exists(outfname) && infname != outfname)
    stop("Output file ", outfname, " already exists")

suppressPackageStartupMessages(library("niftir"))
img <- read.nifti.image(infname)
hdr <- read.nifti.header(infname)

if (length(hdr$dim)!=5)
    stop("Image has wrong dimensions")

img <- img[,,,1,]
hdr$dim <- hdr$dim[-4]
hdr$pixdim <- hdr$pixdim[-4]

if (infname == outfname)
    invisible(file.remove(infname))

invisible(write.nifti(img, hdr, outfile=outfname))
