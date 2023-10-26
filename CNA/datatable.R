#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V="Version: 1.1"
D="Depends: R (>= 3.4.0), optparse, data.table"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-e", "--dt_fun"), type="character", default="NA",
     help="data.table function (required). Possible value is 'dcast(...)' [%default].", metavar="character"),
    # fread
     make_option(c("--out_quote"), type="character", default="FALSE",
     help="Should output be surrounded by double quotes, 'auto'|TRUE|FALSE [%default]", metavar="character"),
    make_option(c("--str2fac"), type="logical", default="FALSE", action = "store_true",
     help="Should strings be converted to factors [%default].", metavar="logical"),
    make_option(c("--col_sep"), type="character", default="\t",
     help="Separator between columns [\\t] \n\t\t'auto' represents '[,\\t |;:]'. See also 'data.table::fread'.", metavar="character"),
    make_option(c("--row_num"), type="integer", default="-1",
     help="Number of rows to read [%default] \n\t\t'-1' means all. '0' is a special case that just returns the column names.", metavar="integer"),
     make_option(c("--col_name"), type="character", default="TRUE",
     help="The first data line contain column names, 'auto'|TRUE|FALSE [%default] \n\t\tDefaults according to whether every non-empty field on the first data\n\t\tline is type character.", metavar="character"),
    make_option(c("--na_str"), type="character", default="NA",
     help="String(s) which are to be interpreted as NA values [%default]", metavar="character"),
    make_option(c("--row_skip"), type="character", default="0",
     help="Row can be taken as the first data row [%default] \n\t\tskip='string' searches for 'string' in the file and starts on that row.", metavar="character"),
    make_option(c("--row_name"), type="character", default="1",
     help="Column number (or name) can be taken as the row names [%default]", metavar="character")
)

opt_parser <- OptionParser(usage="usage: %prog [options] <input.tsv> \n\t<input.tsv> is a TSV file (separator: '\\t', stdin: '-').", option_list=option_list, description = paste(V, D, sep="\n"))
opt <- parse_args(opt_parser, args=a, positional_arguments=TRUE)

e <- opt$options$dt_fun

oq <- opt$options$out_quote
sf <- opt$options$str2fac
cs <- opt$options$col_sep
rn <- opt$options$row_num
cn <- opt$options$col_name
na <- opt$options$na_str
rs <- opt$options$row_skip
ra <- opt$options$row_name

if(oq != 'auto') oq=as.logical(oq)
if(cn != 'auto') cn=as.logical(cn)
if(grepl("^[[:digit:]]+$", rs)) rs=as.numeric(rs)
if(grepl("^[[:digit:]]+$", ra)) ra=as.numeric(ra)

rv0 <- as.numeric(R.version$major)
rv1 <- as.numeric(R.version$minor)

if ((rv0<3)|((rv0==3)&(rv1<5))) {
    d <- ifelse(opt$args[1]=='-', 'file:///dev/stdin', opt$args[1])
} else {
    d <- ifelse(opt$args[1]=='-', 'cat /dev/stdin', opt$args[1])
}
# read
d <- fread(d, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, stringsAsFactors=sf, data.table=TRUE) # D_Table

# main
od <- eval(parse(text=e))

# output
fwrite(od, quote=oq, sep=cs, na=na, col.names=cn)

### THE END ###
