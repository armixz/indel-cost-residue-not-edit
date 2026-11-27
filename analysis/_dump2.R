args <- commandArgs(trailingOnly=TRUE)
x <- readRDS(args[1])
cat("class:", paste(class(x),collapse=","), " length:", length(x), "\n")
str(x, max.level=2, list.len=20)
