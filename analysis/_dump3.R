args <- commandArgs(trailingOnly=TRUE)
x <- readRDS(args[1])
d <- x[[1]]
write.csv(as.data.frame(d), args[2], row.names=FALSE)
cat("rows", nrow(d), "\n")
print(table(d$mut_type))
print(table(d$domain))
