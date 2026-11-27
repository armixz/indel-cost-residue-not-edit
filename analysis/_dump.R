args <- commandArgs(trailingOnly=TRUE)
x <- readRDS(args[1])
cat("class:", class(x), "\n")
if (is.list(x) && !is.data.frame(x)) {
  cat("list names:", paste(names(x), collapse=", "), "\n")
  for (n in names(x)) {
    e <- x[[n]]
    cat("---", n, "class", class(e)[1], "dim", paste(dim(e), collapse="x"), "\n")
    if (is.data.frame(e)) {
      cat("   cols:", paste(names(e), collapse=","), "\n")
      write.csv(e, file.path(args[2], paste0(n, ".csv")), row.names=FALSE)
    }
  }
} else if (is.data.frame(x)) {
  cat("cols:", paste(names(x), collapse=","), "\n")
  write.csv(x, file.path(args[2], "table.csv"), row.names=FALSE)
}
