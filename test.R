path_code <- paste0(dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))]))), "/")
print(path_code)
