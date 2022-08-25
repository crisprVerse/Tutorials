library(rmarkdown)
library(utils)
files <- list.files(pattern=".Rmd", recursive=TRUE, full.names=TRUE)
render(files[1], output_format='pdf_document')