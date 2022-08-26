library(rmarkdown)
library(utils)
files <- list.files("../", pattern=".Rmd", recursive=TRUE, full.names=TRUE)
for (i in 1:length(files)){
    render(files[i], output_format='pdf_document')
    print(i)
}
