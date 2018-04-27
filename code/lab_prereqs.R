# Install prerequisites for Micro360 lab
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("dada2")
biocLite("phyloseq")
biocLite("Rqc")
install.packages(c("knitr", "devtools", "kableExtra", "ggplot2", "dplyr", "vegan", "Rmisc", "caTools"))

