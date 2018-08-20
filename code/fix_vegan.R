remove.packages(c(“phyloseq”, “vegan”)
library(devtools)
install_version("vegan", version = "2.4-6", repos = "http://cran.us.r-project.org”)
library(biocInstaller)
biocLite("phyloseq")

