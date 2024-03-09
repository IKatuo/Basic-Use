library(data.table)
library(TwoSampleMR)
library(stringr)
setwd("C://Users/95476/Desktop/mQTL/step1")
file <- list.files(path = "./tissue_gene")
file <- as.list(file)
for (i in c(1:length(file))) {
  name <- file[i]
  name <- as.character(name)
  new <- readxl::read_xlsx(paste0("tissue_gene/",name))
  new$eaf <- NULL
  for (j in c(1:nrow(new))){
         if (new$ref_factor[j]==1) {new$eaf[j] <- new$maf[j]} else {new$eaf[j] <- 1-new$maf[j]}
  }
  writexl::write_xlsx(new,path = paste0("new/",name))
}
