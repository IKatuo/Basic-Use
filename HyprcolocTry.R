#需要先找locus，不能全基因组coloc##
#将aging outcomes按snp进行merge；各Trait的beta，se读取为matrix

library(dplyr)
library(data.table)
library(hyprcoloc)
#获取file list
trait_folder_path<-"/Users/sunyuqi/Desktop/生命科学竞赛/原始数据/longevity_outcome/重命名out/telomere_length"
trait_file_path<-list.files(trait_folder_path,full.names = TRUE)

# 初始化一个空列表来存储数据框
data_list <- list()
dat_hypr<-list()
i=1
# 循环读取每个txt文件
for (file in trait_file_path) {
  data <- fread(file, fill=T, header = TRUE) 
  dat_hypr$SNP<-data$SNP
  dat_hypr$beta<-data$beta.outcome
  dat_hypr$se<-data$se.outcome
  a<-paste0("beta",i)
  b<-paste0("se",i)
  names(dat_hypr)<-c("SNP",a,b)
  data_list[[i]] <- dat_hypr
  i=i+1
  rm(data)
  dat_hypr<-list()
}

# 初始化交集数据框
snp_intersect <- data_list[[1]]

# 循环取交集并提取beta数据，假设有6个trait##
m<-2
for (m in 2:6) {
  # 取交集
  snp_intersect <- merge(snp_intersect, data_list[[m]], by = "SNP")
  m<-m+1
}

# 选择snp和beta列，并转换为矩阵
beta_matrix <- as.matrix(snp_intersect[c("SNP", paste0("beta", 1:6))])
ses<-as.matrix(snp_intersect[c("SNP", paste0("se", 1:6))])
rsid <- snp_intersect$SNP
traits <- paste0("T", 1:6)
# Colocalization analysis  
hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)
#接下来可以画图，如果有需要的话#

