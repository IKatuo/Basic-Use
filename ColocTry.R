###使用说明：除了文件路径，和$后面的内容，其他什么都不要动就可以直接运行##
##不要更改群内对应的文件夹结构和文件名称##
library(coloc)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(readxl)
library(writexl)
library(Get_MR)
###aging data formatted###
##get input for coloc##
# 读取outcome和exposure合并后的FDR corrected results表格
data <- read_xlsx("/Users/sunyuqi/Desktop/生命科学竞赛/原始数据/step2+3/step2/FDRGathered.xlsx")

# 设置两个文件夹的路径，一个用于outcome文件，一个用于exposure文件
out_folder_path <- "/Users/sunyuqi/Desktop/生命科学竞赛/原始数据/longevity_outcome/重命名out"
exp_folder_path <- "/Users/sunyuqi/Desktop/生命科学竞赛/原始数据/step2+3/step2/cg"
m<-0
# 遍历data中的每一行
for (i in 1:nrow(data)) {
  # 获取当前行的outcome和exposure值
  outcome <- data$outcome[i]
  exposure <- data$exposure[i]
  m<-m+1
  print(m)
  # 构建.txt文件的完整路径
  out_file_path <- file.path(out_folder_path, paste0(outcome, ".txt"))
  
  # 如果.txt文件存在，则读取它
  if (file.exists(out_file_path)) {
    file2 <- fread(out_file_path) # 假设.txt文件是制表符分隔的
    # 在这里处理file2，例如与data中的其他行合并等
  } else {
    warning(paste("No outcome file found for outcome:", outcome))
  }
  #get eaf from 1000G
  if(is.na(file2$eaf.outcome)){
    dat<-file2
    path<-"/Users/sunyuqi/Desktop/fileFrequency.frq"
    get_eaf_from_1000G(dat, path, type = "outcome")
  }
  # 构建表格文件所在文件夹的完整路径
  exp_subfolder_path <- file.path(exp_folder_path, exposure)
  
  # 如果文件夹存在，则列出其中的所有文件
  if (dir.exists(exp_subfolder_path)) {
    # 获取文件夹中所有表格文件的完整路径
    exp_file_paths <- list.files(path = exp_subfolder_path, full.names = TRUE)
    
    # 遍历每个表格文件
    for (file in exp_file_paths) {
      # 读取表格文件
      file1 <- read_xlsx(file)
      # 在这里处理file1，例如与data中的其他行或file2合并等
      #直接获取input，注意by=“SNP”
      input <- merge(file1, file2, by="SNP", all=FALSE)
      head(input)
      input_path<-paste0("/Users/sunyuqi/Desktop/生命科学竞赛/原始数据/step2+3/step2/coloc_input/",outcome,"_",exposure,basename(file))
      # 将input保存为新的txt文件
      write_xlsx(input,input_path) # 使用制表符作为分隔符，根据你的数据调整
      rm(file1)
    }
  } else {
    warning(paste("No folder found for exposure:", exposure))
  }
  rm(file2)
}

# 在这里可以对file1和file2进行进一步处理，比如合并等
# 例如：combined_data <- merge(file1, file2, by = "common_column")



##coloc function for quantitative trait## 
coloc.analysis.quant <- function(snp,beta1,beta2,se1,se2,MAF1,MAF2,N1,N2){
  
  #Convert the inputs in order to run in coloc function.
  #type, quant (quantitative) for pQTL study 
  dataset1 <- list(snp=snp,beta=beta1, varbeta=se1^2, MAF=MAF1,type="quant", N=N1)
  
  
  #type, quant (quantitative) for quantitative trait study 
  dataset2 <- list(snp=snp,beta=beta2, varbeta=se2^2,MAF=MAF2, type="quant",N=N2)
  
  
  #Run the coloc analysis, setting the prior probabilities for association with each trait (p1, p2) and both traits together (p12) as 1E-5.
  #p1 prior probability a SNP is associated with trait 1, default 1e-4
  #p2 prior probability a SNP is associated with trait 2, default 1e-4
  #p12 prior probability a SNP is associated with both traits, default 1e-5
  result <- coloc.abf(dataset1, dataset2, p1=1e-4, p2=1e-4, p12=1e-5)  
  
  #Format the data to save out.
  
  #List into data frame.
  df <- data.frame(matrix(unlist(result$summary), nrow=1, byrow=T))
  df <- df[complete.cases(df), ] 
  #Label the columns in the data frame.
  names(df) <- c("nsnps", "PP.H0.abf",    "PP.H1.abf",    "PP.H2.abf",    "PP.H3.abf",    "PP.H4.abf")
  
  #Make the filename and save out.
  return(df)
  
}

#prepare MAF reference from 1000G EUR#
MAF_ref<-fread("/Users/sunyuqi/Desktop/fileFrequency.frq")
#call the input file
input_file_list<-list.files("/Users/sunyuqi/Desktop/生命科学竞赛/原始数据/step2+3/step2/coloc_input（1800-1877）",full.names = TRUE)
#pair count
v=0
for (file in input_file_list){
df <- read_xlsx(file)
if(nrow(df)==0){
  print(paste0("input is NULL",basename(file)))
  }else{
#get eaf from 1000G
v=v+1
print(v)

#check input complete
n=0
df$MAF1<-NULL
df$MAF2<-NULL
if(!is.na(df$eaf.outcome[1])){
  #get maf
  for(j in 1:nrow(df)){
    if(df$eaf.exposure[j]>0.5){
      df$MAF1[j]=1-df$eaf.exposure[j]
      }else{df$MAF1[j]=df$eaf.exposure[j]}
    if(df$eaf.outcome[j]>0.5){
      df$MAF2[j]=1-df$eaf.outcome[j]
    }else{df$MAF2[j]=df$eaf.outcome[j]}
    }
  }else{
  for(j in 1:nrow(df)){
    if(df$eaf.exposure[j]>0.5){
      df$MAF1[j]=1-df$eaf.exposure[j]
      }else{df$MAF1[j]=df$eaf.exposure[j]
      }
  }
    indices <- match(df$SNP, MAF_ref$SNP)
      # 如果找到匹配项，则从MAF_ref中提取MAF值
      df$MAF2 <- MAF_ref$MAF[indices]
}
n<n+1
print(n)

df <- df[!duplicated(df$SNP),]
df <- df[which(df$MAF1>0),]
df <- df[which(df$MAF2>0),]
#sample size for two studies
N1<-df$samplesize.exposure
N2<-df$samplesize.outcome

df<-as.data.frame(df)
result <- coloc.analysis.quant(df$SNP,df$beta.exposure,df$beta.outcome,df$se.exposure,df$se.outcome,df$MAF1,df$MAF2,N1,N2)
write_xlsx(result,paste0("/Users/sunyuqi/Desktop/生命科学竞赛/原始数据/step2+3/step2/coloc_result/",basename(file)))
}
}
