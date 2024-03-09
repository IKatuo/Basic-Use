library(dplyr)
# 设置文件夹路径
parent_folder <- "/Users/sunyuqi/Desktop/生命科学竞赛/原始数据/step2+3/step2/res"
gathered_folder <- "/Users/sunyuqi/Desktop/生命科学竞赛/原始数据/step2+3/step2/Gathered"
# 获取子文件夹列表
subfolders <- list.dirs(parent_folder, recursive = FALSE)

# 循环处理每个子文件夹
for (subfolder in subfolders) {
  # 获取子文件夹名称
  subfolder_name <- basename(subfolder)
  
  # 获取子文件夹中所有表格文件
  files <- list.files(subfolder, full.names = TRUE)
  
  # 创建一个空的结果数据框
  result <- data.frame()
  
  # 循环读取每个表格文件
  for (file in files) {
    # 读取表格文件
    data <- read_xlsx(file) 
    if (nrow(data) > 0) {
      file_name<-basename(file)
      data<-data %>% arrange(pval)
      data$FDR<-(nrow(data)/seq_len(nrow(data)))*data$pval
      filtered_data<-data %>% filter(FDR<=0.05)
      save_path<-file.path("/Users/sunyuqi/Desktop/生命科学竞赛/原始数据/step2+3/step3/FDR",paste0(file_name))
      if (nrow(filtered_data) > 0) {
        write_xlsx(filtered_data, save_path)}
    }
}}

