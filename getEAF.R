#get eaf from 1000G
if(is.na(file2$eaf.outcome)){
  dat<-file2
  path<-"/Users/sunyuqi/Desktop/fileFrequency.frq"
  get_eaf_from_1000G(dat, path, type = "exposure")
}