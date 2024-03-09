with open('/Users/sunyuqi/Desktop/大创/DNAm/HbA1c/GCST90014006_buildGRCh38.tsv','r')as f:
    lines= f.readlines()
    head_lines = [line for line in lines if line.split()[0] == '#chromosome']
    filtered_lines=[line for line in lines if line.split()[0]== '6'and (line.split()[1]>161582957 and line.split()[1]<162582957)]
with open('/Users/sunyuqi/Desktop/大创/DNAm/HbA1c/GLP1R_500kb_chr6selected.txt','w')as file:
    file.writelines(head_lines)
    file.writelines(filtered_lines)





