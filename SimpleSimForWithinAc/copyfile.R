#!~/bin/Rscript
library(data.table)

ldr = list.dirs(path = getwd(),full.names = TRUE)

ldr = setdiff(ldr,getwd())

filedir = getwd()

allfile = list.files(filedir)

allfile = allfile[allfile%flike%".R"]

for(i in 1:length(ldr)){
file.copy(from = paste(filedir,allfile,sep="/"),to = ldr[i])
}

