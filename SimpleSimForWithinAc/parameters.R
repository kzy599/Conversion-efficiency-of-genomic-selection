ptwd = "/home/kangziyi/DRcand/checkwac/"
curDirname = substr(getwd(),nchar(ptwd)+1,nchar(getwd()))
usePed = FALSE
Fonly = FALSE
useRand = TRUE
nGeneration = 20
TD = 25 # 25 45 65 90
nCrosses = 50
nSire = 25
nDam = 50
nSelection = 100
nProgenyPerCross = 300
# 6 candidates within each family, total 300
nHDsib = 1
nChr = 10
ChrSize = (2.6 * 10^9) / 44
MutRate = 2.5E-7#
RecRate = 1.67E-8#，RecRate=GenLen/ChrSize

nT = c(1,0.1,0.3,0.5,0.7,0.9)


if(curDirname!="BasePoP"){

 ncand = 6*nCrosses

 nSR = str_extract_all(curDirname, "\\d+") %>% unlist()
 
 if(length(nSR)==1){nTra = 1}else{
  nTra = paste(nSR[1],nSR[2],sep = ".") %>% as.numeric()

  nTra = which(nT==nTra)
 }


}
