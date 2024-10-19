ptwd = "/home/mengxh/kangzy/WssvCode/WssvCodeFinal/"
curDirname = substr(getwd(),nchar(ptwd)+1,nchar(getwd()))

nGeneration = 20

nCrosses = 50
nSire = 25
nDam = 50
nSelection = 100
nProgenyPerCross = 300
# 6 candidates within each family, total 300
ncand = 300
nHDsib = 1
nChr = 10
ChrSize = (2.6 * 10^9) / 44
MutRate = 2.5E-7#
RecRate = 1.67E-8#，RecRate=GenLen/ChrSize

lowdensity = c(3, 12, 23, 46, 114, 228, 455, 682) * nChr

if(curDirname!="BasePoP"){

#Its important to set a correct parameters in cratedirandpbs.sh
if(curDirname%flike%"pop"){gsmode = "Pop"}else{gsmode = "Fam"}

if(curDirname%flike%"IT"){imp = TRUE}else{imp = FALSE}

nSR = str_extract_all(curDirname, "\\d+")

nSR = unlist(nSR)

lowpenal = as.character(as.numeric(nSR[1])*nChr)

nReference = as.numeric(nSR[2])
#

if(imp == TRUE){
    Im = "IT"
}else{
    Im = "IF"
}

}
