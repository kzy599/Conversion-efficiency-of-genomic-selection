for(r in c(1:10)){
rm(list=setdiff(ls(),"r"))
gc()
#prepare for simulation
#load R package
source("package.R")

#globe parameters
source("parameters.R")

basewd = paste(ptwd,"BasePoP",sep = "")

setwd(basewd)

founderPop = quickHaplo(nInd=1000, nChr=10, segSites=2000, inbred=FALSE)

tvar = c(1,1,1,1,1,1)
tmean = c(0,0,0,0,0,0)
tcor = diag(ncol=6,nrow=6)
tcor[1,2:6] = c(0.1,0.3,0.5,0.7,0.9)
tcor[2:6,1] = c(0.1,0.3,0.5,0.7,0.9)
th2 = rep(0.1,6)
#set genetic parameters
SP = SimParam$new(founderPop)
SP$restrSegSites(minSnpFreq = 0.05,overlap = TRUE)
SP$addTraitA(100,mean=tmean,var = tvar,corA =tcor)
SP$setVarE(h2=th2)#遗传力和偏差
SP$addSnpChip(114,minSnpFreq = 0.05)#55kSNP芯片
SP$setSexes("yes_sys")#按照一个雌一个雄来分配个体的性别
pop_founder = newPop(founderPop, simParam=SP)

pop <- selectCross(pop_founder,
                   nFemale = nDam,nMale = nSire,
                   nCrosses = nCrosses,nProgeny = nProgenyPerCross,
                   use = "rand",
                   simParam = SP)

rm(.Random.seed)

save.image(file = paste(basewd,"/","basePOP",r,".rda",sep=""))

}