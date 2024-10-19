#!~/bin/Rscript
#globe script
#for(r in c(1:10)){
#rm(list=setdiff(ls(),"r"))
#gc()

#prepare for simulation
#load R package
source("package.R")

tempname = substr(ngenpop,1,nlen)

if(tempname == "allpop"){
 nGeneration  = 20 #burn-in generation
}else if(tempname == "allpop3"){
    nGeneration = 3 #burn-in generation
}

nFemale = 50

nMale = 25

nCrosses = 50

nCand = 300 

nProgenyPerCross = 200
#creat base and G0 population
source("population.R")
#}