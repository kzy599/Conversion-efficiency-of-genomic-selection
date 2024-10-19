for(r in 1:10){
    rm(list = setdiff(ls(),"r"))
    gc()
    
    source("package.R")

    
    source("parameters.R")
    
    basewd = paste(ptwd,"BasePoP",sep = "")
    
    load(paste(basewd,"/","basePOP",r,".rda",sep=""))

    #founderPop = quickHaplo(nInd=1000, nChr=10, segSites=2500, inbred=FALSE)
    #SP = SimParam$new(founderPop)
    #SP$addTraitA(mean = 0,var = 1,nQtlPerChr = 100)
    #SP$setVarE(h2=0.15)
    #SP$addSnpChip(1250)
    #SP$setSexes("yes_sys")

    #pop_founder = newPop(founderPop, simParam=SP)
    
    #source("makelowpenal.R")


    #pop <- selectCross(pop_founder,
     #                  nFemale = nDam,nMale = nSire,
      #                 nCrosses = nDam,nProgeny = nProgenyPerCross,
       #                use = "rand",)
    
    source("package.R")

    
    source("parameters.R")



output <- data.frame(generation=0:(nGeneration),
                     mean_gv=numeric(nGeneration+1),
                     inbreeding = numeric(nGeneration+1),
                     inbreeding_plink = numeric(nGeneration+1),
                     LD=numeric(nGeneration+1),
                     LDscore=numeric(nGeneration+1),
                     Ne = numeric(nGeneration+1),
                     Va = numeric(nGeneration+1),
                     Vp = numeric(nGeneration+1),
                     Vg = numeric(nGeneration+1),
                     mean_pheno=numeric(nGeneration+1),           
                     genicVa = numeric(nGeneration+1),
                     genicVg = numeric(nGeneration+1),
                     h2 = numeric(nGeneration+1),
                     nCoancestor = numeric(nGeneration+1),
                     acimp = numeric(nGeneration+1),
                     accuracy = numeric(nGeneration+1))

    
    source("perform.R")

}