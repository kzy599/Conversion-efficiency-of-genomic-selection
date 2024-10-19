for(g in 0:nGeneration){
gc()
#===========collect data====================
#=========allocate the training and valid set================
pop_cand_F = selectWithinFam(pop = pop,nInd = nSelection/2,use = "rand",sex = "F")
pop_cand_M = selectWithinFam(pop = pop,nInd = nSelection/2,use = "rand",sex = "M")
pop_cand = c(pop_cand_F,pop_cand_M)
 
if(usePed){

  if(!exists("candidate")){
    pop_par = pop_founder[unique(c(pop@mother,pop@father))]
    ped_ne = data.table(id = pop_par@id, sir = NA, dam = NA)
    colnames(ped_ne) = c("id","sir","dam")
  }else{
    pop_par = candidate[unique(c(pop@mother,pop@father))]
    ped_temp = data.table(id = pop_par@id, sir = NA, dam = NA)
    colnames(ped_temp) = c("id","sir","dam")
    ped_ne = rbind(ped_ne,ped_temp)
  }

}
needfam= unique(pop_cand@mother)


cand_list = lapply(needfam,function(x){

        pop_v = pop_cand[pop_cand@mother %in% x ]
        
        if(Fonly){
          tempgv = rep(mean(pop_v@gv[,1]),pop_v@nInd)
        }else{
          tempgv = pop_v@gv[,nTra] - mean(pop_v@gv[,nTra]) + mean(pop_v@gv[,1])
          }

        pop_v@ebv = matrix(tempgv,ncol = 1)
        
        return(pop_v)
    })




#=============================================================================================

#numberofid= lapply(cand_list,function(x){
 #   x= x@nInd
  #  return(x)
#})

#numberofid = unlist(numberofid)

#if(sum(is.na(numberofid))>0) save.image(file = paste("errorPOP",r,".rda",sep=""))


if(!Fonly){
acforeach= lapply(cand_list,function(x){
    x= cor(x@ebv,x@gv[,1])

    return(x)
})

acforeach = unlist(acforeach)
output$minacforeach[g+1] = min(acforeach)
output$maxacforeach[g+1] = max(acforeach)
output$acforeach[g+1] = mean(acforeach)

output$acall[g+1] = acforall(cand_list = cand_list)
}else{

output$minacforeach[g+1] = 0
output$maxacforeach[g+1] = 0
output$acforeach[g+1] = 0

output$acall[g+1] = 0
}





#===================================================================================================


##===================================allcocate the number of candidate within each candidate family=================================================

cand_list = lapply(cand_list,function(x){
    needcand = ncand/nCrosses
    
    if(Fonly){
    x_F = selectWithinFam(x,nInd = round((needcand)*2/3),sex = "F",use = "rand" )
    
    x_M = selectWithinFam(x,nInd = round((needcand)*1/3),sex = "M",use = "rand" )
    }else{
    x_F = selectWithinFam(x,nInd = round((needcand)*2/3),sex = "F",use = "ebv" )
    
    x_M = selectWithinFam(x,nInd = round((needcand)*1/3),sex = "M",use = "ebv" )
    }

    x = c(x_F,x_M)
    return(x)
})

candidate = mergePops(cand_list)


##==========================================================================================

#===============================calculate the population parameters======================

source("calparameters.R")

#========================================================================================



#======================generate the NP of the next generation============================================
if(g<nGeneration){
  
if(usePed){
  pop <- ocs(
    pop = candidate,
    usePed = usePed,
    ped_parents = ped_ne,
    nCrosses = nCrosses,
    nProgenyPerCross = nProgenyPerCross,
    nFemalesMax = nDam,
    nMalesMax = nSire,
    equalizeFemaleContributions =TRUE,
    equalizeMaleContributions = TRUE,
    targetDegree = TD,
    use = "ebv"
  )
}else if(useRand){
  pop <- selectCross(candidate,
                   nFemale = nDam,nMale = nSire,
                   nCrosses = nCrosses,nProgeny = nProgenyPerCross,
                   use = "rand",
                   simParam = SP)
}else{
pop <- ocs(
    pop = candidate,
    nCrosses = nCrosses,
    nProgenyPerCross = nProgenyPerCross,
    nFemalesMax = nDam,
    nMalesMax = nSire,
    equalizeFemaleContributions =TRUE,
    equalizeMaleContributions = TRUE,
    targetDegree = TD,
    use = "ebv"
  )
}
   

}

}
warnings()


if(TD==45){
  prec = ""
}else if(TD==25){
prec = "l"
}else if(TD==65){
  prec = "m"
}else if(TD == 90){
prec = "h"
}

if(useRand) prec = "n"

if(usePed){

if(Fonly){
fwrite(output,paste(prec,"pop",0,r,".csv",sep = ""),sep = ",")
}else{
 fwrite(output,paste(prec,"spop",nTra,r,".csv",sep = ""),sep = ",")

}

}else{

if(Fonly){
  fwrite(output,paste(prec,"spop",0,r,".csv",sep = ""),sep = ",")
}else{
  fwrite(output,paste(prec,"pop",nTra,r,".csv",sep = ""),sep = ",")
}


}
