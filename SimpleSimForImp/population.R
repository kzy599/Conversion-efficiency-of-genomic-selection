#creadt base poopulation and G0 population
#Historical effective population size 
#histNe=c(305,447,509,588,676,737,842,965,1115,1317,1554,1887,2195,2651,3266,3837,4495,4978,5619,6427,7129,7551,7568,8517,8253)
histNe=c(1115,1317,1554,1887,2195,2651,3266,3837,4495,4978,5619,6427,7129,7551,7568,8517,8253)
#histGen=c(13,15,17,20,23,27,32,38,45,54,65,80,98,122,151,188,238,294,357,454,555,708,784,877,952)
histGen=c(45,54,65,80,98,122,151,188,238,294,357,454,555,708,784,877,952)
BaseNe=1000
MaCSeNFlags = ""
if(length(histNe)>0){
  histNe = histNe/BaseNe
  histGen = histGen/(4*BaseNe)
  for(i in 1:length(histNe)){
    MaCSeNFlags = paste(MaCSeNFlags,"-eN",
                        histGen[i],histNe[i])
  }
}
ChrSize = (2.6 * 10^9) / 44
MutRate = 2.5E-7#突变率
RecRate = 1.67E-8#重组率，RecRate=GenLen(遗传距离)/ChrSize
founderPop = runMacs(nInd = 1000,
                     nChr = 44,
                     segSites = 2700,
                     manualCommand = paste(as.integer(ChrSize),
                                           "-t", MutRate * 4 * BaseNe,
                                           "-r", RecRate * 4 * BaseNe,
                                           MaCSeNFlags),
                     manualGenLen = RecRate * ChrSize)

#easy population
#founderPop = quickHaplo(nInd=1000, nChr=10, segSites=2000, inbred=FALSE)
#set genetic parameters
SP = SimParam$new(founderPop)
SP$restrSegSites(minSnpFreq = 0.05,overlap = TRUE)
#SP$addTraitA(100,mean=c(0,0),var = c(5.39,5.39),corA = matrix(c(1,0.8,0.8,1),nrow = 2))
#SP$addTraitA(100,mean=c(0,0),var = c(5.39,5.39),corA = matrix(c(1,0.5,0.5,1),nrow = 2))
SP$addTraitA(100,mean=c(0,0),var = c(1,1),corA = matrix(c(1,0,0,1),nrow = 2))
SP$setVarE(h2=c(0.41,0.10))#
SP$addSnpChip(80,minSnpFreq = 0.05)#
SP$setSexes("yes_sys")#
pop_founder = newPop(founderPop, simParam=SP)
pop <- selectCross(pop_founder,
                   nFemale = nFemale,nMale = nMale,
                   nCrosses = nCrosses,nProgeny = nProgenyPerCross,
                   use = "rand",
                   simParam = SP)

for(g in c(0:nGeneration)){

pop_cand_F = selectWithinFam(pop = pop,nInd = nProgenyPerCross/4,use = "rand",sex = "F")
pop_cand_M = selectWithinFam(pop = pop,nInd = nProgenyPerCross/4,use = "rand",sex = "M")
pop_cand = c(pop_cand_F,pop_cand_M)
pop_test = pop[setdiff(pop@id,pop_cand@id)]


#==========collect parents data==========
if(!exists("candidate")){
 pop_parents = pop_founder[unique(c(pop@mother,pop@father))]

 dt_parents = data.table(id = pop_parents@id, sir = NA, dam = NA, 
           fam = NA,nGeneration = "founder",suvtime = NA,BW = NA)

 colnames(dt_parents) = c("id","sir","dam","fam","nGeneration","suvtime","BW")
 
 ped_ne = dt_parents[,1:3]
}else{
 pop_parents = candidate[unique(c(pop@mother,pop@father))]

 dt_parents1 = data.table(id = pop_parents@id, sir = pop_parents@father, dam = pop_parents@mother, 
           fam = paste(pop_parents@father,pop_parents@mother,sep = ""),nGeneration = (g-1),suvtime = NA,BW = NA)

 colnames(dt_parents1) = c("id","sir","dam","fam","nGeneration","suvtime","BW")

 dt_parents = rbind(dt_parents,dt_parents1)

 ped_ne = dt_parents[,1:3]
}

#==================================================================


#==============================collect cand and test data==============================
if(!exists("alldt")){
alldt_test = data.table(id = pop_test@id, sir = pop_test@father, dam = pop_test@mother, 
           fam = paste(pop_test@father,pop_test@mother,sep = ""),nGeneration = g,suvtime = pop_test@pheno[,2],BW = pop_test@pheno[,1])

colnames(alldt_test) = c("id","sir","dam","fam","nGeneration","suvtime","BW")
}else {
alldt_test1 = data.table(id = pop_test@id, sir = pop_test@father, dam = pop_test@mother, 
           fam = paste(pop_test@father,pop_test@mother,sep = ""),nGeneration = g,suvtime = pop_test@pheno[,2],BW = pop_test@pheno[,1])

colnames(alldt_test1) = c("id","sir","dam","fam","nGeneration","suvtime","BW")

alldt_test = rbind(alldt_test,alldt_test1)

if(g>=4){
   alldt_test = alldt_test[nGeneration%in%c((g-4):g),]
}

}
dt_cand = data.table(id = pop_cand@id, sir = pop_cand@father, dam = pop_cand@mother, 
           fam = paste(pop_cand@father,pop_cand@mother,sep = ""),nGeneration = g,suvtime = NA,BW = NA)
colnames(dt_cand) = c("id","sir","dam","fam","nGeneration","suvtime","BW")

alldt = rbind(alldt_test,dt_cand)
#========================================================================================================================
##collect extra reference data
if(exists("pop_gaint")){
alldt_extra = data.table(id = pop_gaint@id, sir = pop_gaint@father, dam = pop_gaint@mother, 
           fam = paste(pop_gaint@father,pop_gaint@mother,sep = ""),nGeneration = g,suvtime = pop_gaint@pheno[,2],BW = pop_gaint@pheno[,1])

colnames(alldt_extra) = c("id","sir","dam","fam","nGeneration","suvtime","BW")
}

#

#==============================collect cand and test data==============================

if(!exists("genodt_t")){
  genodt_t = pop_test
}else{
  genodt_t = c(genodt_t,pop_test)
}

if(g>=4){
   genodt_t= genodt_t[genodt_t@id %in% alldt$id]
}

genodt = c(genodt_t,pop_cand)

genodt = pullSnpGeno(genodt)
#============================================================================================================


#===========================create next generation====================================
if(g<nGeneration){

  needcand = nCand/nCrosses
  
  pop_cand@ebv = matrix((pop_cand@pheno[,1] + pop_cand@pheno[,2])/2, ncol = 1)

  
  candidate_F = selectWithinFam(pop_cand,nInd = round((needcand)*2/3),sex = "F",use = "ebv" )
  
  candidate_M = selectWithinFam(pop_cand,nInd = round((needcand)*1/3),sex = "M",use = "ebv" )

  candidate = c(candidate_F,candidate_M)


  if(g == (nGeneration - 1)){


    pop_gaint <- selectCross(candidate,
                   nFemale = nFemale,nMale = nMale,
                   nCrosses = nCrosses,nProgeny = (5000+nProgenyPerCross),
                   use = "rand",
                   simParam = SP)
    
    pop_temp_F = selectWithinFam(pop_gaint,nInd = (nProgenyPerCross/2),sex = "F",use = "rand")

    pop_temp_M = selectWithinFam(pop_gaint,nInd = (nProgenyPerCross/2),sex = "M",use = "rand")

    pop = c(pop_temp_F,pop_temp_M)

    pop_gaint = pop_gaint[!pop_gaint@id %in% pop@id]

  }else{

    pop <- selectCross(candidate,
                   nFemale = nFemale,nMale = nMale,
                   nCrosses = nCrosses,nProgeny = nProgenyPerCross,
                   use = "rand",
                   simParam = SP)

  }
}

#===============================================================

}

#===================================four generations data for pop============================
mapdt = makemap()

 genodtput = makeped(z = genodt)

 fwrite(genodtput,file="hib.ped",col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

 fwrite(mapdt, file = "hib.map",col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

 fwrite(alldt,file = "pheno.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

 system("plink --file hib --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib")
 
 system("hiblup --make-xrm --threads 32 --bfile hib --add --out hib")
 
 system(paste("hiblup --mme --pheno pheno.csv --pheno-pos 6 --xrm hib.GA --vc-priors ",varA(pop)[2,2],",",(varP(pop)[2,2]-varA(pop)[2,2]),
               " --pcg --threads 32 --out hib",sep = ""))

 hiebv_ST<- fread("hib.rand",sep = "\t")

 system(paste("hiblup --mme --pheno pheno.csv --pheno-pos 7 --xrm hib.GA --vc-priors ",varA(pop)[1,1],",",(varP(pop)[1,1]-varA(pop)[1,1]),
               " --pcg --threads 32 --out hib",sep = ""))

 hiebv_BW<- fread("hib.rand",sep = "\t")

#=============================================================================



#=============save the workplace===========================
#save.image(file = paste("allpop",r,".rda",sep=""))

save.image(file = ngenpop )
