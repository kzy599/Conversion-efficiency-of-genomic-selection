for(g in 0:nGeneration){
gc()
#===========collect data====================
#=========allocate the training and valid set================
pop_cand_F = selectWithinFam(pop = pop,nInd = nSelection/2,use = "rand",sex = "F")
pop_cand_M = selectWithinFam(pop = pop,nInd = nSelection/2,use = "rand",sex = "M")
pop_cand = c(pop_cand_F,pop_cand_M)


pop_test = pop[setdiff(pop@id,pop_cand@id)]

if((nReference*nCrosses)!=pop_test@nInd){
pop_test_F = selectWithinFam(pop = pop_test,nInd = nReference/2,use = "rand",sex = "F")
pop_test_M = selectWithinFam(pop = pop_test,nInd = nReference/2,use = "rand",sex = "M")
pop_test = c(pop_test_F,pop_test_M)
}
#==========================================

#==========collect pedigree data for calpopparameters==========

if(!exists("candidate")){
 pop_parents = pop_founder[unique(c(pop@mother,pop@father))]

 dt_parents = data.table(id = pop_parents@id, sir = NA, dam = NA, 
           fam = NA,nGeneration = "founder",suvtime = NA)

 colnames(dt_parents) = c("id","sir","dam","fam","nGeneration","suvtime")
 
 ped_ne = dt_parents[,1:3]
}else{
 pop_parents = candidate[unique(c(pop@mother,pop@father))]

 dt_parents1 = data.table(id = pop_parents@id, sir = pop_parents@father, dam = pop_parents@mother, 
           fam = paste(pop_parents@father,pop_parents@mother,sep = ""),nGeneration = (g-1),suvtime = NA)

 dt_parents = rbind(dt_parents,dt_parents1)

 ped_ne = dt_parents[,1:3]
}

#==================================================================



#pedigree and phenotype

#if(gsmode == "Pop"){
if(!exists("alldt_test")){
alldt_test = data.table(id = pop_test@id, sir = pop_test@father, dam = pop_test@mother, 
           fam = paste(pop_test@father,pop_test@mother,sep = ""),nGeneration = g,suvtime = pop_test@pheno)

colnames(alldt_test) = c("id","sir","dam","fam","nGeneration","suvtime")
}else {
alldt_test1 = data.table(id = pop_test@id, sir = pop_test@father, dam = pop_test@mother, 
           fam = paste(pop_test@father,pop_test@mother,sep = ""),nGeneration = g,suvtime = pop_test@pheno)

colnames(alldt_test1) = c("id","sir","dam","fam","nGeneration","suvtime")

alldt_test = rbind(alldt_test,alldt_test1)

if(g>=3){
   alldt_test = alldt_test[nGeneration%in%c((g-3):g),]
}

}
dt_cand = data.table(id = pop_cand@id, sir = pop_cand@father, dam = pop_cand@mother, 
           fam = paste(pop_cand@father,pop_cand@mother,sep = ""),nGeneration = g,suvtime = NA)
colnames(dt_cand) = c("id","sir","dam","fam","nGeneration","suvtime")
#}

#if(gsmode == "Fam"){
 #   alldt_test = data.table(id = pop_test@id, sir = pop_test@father, dam = pop_test@mother, 
  #         fam = paste(pop_test@father,pop_test@mother,sep = ""),nGeneration = g,suvtime = pop_test@pheno)
           
  #  colnames(alldt_test) = c("id","sir","dam","fam","nGeneration","suvtime")

   # dt_cand = data.table(id = pop_cand@id, sir = pop_cand@father, dam = pop_cand@mother, 
    #       fam = paste(pop_cand@father,pop_cand@mother,sep = ""),nGeneration = g,suvtime = NA)
    #colnames(dt_cand) = c("id","sir","dam","fam","nGeneration","suvtime")
#}

alldt = rbind(alldt_test,dt_cand)


#============genotype=============================
if(gsmode == "Pop"){

if(!exists("genodt_t")){
     genodt_t = pop_test
    }else{
     genodt_t = c(genodt_t,pop_test)
    }

if(g>=3){
   genodt_t= genodt_t[genodt_t@id %in% alldt$id]
}
genodt = c(genodt_t,pop_cand)


}else{

 genodt = pop

}

#if(gsmode == "Pop" | accumulate){

#if(!exists("genodt_t")){
 #    genodt_t = pop_test
  #  }else{
   #  genodt_t = c(genodt_t,pop_test)
   # }

#if(g>=4){
 #  genodt_t= genodt_t[genodt_t@id %in% alldt$id]
#}
#genodt = c(genodt_t,pop_cand)


#}else{

 #genodt = pop

#}

if(lowpenal %in% as.character(lowdensity)){
#============impute========================================================================
if(imp == TRUE){
#correct the numberof snp based on the penal we chosed
#genodt[,!colnames(genodt)%in%snpPenal[[lowpenal]]$id] = 9
if(!exists("candidate")){
pop_p = pop_parents

}else{
pop_p = c(pop_p,pop_parents)
}

if(nHDsib!=0){

if(!exists("HDsib")){

    HDsib = selectWithinFam(pop = pop_test,nInd = nHDsib, use = "rand")

}else{

    HDsib = c(HDsib,selectWithinFam(pop = pop_test,nInd = nHDsib, use = "rand"))
    
}

}

if(g>=3){
    HDsib = HDsib[HDsib@id %in% genodt@id]
}
#While the information from the previous 4 generations was initially deemed sufficient for 
#genotype imputation in the current dataset, 
#we opted to utilize data from the previous 8 generations. 
#This choice was made because the population mode contains data from the most recent 4 generations, 
#which necessitates access to information from earlier generations. 
#To ensure consistency in generation length between both modes of imputation, 
#we employed information from the past 8 generations for both modes
if(g>6){
   needparents= dt_parents[nGeneration%in%c((g-7):(g-1)),]
   pop_p = pop_p[pop_p@id%in%needparents$id]
}else{
   needparents = dt_parents
}

#仅用单世代填充，然后累代填充后的数据叠加，减少运算时间
#pop_p = pop_parents
#genodt = pop
impgeno = c(pop_p,genodt)

impgeno = pullSnpGeno(impgeno)

if(exists("HDsib")){


    impgeno[rownames(impgeno)%in%setdiff(genodt@id,HDsib@id),!colnames(impgeno)%in%snpPenal[[lowpenal]]$id] = 9

}else{

    impgeno[rownames(impgeno)%in%genodt@id,!colnames(impgeno)%in%snpPenal[[lowpenal]]$id] = 9
}

lociname = colnames(impgeno)

if(gsmode == "Fam"){
 impdt= rbind(needparents,alldt[nGeneration == g,])
}else{
 impdt = rbind(needparents,alldt)
}

#if(gsmode == "Pop"| accumulate){
 #impdt = rbind(needparents,alldt)
#}else{
 #impdt= rbind(needparents,alldt[nGeneration == g,])
#}

#如果仅用单世代填充
#needparents = data.table(id = pop_parents@id,sir = 0, dam = 0)
#colnames(needparents) = c("id","sir","dam")
#currentdt = data.table(id = pop@id,sir = pop@father,dam = pop@mother)
#colnames(currentdt) = c("id","sir","dam")
#impdt = rbind(needparents,currentdt)
fwrite(as.data.table(cbind(rownames(impgeno),impgeno)),quote = FALSE,file = "genotype.txt",sep=" ")

fwrite(impdt[,c(1:3)],file = "pedigree.txt",sep = " ",quote = FALSE,na = "0")

system("AlphaFamImpute -genotypes genotype.txt -pedigree pedigree.txt -out imp -iothreads 8")

#system("AlphaImpute2 -genotypes genotype.txt -pedigree pedigree.txt -out imp -maxthreads 16")

#system("AlphaImpute2 -genotypes genotype.txt -pedigree pedigree.txt -out imp -ped_only -final_peeling_threshold 0.98 -maxthreads 16")

impgeno= fread("imp.genotypes",sep = " ")

idforimp = impgeno$V1

impgeno = impgeno[,-1]

impgeno = as.matrix(impgeno)

#仅使用单世代填充时注意这个填充缺失基因型的步骤，仅用了单世代均值
impgeno = apply(impgeno,2,function(x){
    
    x[which(!x%in%c(0,1,2))] = NA
    
    avrgeno = mean(x,na.rm = TRUE)

    x[is.na(x)] = avrgeno

    return(x)
})

#If the type of the data is data.fram, the row.names need to be row.names
rownames(impgeno) = idforimp

colnames(impgeno) = lociname

if(exists("HDsib")){

    output$acimp[g+1] = calacforimp(pop = genodt[setdiff(genodt@id,HDsib@id)],impgeno = impgeno)

}else{
    
    output$acimp[g+1] = calacforimp(pop = genodt,impgeno = impgeno)
}

#如果单世代
#if(!exists("impgeno_contain")){
#impgeno_contain = impgeno[rownames(impgeno)%in%genodt@id,]
#}else{
#impgeno_contain = rbind(impgeno_contain,impgeno[rownames(impgeno)%in%genodt@id,])
#}
#if(g>=3){
#impgeno_contain = impgeno_contain[rownames(impgeno_contain) %in% alldt$id,]
#}


if(exists("HDsib")){

    temp_geno = pullSnpGeno(HDsib)

    #如果单世代填充，注销下面这行代码
    genodt = impgeno[rownames(impgeno)%in%genodt@id,]
    
    #如果单世代填充
    #genodt = impgeno_contain
    #如果是单世代填充，注意原本的impgeno_contain也会改变，但是不影响，因为就是要使用原本的高密度子代
    genodt[match(rownames(temp_geno),rownames(genodt)),] = temp_geno

}else{

    #如果单世代填充,只使用注释的代码
    #genodt = impgeno_contain
    genodt = impgeno[rownames(impgeno)%in%genodt@id,]
}

}else{

   genodt = pullMarkerGeno(pop = genodt, markers = snpPenal[[lowpenal]]$id )

   #genodt = pullSnpGeno(genodt)
   #genodt = genodt[,colnames(genodt)%in%snpPenal[[lowpenal]]$id]
}

}else{
    
    if(lowpenal == "0"){
        genodt = 0
    }else{
        genodt = pullSnpGeno(genodt)
    }
    
}

#==============================================================================


#===========calculate the ebv of indivdiuals of valid set================
if(gsmode == "Fam" & is.matrix(genodt)){

    if(imp==FALSE){
        if(g>6){
                needparents= dt_parents[nGeneration%in%c((g-7):(g-1)),]
            }else{
                needparents= dt_parents
            }
    }

    alldt = rbind(needparents, alldt)

    fwrite(alldt,file = "pheno.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")
    
    fwrite(alldt[,1:3],file = "ped.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

    system('hiblup --make-xrm --threads 32 --pedigree ped.csv --add --out hib')

    system(paste("hiblup  --mme --pheno pheno.csv --pheno-pos 6 --xrm hib.PA --vc-priors ",varA(pop),",",(varP(pop)-varA(pop)),
               " --pcg --threads 32 --out hib",sep = ""))

    hiebv<- fread("hib.rand",sep = "\t")

    alldt[,ebv:= hiebv$hib.PA[match(alldt$id,hiebv$ID)]]

    needfam= unique(alldt[nGeneration==g,fam])

    cand_list = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,.(id,suvtime)]
        pheno_t = pheno_a[id%in%pop_test@id,]
        
        pop_v = pop_cand[pop_cand@id%in%pheno_a$id]

        #geno_t= genodt[match(pheno_t$id,rownames(genodt)),]
        #geno_t = geno_t - 1
        #geno_t[geno_t==8] = NA
        
        #geno_v= genodt[match(pop_v@id,rownames(genodt)),]
        #geno_v = geno_v - 1 
        #geno_v[geno_v==8] = NA

        #if(!all(rownames(geno_t)==pheno_t$id)) stop("There is an error in calebv of gsmode fam")

        #ans <- kinship.BLUP(y=pheno_t$suvtime,G.train=geno_t,G.pred=geno_v)

       
        #sebv= ans$g.pred[match(pop_v@id,names(ans$g.pred))]

       genofam= genodt[match(c(pop_v@id,pheno_t$id),rownames(genodt)),]
       #genofam= genodt[match(c(pop_v@id,pheno_t$id,searchParentsFam(pop = pop_v)),rownames(genodt)),]
       
       if(chkpoly(geno_mt = genofam)){
       
       mapdt = makemap()

       genodtput = makeped(z = genofam)

       coutfam = which(needfam == x)

       fwrite(genodtput,file=paste("hib",coutfam,".ped",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

       fwrite(mapdt, file = paste("hib",coutfam,".map",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

       fwrite(alldt,file = paste("pheno",coutfam,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

       system(paste("plink --file hib",coutfam," --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --make-xrm --threads 32 --bfile hib",coutfam," --add --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 6 --xrm hib",coutfam,".GA --vc-priors ",varA(pop),",",(varP(pop)-varA(pop)),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv)[2] = "hib.GA"
       
        sebv = hiebv$hib.GA[match(pop_v@id,hiebv$ID)]

        mebv = alldt$ebv[match(pop_v@mother,alldt$id)]
        
        febv = alldt$ebv[match(pop_v@mother,alldt$id)]

        pop_v@ebv = matrix(sebv + ((mebv+febv)/2),ncol = 1)

       }else{
        
        mebv = alldt$ebv[match(pop_v@mother,alldt$id)]
        
        febv = alldt$ebv[match(pop_v@mother,alldt$id)]

        pop_v@ebv = matrix(((mebv+febv)/2),ncol = 1)
       }
        return(pop_v)

    })

}else if(gsmode == "Pop"){

if(chkpoly(geno_mt = genodt)){
 mapdt = makemap()

 genodtput = makeped(z = genodt)

 fwrite(genodtput,file="hib.ped",col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

 fwrite(mapdt, file = "hib.map",col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

 fwrite(alldt,file = "pheno.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

 system("plink --file hib --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib")
 
 system("hiblup --make-xrm --threads 32 --bfile hib --add --out hib")
 
 system(paste("hiblup --mme --pheno pheno.csv --pheno-pos 6 --xrm hib.GA --vc-priors ",varA(pop),",",(varP(pop)-varA(pop)),
               " --pcg --threads 32 --out hib",sep = ""))

 hiebv<- fread("hib.rand",sep = "\t")

 pop_cand@ebv=matrix(hiebv$hib.GA[match(pop_cand@id,hiebv$ID)],ncol = 1)

}else{
    if(g>6){
                needparents= dt_parents[nGeneration%in%c((g-7):(g-1)),]
            }else{
                needparents= dt_parents
            }

    alldt = rbind(needparents, alldt)

    fwrite(alldt,file = "pheno.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")
    
    fwrite(alldt[,1:3],file = "ped.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

    system('hiblup --make-xrm --threads 32 --pedigree ped.csv --add --out hib')

    system(paste("hiblup  --mme --pheno pheno.csv --pheno-pos 6 --xrm hib.PA --vc-priors ",varA(pop),",",(varP(pop)-varA(pop)),
               " --pcg --threads 32 --out hib",sep = ""))

    hiebv<- fread("hib.rand",sep = "\t")

    pop_cand@ebv=matrix(hiebv$hib.PA[match(pop_cand@id,hiebv$ID)],ncol = 1)
}
 
 

  needfam= unique(alldt[nGeneration==g,fam])

  cand_list = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,id]

        pop_v = pop_cand[pop_cand@id %in% pheno_a ]

        return(pop_v)
    })

}


#=============================================================================================

#numberofid= lapply(cand_list,function(x){
 #   x= x@nInd
  #  return(x)
#})

#numberofid = unlist(numberofid)

#if(sum(is.na(numberofid))>0) save.image(file = paste("errorPOP",r,".rda",sep=""))

if(chkpoly(geno_mt = genodt)){

    acforeach= lapply(cand_list,function(x){
    x= cor(x@ebv,x@gv)
    return(x)
})

acforeach = unlist(acforeach)
output$minacforeach[g+1] = min(acforeach)
output$maxacforeach[g+1] = max(acforeach)
output$acforeach[g+1] = mean(acforeach)

output$acwfall[g+1] = acwfforall(cand_list = cand_list)

}else{
output$minacforeach[g+1] = 0
output$maxacforeach[g+1] = 0
output$acforeach[g+1] = 0
output$acwfall[g+1] = 0
}

output$acall[g+1] = acforall(cand_list = cand_list)




#===================================================================================================


##===================================allcocate the number of candidate within each candidate family=================================================

cand_list = lapply(cand_list,function(x){
    needcand = ncand/nCrosses
    if(chkpoly(geno_mt = genodt)){
    
    x_F = selectWithinFam(x,nInd = round((needcand)*2/3),sex = "F",use = "ebv" )
    
    x_M = selectWithinFam(x,nInd = round((needcand)*1/3),sex = "M",use = "ebv" )

    }else{
    x_F = selectWithinFam(x,nInd = round((needcand)*2/3),sex = "F",use = "rand" )
    
    x_M = selectWithinFam(x,nInd = round((needcand)*1/3),sex = "M",use = "rand" )
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
  
  if(chkpoly(geno_mt = genodt)){

   pop <- ocs(
    pop = candidate,
    nCrosses = nCrosses,
    nProgenyPerCross = nProgenyPerCross,
    nFemalesMax = nDam,
    nMalesMax = nSire,
    equalizeFemaleContributions =TRUE,
    equalizeMaleContributions = TRUE,
    targetDegree = 45,
    use = "ebv"
  )

  }else{
    pop <- ocs(
    pop = candidate,
    usePed = TRUE,
    ped_parents = ped_ne,
    nCrosses = nCrosses,
    nProgenyPerCross = nProgenyPerCross,
    nFemalesMax = nDam,
    nMalesMax = nSire,
    equalizeFemaleContributions =TRUE,
    equalizeMaleContributions = TRUE,
    targetDegree = 45,
    use = "ebv"
  )
  }
}
}
warnings()

fwrite(output,paste(gsmode,lowpenal,Im,nReference,r,".csv",sep = ""),sep = ",")