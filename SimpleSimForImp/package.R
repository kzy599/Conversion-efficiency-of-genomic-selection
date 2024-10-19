library(dplyr)
library(purrr)
library(AlphaSimR)
library(rrBLUP)
library(data.table)
library(visPedigree)
library(optiSel)
library(sampling)
library(AlphaMME)
library(AlphaLearn)
 #make .ped file required by plink
  makeped = function(z){
    z[z==2]= 22
    z[z==1]= 12
    z[z==0]= 11
    z[z==9]= -9
    needname = rownames(z)
    z = as.data.table(cbind(needname,z))
    return(z)
  }

   #make .map file required by plink 
  makemap = function(...){
  mapck = getSnpMap()
  #if(lowpenal %in% as.character(lowdensity)){
  #  if(imp == FALSE){
   #  mapck = mapck[snpPenal[[lowpenal]]$id,]
   # }
  #}
  mapck$id = mapck$chr
  mapck$chr = rownames(mapck)
  mapck$site = mapck$pos
  mapck$pos = rep(0,nrow(mapck))
  return(mapck)
  }



  makeped_cal = function(pop){
    
    z=pullSnpGeno(pop)
    
    z[z==2]=22
    z[z==1]=12
    z[z==0]=11

    #z = as.data.table(cbind(pop@id,z))
    z = as.data.table(cbind(pop@id,z))
    return(z)
  }

makemap_cal = function(...){
  mapck = getSnpMap()
  mapck$id = mapck$chr
  mapck$chr = rownames(mapck)
  mapck$site = mapck$pos
  mapck$pos = rep(0,nrow(mapck))
  return(mapck)
  }


calNeped = function(ped,keep){
  pedig <- prePed(ped)
  pKin   <- pedIBD(pedig, keep.only = keep)
  Summary <- summary(pedig)
  id     <- keep
  x      <- Summary[Indiv %in% id]$equiGen
  N      <- length(x)
  n      <- (matrix(x, N, N, byrow = TRUE) + matrix(x, N, N, byrow = FALSE)) / 2
  deltaC <- 1 - (1 - pKin[id, id]) ^ (1 / n)
  Ne   <- 1 / (2 * mean(deltaC))
  return(Ne)
  
}

calInbped = function(ped,keep){
  
  Pedig <- prePed(ped, keep=keep)

  Res   <- pedInbreeding(Pedig)
  
  inbreeding <- mean(Res$Inbr[Res$Indiv %in% keep])
  
  return(inbreeding)
  
}

calinbya2 <- function(snp_012_dt,snpfre){
  calformula <- function(x, snpfre){
    y <- mean((x^2 - ((1 + 2 * snpfre) * x) + 2 * (snpfre^2)) / (2 * snpfre * (1 - snpfre)), na.rm = TRUE)
    return(y)
  }
  inb <- apply(snp_012_dt, 1, calformula,snpfre = snpfre)
  return(inb)
}



scoreGreedy <- function(start, end, maf, pos){
  return(
    maf * (end - start - abs((2 * pos) - end - start))
  )
}
greedyChooseLoci <- function(num, genos, map){
  colnames(num) <- c("chr", "num", "length")
  # calc minor allele freq
  maf <- (colSums(genos) / (2 * nrow(genos)))
  maf[maf > 0.5] <- 1 - maf[maf > 0.5]
  map <- map %>% left_join(data.frame(maf = maf, id = names(maf)), by = "id") 
  
  panel <- data.frame()
  for(i in 1:nrow(num)){ # for each chr
    cands <- map %>% filter(chr == num$chr[i])
    if(nrow(cands) < num$num[i]){
      stop("Not enough SNPs in chromosome ", num$chr[i])
    }
    # calculate scores at the start
    cands <- cands %>% 
      mutate(score = scoreGreedy(start = 0, end = num$length[i], maf = maf, pos = pos),
             lastStart = 0,
             lastEnd = num$length[i])
    numSNPs <- 0
    while(TRUE){ # for each desired SNP
      # choose SNP
      temp <- which.max(cands$score)
      chosen <- cands %>% slice(temp)
      panel <- panel %>% bind_rows(chosen)
      cands <- cands[-temp,]
      numSNPs <- numSNPs + 1
      if(numSNPs == num$num[i]) break
      
      # recalculate scores for affected SNPs
      # only those whose interval used for the last calculation are affected
      # note that there is only one interval to update b/c the new SNP can only
      # have been in one interval
      toUpdate <- cands %>% filter(lastStart < chosen$pos & lastEnd > chosen$pos) %>%
        select(lastStart, lastEnd) %>% distinct() # this is some bs to help vectorize operations for R
      tempBool <- cands$lastStart == toUpdate$lastStart & cands$pos < chosen$pos
      cands$score[tempBool] <- scoreGreedy(start = toUpdate$lastStart,
                                           end = chosen$pos,
                                           maf = cands$maf[tempBool],
                                           pos = cands$pos[tempBool])
      cands$lastEnd[tempBool] <- chosen$pos
      tempBool <- cands$lastStart == toUpdate$lastStart & cands$pos > chosen$pos
      cands$score[tempBool] <- scoreGreedy(start = chosen$pos,
                                           end = toUpdate$lastEnd,
                                           maf = cands$maf[tempBool],
                                           pos = cands$pos[tempBool])
      cands$lastStart[tempBool] <- chosen$pos
    }
  }
  return(panel %>% arrange(chr, site))
}

calcandnum= function(facc,ncand){
  neednum = ncand
  if(ncol(facc)!=2){
    warnings("The facc must be a type of data.table with 2 cols that first col is the name of family and the second is accuracy")
  }
  if(!is.data.table(facc)){
    facc= as.data.table(facc)
  }
  
  colnames(facc) = c("fname","acc")
  
  if(!is.numeric(facc$acc)){
    facc$acc = as.numeric(facc$acc)
  }
  
  setorder(facc,-acc)
  fnum = data.table(fname = facc$fname, candnum = numeric(length = length(facc$fname)))
  #========calculating the number of individuals within each family==============
  
  nfam = length(facc$fname)
  
  nloop = 0
  
  while(TRUE){
    if(ncand>nfam){
      candnum = (ncand/sum(facc$acc[1:nfam]))*facc$acc[1:nfam]
      candnum = floor(candnum)  
    }else if(ncand == 1){
      candnum = c(1,numeric(length = nfam-ncand))
    }else{
      candnum = (ncand/sum(facc$acc[1:ncand])) * facc$acc[1:ncand]
      #There is an bug in floor if the type of candnum is numeric and equal 1,
      #floor(candnum) = 0
      #as.intger(cannum) = 0
      candnum = floor(candnum) 
      candnum = c(candnum,numeric(length = nfam-ncand))
    }
    
    fnum$candnum= fnum$candnum + candnum
    
    ncand = ncand - sum(candnum)
    
    nloop = nloop + 1

    if(nloop>1000){
      stop("There is an error in calcandnum")
    }

    if(sum(fnum$candnum) == neednum) break
  }
  #========================ending================================================
return(fnum)
}

mergelist = function(x){

  for(l in 1:length(x)){
   if(l==1){
    canlist = x[[l]]
   }else{
    canlist = rbind(canlist,x[[l]])
   }
  }

  return(canlist)
}

reimpute = function(impgeno){

  for(rn in 1:nrow(impgeno)){
    x = impgeno[rn,]
    colid= names(which(x==9))
    if(length(colid) != 0 ){
     for(rimp in 1:length(colid)){

        needg = impgeno[,colid[rimp]]

        x[colid[rimp]] = round(mean(needg[needg!=9]))
    }
    }
    if(rn == 1){
        ckck = x
    }else{
        ckck = rbind(ckck,x)
    }
}
return(ckck)
}


acforall = function(cand_list){
    allcand = mergePops(cand_list)
    ac = cor(allcand@gv,allcand@ebv)
    return(ac)
}

calacforimp = function(pop,impgeno){
  
  trueGenos = pullSnpGeno(pop)
  imputeCalls = impgeno[match(rownames(trueGenos),rownames(impgeno)),]

  #remove non-variant SNP
  temp <- (apply(trueGenos, 2, n_distinct) > 1) & (apply(imputeCalls, 2, n_distinct) > 1)
	
  trueGenos <- trueGenos[,temp]
	
  imputeCalls <- imputeCalls[,temp]
	
  rm(temp)


  acimp= mean(sapply(1:ncol(trueGenos), function(x) cor(trueGenos[,x], imputeCalls[,x])))

  return(acimp)
}


searchParentsFam = function(pop){
  count = 0
  while(TRUE){

    if(count == 0){
         pid = unique(c(pop@mother,pop@father))
    }else{
         pid = unique(c(famdt$dam,famdt$sir))
    }

   pfam = unique(needparents[id %in% pid,fam])

   if(length(pfam)==1){if(is.na(pfam)) break }

   if(!any(pfam %in% unique(alldt$fam)))break

   famdt = alldt[fam %in% pfam ,]

   if(count == 0){
     out = famdt
     
   }else{
    out = rbind(out,famdt)
   }

   count = count +1

   if(count > 2000) stop("There is an error in searchParentsFam() ")

  }
  
  if(count == 0){
    outid = c()
  }else{
    outid = out$id
  }

  return(outid)

}

caleachfamsnp = function(needfam){

  nSnp = function(genofam){
    snpreq = apply(genofam,2,function(x){
    freq = sum(x,na.rm = TRUE)/(2*length(x))
    return(freq)
     })
    snpreq[snpreq>0.5] = 1 - snpreq[snpreq>0.5]
    MAF1= sum(snpreq>0.01)
    MAF5 = sum(snpreq>0.05)
    mafdt= data.table(MAF1,MAF5)
    colnames(mafdt) = c("MAF1","MAF5")
    return(mafdt)
}
   for(ff in 1:length(needfam)){
    
    needid = alldt[fam==needfam[ff],id]
    
    genofam= genodt[match(needid,rownames(genodt)),]
    
    maf_temp = nSnp(genofam)
    
    genfam_temp = data.table(nGeneration = g,fam = needfam[ff])
    
    colnames(genfam_temp) = c("nGeneration","fam")
    
    out_temp = cbind(genfam_temp,maf_temp)

    if(ff == 1){out = out_temp}else{out = rbind(out,out_temp)}

   }

   return(out)
}




calacforeach = function(cand_list,pname,ttype){

if(ttype == "BW") tcount = 1

if(ttype == "ST") tcount = 2

acforeach= lapply(cand_list,function(x){

    if(ncol(x@ebv)==2){ecount = tcount}else{ecount = 1}
    x= cor(x@ebv[,ecount],x@gv[,tcount])
    return(x)
})

acforeach = unlist(acforeach)
pop_temp = mergePops(cand_list)

if(ncol(pop_temp@ebv)==2){ecount = tcount}else{ecount = 1}
accross= cor(pop_temp@ebv[,ecount],pop_temp@gv[,tcount])

out = data.table(pname = pname,ttype = ttype,minac = min(acforeach),maxac = max(acforeach),meanac = mean(acforeach),accross = accross)

colnames(out) = c("pname","ttype","min","max","mean","accross")

return(out)
}

calactrait = function(hiebv,pname,ttype){

  pop_cand@ebv=matrix(hiebv$hib.GA[match(pop_cand@id,hiebv$ID)],ncol = 1)

  needfam= unique(alldt[nGeneration==g,fam])

  cand_list = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,id]

        pop_v = pop_cand[pop_cand@id %in% pheno_a ]

        return(pop_v)
    })

  acout = calacforeach(cand_list = cand_list, pname = pname,ttype = ttype)

  return(acout)
}



calidenticalfamsnp = function(needfam){

  calnsnp = function(genofam){
    snpreq = apply(genofam,2,function(x){
    freq = sum(x,na.rm = TRUE)/(2*length(x))
    return(freq)
     })
    snpreq[snpreq>0.5] = 1 - snpreq[snpreq>0.5]

    return(snpreq)
    }

  nidenticalSnp = function(genofam_cand,genofam_test){
    
snpreq_cand = calnsnp(genofam = genofam_cand)

snpreq_test = calnsnp(genofam = genofam_test)

locicand1  = names(snpreq_cand)[which(snpreq_cand>0.01)]

locicand5  = names(snpreq_cand)[which(snpreq_cand>0.05)] 

locitest1  = names(snpreq_test)[which(snpreq_test>0.01)]

locitest5  = names(snpreq_test)[which(snpreq_test>0.05)]

allloci1 = c(locicand1,locitest1)

allloci5 = c(locicand5,locitest5)

alllociall1 = c(locicand1,locialltest1)

alllociall5 = c(locicand5,locialltest5)

mafdt= data.table(MAF1cand = length(locicand1),MAF5cand =length(locicand5),
                  MAF1test = length(locitest1),MAF5test =length(locitest5),  
                  MAF1testall = length(locialltest1),MAF5testall =length(locialltest5),  
                  nsameloci1 = length(allloci1)- n_distinct(allloci1),
                  nsameloci5 = length(allloci5) - n_distinct(allloci5),
                  nsamelociall1 = length(alllociall1)- n_distinct(alllociall1),
                  nsamelociall5 = length(alllociall5) - n_distinct(alllociall5))
  colnames(mafdt) = c("MAF1cand", "MAF5cand", "MAF1test"," MAF5test",
                       "MAF1alltest"," MAF5alltest",
                       "nsameloci1","nsameloci5","nsamelociall1","nsamelociall5")
  return(mafdt)

}

alltestid = alldt[!id %in% pop_cand@id,id]

genofam_alltest = genodt[match(alltestid,rownames(genodt)),]

snpreq_alltest = calnsnp(genofam = genofam_alltest)


locialltest1  = names(snpreq_alltest)[which(snpreq_alltest>0.01)]

locialltest5  = names(snpreq_alltest)[which(snpreq_alltest>0.05)]


   for(ff in 1:length(needfam)){
    
    needid = alldt[fam==needfam[ff],id]

    candid = needid[needid %in% pop_cand@id]

    testid = needid[needid %in% pop_test@id]
    
    genofam_cand = genodt[match(candid,rownames(genodt)),]
    
    genofam_test = genodt[match(testid,rownames(genodt)),]

    maf_temp = nidenticalSnp(genofam_cand = genofam_cand,genofam_test = genofam_test)
    
    genfam_temp = data.table(nGeneration = g,fam = needfam[ff])
    
    colnames(genfam_temp) = c("nGeneration","fam")
    
    out_temp = cbind(genfam_temp,maf_temp)

    if(ff == 1){out = out_temp}else{out = rbind(out,out_temp)}

   }

   return(out)
}

searchrealationsFam = function(genodt,pop,pheno_t){

 
 genotemp = genodt[!rownames(genodt)%in% pheno_t$id,]
 
 idtemp = rownames(genotemp)
 
 Gmat = calcG(genotemp)

 rownames(Gmat) = idtemp
 colnames(Gmat) = idtemp
 
 rowbool = rownames(Gmat) %in% pop@id
 
 colbool = !colnames(Gmat) %in% pop@id

 Gmat = Gmat[rowbool, colbool]

 relcoe = apply(Gmat,2,function(x){
  relcoe = mean(x,na.rm = TRUE)
  return(relcoe)
 })

 namebool = which(relcoe >=0.25 & relcoe <= 0.75)

 idfinal = names(relcoe)[namebool]

 return(idfinal)

}