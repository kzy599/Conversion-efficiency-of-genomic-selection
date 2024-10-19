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
library(stringr)
 #make .ped file required by plink
  makeped = function(z){
    z[!z%in%c(0,1,2)] = -9
    z[z==2]= 22
    z[z==1]= 12
    z[z==0]= 11
    needname = rownames(z)
    z = as.data.table(cbind(needname,z))
    return(z)
  }

   #make .map file required by plink 
  makemap = function(...){
  mapck = getSnpMap()
  if(lowpenal %in% as.character(lowdensity)){
    if(imp == FALSE){
     mapck = mapck[snpPenal[[lowpenal]]$id,]
    }
  }
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
    ac = cor(allcand@gv[,1],allcand@ebv)
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

chkpoly = function(geno_mt){
  if(!is.matrix(geno_mt)){

    Boolre = FALSE

  }else{

    count_temp = apply(geno_mt,2,function(x){
    
    return(n_distinct(x))
    
  })
  
  #if poly
  Boolre = sum(count_temp > 1) >= 10 

  }
  return(Boolre)
}

source("ocs.R")