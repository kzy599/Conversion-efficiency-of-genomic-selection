snpMap <- getSnpMap()
snpMap$pos<- snpMap$pos* (2.6 * 10^9) / 44#单条染色体的物理距离
numloci<- c(100, 250, 500, 758, 1000, 2000, 5000, 10000)

#orignial design
empNum <- lapply(numloci, function(x){
  dfOut <- data.frame(chr = 1:44, num = round(x / 44),
                      length = (2.6 * 10^9) / 44) # proportional to chr length
  # account for rounding error (in same way as above, so numbers match w/ largest panel)
  dfOut$num[which.max(dfOut$num)] <- dfOut$num[which.max(dfOut$num)] + x - sum(dfOut$num)
  
  # make sure enough loci and reallocate as needed
  toRealloc <- c()
  numSnps <- c()
  for(i in 1:nrow(dfOut)){ # for each chr
    tmp <- sum(snpMap$chr == dfOut$chr[i]) # number of snps available
    numSnps <- c(numSnps, tmp)
    toRealloc <- c(toRealloc, dfOut$num[i] - tmp)
    rm(tmp)
  }
  if(any(toRealloc > 0)){
    if(sum(toRealloc) > 0) stop("Not enough SNPs")
    warning("Reallocating: Not enough SNPs in one or more chromosomes")
    totalSnps <- sum(dfOut$num) # for checking
    numToReall <- sum(toRealloc[toRealloc > 0])
    # distribute proportion to excess SNPs, making sure the total matches
    # and the limit on all other chr isn't exceeded
    ch <- c()
    for(i in 1:length(toRealloc)) if(toRealloc[i] < 0) ch <- c(ch, rep(i, -1 * toRealloc[i]))
    addSnps <- sample(ch, numToReall, replace = FALSE)
    for(i in 1:length(toRealloc)){
      if(toRealloc[i] > 0){
        toRealloc[i] <- toRealloc[i] * -1
      } else {
        toRealloc[i] <- sum(addSnps == i)
      }
    }
    dfOut$num <- dfOut$num + toRealloc
    if(totalSnps != sum(dfOut$num)) stop("Error reallocating. Totals don't match.")
  }
  
  return(dfOut)
})

#More advanced design   author: Kang ZY
snpNum = lapply(numloci,function(x){
  chr_snp_num= data.frame(chr = c(1:44),num = round(x/44), length =(2.6 * 10^9) / 44 )

  reallnum = sum(chr_snp_num$num) - x

  if(reallnum > 0){
   chrid = sample(chr_snp_num$chr,size = reallnum)

   chr_snp_num[chr%in%chrid,num:=(num - 1)]

  }else if(reallnum < 0){
    chrid = sample(chr_snp_num$chr,size = abs(reallnum))
    chr_snp_num[chr%in%chrid,num:=(num+1)]
  }

  for(i in 1:nrow(chr_snp_num)){ # for each chr
    tmp <- sum(snpMap$chr == chr_snp_num$chr[i]) # number of snps available
    numSnps <- c(numSnps, tmp)
    toRealloc <- c(toRealloc, chr_snp_num$num[i] - tmp)
    rm(tmp)
  }
  if(any(toRealloc > 0)){
    if(sum(toRealloc) > 0) stop("Not enough SNPs")
    warning("Reallocating: Not enough SNPs in one or more chromosomes")
    totalSnps <- sum(chr_snp_numt$num) # for checking
    numToReall <- sum(toRealloc[toRealloc > 0])
    # distribute proportion to excess SNPs, making sure the total matches
    # and the limit on all other chr isn't exceeded
    ch <- c()
    for(i in 1:length(toRealloc)) if(toRealloc[i] < 0) ch <- c(ch, rep(i, -1 * toRealloc[i]))
    addSnps <- sample(ch, numToReall, replace = FALSE)
    for(i in 1:length(toRealloc)){
      if(toRealloc[i] > 0){
        toRealloc[i] <- toRealloc[i] * -1
      } else {
        toRealloc[i] <- sum(addSnps == i)
      }
    }
    chr_snp_num$num <- chr_snp_num$num + toRealloc
    if(totalSnps != sum(chr_snp_num$num)) stop("Error reallocating. Totals don't match.")
  }

  return(chr_snp_num)
})


numneed<- empNum[[1]]

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


#another format for loop author: Kang ZY
greedyChooseLoci <- function(num, genos, map){
  colnames(num) <- c("chr", "num", "length")
  # calc minor allele freq
  maf <- (colSums(genos) / (2 * nrow(genos)))
  maf[maf > 0.5] <- 1 - maf[maf > 0.5]
  map <- map %>% left_join(data.frame(maf = maf, id = names(maf)), by = "id") 
  
   panel_chr = lapply(c(1:44),function(x){
    panel <- data.frame()
    cands <- map %>% filter(chr == x)
    if(nrow(cands) < num$num[x]){
      stop("Not enough SNPs in chromosome ", num$chr[x])
    }
    # calculate scores at the start
    cands <- cands %>% 
      mutate(score = scoreGreedy(start = 0, end = num$length[x], maf = maf, pos = pos),
             lastStart = 0,
             lastEnd = num$length[x])
    numSNPs <- 0
    while(TRUE){ # for each desired SNP
      # choose SNP
      temp <- which.max(cands$score)
      chosen <- cands %>% slice(temp)
      panel <- panel %>% bind_rows(chosen)
      cands <- cands[-temp,]
      numSNPs <- numSNPs + 1
      if(numSNPs == num$num[x]) break
      
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
    return(panel %>% arrange(chr, site))
  })
panel = panel_chr %>% reduce(full_join)
return(panel %>% arrange(chr, site))
}


ckckc= greedyChooseLoci(numneed,genos = pullSnpGeno(pop_founder),map = snpMap)