snpNum = lapply(lowdensity,function(x){
  chr_snp_num= data.frame(chr = c(1:nChr),num = round(x/nChr), length = ChrSize )

  reallnum = sum(chr_snp_num$num) - x

  if(reallnum > 0){
   chrid = sample(chr_snp_num$chr,size = reallnum)

   chr_snp_num$num[chr_snp_num$chr%in%chrid] = chr_snp_num$num[chr_snp_num$chr%in%chrid] - 1

  }else if(reallnum < 0){
    chrid = sample(chr_snp_num$chr,size = abs(reallnum))
    chr_snp_num$num[chr_snp_num$chr%in%chrid] = chr_snp_num$num[chr_snp_num$chr%in%chrid] + 1
  }

  return(chr_snp_num)
})

snpPenal= lapply(snpNum,greedyChooseLoci,genos = pullSnpGeno(pop_founder),map = getSnpMap())

names(snpPenal) = as.character(lowdensity)
