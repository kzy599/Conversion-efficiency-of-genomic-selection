
pop_select <- selectWithinFam(pop = pop_cand,nInd = 10,use = "rand",simParam = SP)

keep<- pop_select@id

#pop_LD<- pop_select
pop_Inb <-mergePops(list(pop_founder,pop_select)) 

if(!exists("snpfre_v")){
  snp_012_dt = pullSnpGeno(pop_founder)
  snpfre_v <- apply(snp_012_dt, 2, function(x){
  single_snpfre_s <- sum(x, na.rm = TRUE)/(2*length(x))
  return(single_snpfre_s)
})
}

snp_012_dt = pullSnpGeno(pop_select)


inbreeding_ya2 = mean(calinbya2(snp_012_dt = snp_012_dt, snpfre = snpfre_v))

ped_calparameters <- rbind(ped_ne,alldt[,1:3])

inbreeding_ped = calInbped(ped = ped_calparameters,keep = keep)

if(g==0){
    Ne = 75
}else{
   Ne = calNeped(ped = ped_calparameters,keep = keep)
}

#

peddt = makeped_cal(pop = pop_select)

mapdt = makemap_cal()

fwrite(peddt,file="popLD.ped",col.names = FALSE,sep = "\t")

fwrite(mapdt, file = "popLD.map",col.names = FALSE,sep = " ")

system("plink --file popLD  --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out popLD")

system("gcta64 --bfile popLD --autosome-num 44 --ld-score --ld-wind 10000 --ld-rsq-cutoff 0.01 --out popLD")

LDfile<- fread(input="popLD.score.ld",sep = " ")
LD <- mean(LDfile$mean_rsq)
LDscore <- mean(LDfile$ldscore)


ped_calparameters<- visPedigree::tidyped(ped_calparameters,cand = keep)

output$nCoancestor[g+1]<- length(ped_calparameters[Gen==1,Ind])

output$accuracy[g+1] <- cor(candidate@gv,candidate@ebv)

output$inbreeding[g+1] <- inbreeding_ped

output$inbreeding_plink[g+1] <- inbreeding_ya2

output$LD[g+1] <- LD

output$LDscore[g+1] <- LDscore

output$Ne[g+1] <- Ne

output$Va[g+1] <- varA(pop,simParam = SP)

output$Vg[g+1] <- varG(pop)

output$Vp[g+1] <- varP(pop)

output$mean_pheno[g+1] <- mean(pop@pheno)

output$mean_gv[g+1] <- mean(pop@gv)

output$genicVa[g+1] <- genicVarA(pop)

output$genicVg[g+1] <- genicVarG(pop)

output$h2[g+1] = varA(pop)/varP(pop)