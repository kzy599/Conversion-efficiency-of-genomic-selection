rm(list = ls())
gc()

allfiles = list.files()

r = 1

ngenpop = "allpop3"

nlen = nchar(ngenpop)

ngenpop = paste(ngenpop,r,".rda",sep = "")

if(any(allfiles == ngenpop)){load(ngenpop)}else{source("creatpop.R")}

source("package.R")


rr = FALSE

qc = FALSE

if(rr){
    qc = FALSE
}
#===================================four generations data for pop============================
#fourpop_ac = data.table(pname = "fourpop",minac = min(acforeach),maxac = max(acforeach),meanac = mean(acforeach))

#colnames(fourpop_ac) = c("pname","min","max","mean")

fourpop_ac_BW= calactrait(hiebv = hiebv_BW,pname = "fourpop",ttype = "BW")

fourpop_ac_ST = calactrait(hiebv = hiebv_ST,pname = "fourpop",ttype = "ST")
#==================================================================================================





#===================================one generation data for pop============================
onepopgeno = genodt[rownames(genodt)%in%alldt[nGeneration==g,id],]

onepopgenoput = makeped(z = onepopgeno)

 fwrite(onepopgenoput,file="hib.ped",col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

 fwrite(mapdt, file = "hib.map",col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

 fwrite(alldt[nGeneration == g,],file = "pheno.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

 system("plink --file hib --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib")
 
 system("hiblup --make-xrm --threads 32 --bfile hib --add --out hib")
 
 system(paste("hiblup --mme --pheno pheno.csv --pheno-pos 6 --xrm hib.GA --vc-priors ",varA(pop)[2,2],",",(varP(pop)[2,2]-varA(pop)[2,2]),
               " --pcg --threads 32 --out hib",sep = ""))

 hiebv_ST<- fread("hib.rand",sep = "\t")

 system(paste("hiblup --mme --pheno pheno.csv --pheno-pos 7 --xrm hib.GA --vc-priors ",varA(pop)[1,1],",",(varP(pop)[1,1]-varA(pop)[1,1]),
               " --pcg --threads 32 --out hib",sep = ""))

 hiebv_BW<- fread("hib.rand",sep = "\t")

 onepop_ac_BW= calactrait(hiebv = hiebv_BW,pname = "onepop",ttype = "BW")

 onepop_ac_ST = calactrait(hiebv = hiebv_ST,pname = "onepop",ttype = "ST")

 #pop_cand@ebv=matrix(hiebv$hib.GA[match(pop_cand@id,hiebv$ID)],ncol = 1)

  #needfam= unique(alldt[nGeneration==g,fam])

  #cand_list = lapply(needfam,function(x){
   #     pheno_a = alldt[fam==x,id]

    #    pop_v = pop_cand[pop_cand@id %in% pheno_a ]

     #   return(pop_v)
    #})

#onepop_ac = calacforeach(cand_list = cand_list, pname = "onepop")
#=====================================================================================================================================





#==================================fam mode============================================

#if(imp==FALSE){
        if(g>7){
                needparents= dt_parents[nGeneration%in%c((g-8):(g-1)),]
            }else{
                needparents= dt_parents
            }
 #   }

    alldt_pblup = rbind(needparents, alldt)

    fwrite(alldt_pblup,file = "pheno.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")
    
    fwrite(alldt_pblup[,1:3],file = "ped.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

    system('hiblup --make-xrm --threads 32 --pedigree ped.csv --add --out hib')

    system(paste("hiblup  --mme --pheno pheno.csv --pheno-pos 6 --xrm hib.PA --vc-priors ",varA(pop)[2,2],",",(varP(pop)[2,2]-varA(pop)[2,2]),
               " --pcg --threads 32 --out hib",sep = ""))

    hiebv_ST<- fread("hib.rand",sep = "\t")


    system(paste("hiblup  --mme --pheno pheno.csv --pheno-pos 7 --xrm hib.PA --vc-priors ",varA(pop)[1,1],",",(varP(pop)[1,1]-varA(pop)[1,1]),
               " --pcg --threads 32 --out hib",sep = ""))

    hiebv_BW<- fread("hib.rand",sep = "\t")

    alldt_pblup[,ebv_ST:= hiebv_ST$hib.PA[match(alldt_pblup$id,hiebv_ST$ID)]]

    alldt_pblup[,ebv_BW:= hiebv_BW$hib.PA[match(alldt_pblup$id,hiebv_BW$ID)]]

    needfam= unique(alldt[nGeneration==g,fam])

if(!rr){
    #=============one generation no qc=========================

    cand_list = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,.(id,suvtime)]
        pheno_t = pheno_a[id%in%pop_test@id,]
        
        pop_v = pop_cand[pop_cand@id%in%pheno_a$id]

       genofam= genodt[match(c(pop_v@id,pheno_t$id),rownames(genodt)),]
       #genofam= genodt[match(c(pop_v@id,pheno_t$id,searchParentsFam(pop = pop_v)),rownames(genodt)),]

       mapdt = makemap()

       genodtput = makeped(z = genofam)

       dtfam = alldt[id %in% rownames(genofam),]

       coutfam = which(needfam == x)

       fwrite(genodtput,file=paste("hib",coutfam,".ped",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

       fwrite(mapdt, file = paste("hib",coutfam,".map",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

       fwrite(dtfam,file = paste("pheno",coutfam,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

       system(paste("plink --file hib",coutfam," --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --make-xrm --threads 32 --bfile hib",coutfam," --add --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 6 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[2,2],",",(varP(pop)[2,2]-varA(pop)[2,2]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_ST<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_ST)[2] = "hib.GA"

       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 7 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[1,1],",",(varP(pop)[1,1]-varA(pop)[1,1]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_BW<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_BW)[2] = "hib.GA"
       
        sebv_BW = hiebv_BW$hib.GA[match(pop_v@id,hiebv_BW$ID)]

        mebv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]
        
        febv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]

        sebv_ST = hiebv_ST$hib.GA[match(pop_v@id,hiebv_ST$ID)]

        mebv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]
        
        febv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]

        pop_v@ebv = cbind(sebv_BW + ((mebv_BW+febv_BW)/2),sebv_ST + ((mebv_ST+febv_ST)/2))
        
        return(pop_v)

    })

    #=============one generation qc=========================
    if(qc){
      cand_list_qc = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,.(id,suvtime)]
        pheno_t = pheno_a[id%in%pop_test@id,]
        
        pop_v = pop_cand[pop_cand@id%in%pheno_a$id]

       genofam= genodt[match(c(pop_v@id,pheno_t$id),rownames(genodt)),]
       #genofam= genodt[match(c(pop_v@id,pheno_t$id,searchParentsFam(pop = pop_v)),rownames(genodt)),]

       mapdt = makemap()

       genodtput = makeped(z = genofam)

       dtfam = alldt[id %in% rownames(genofam),]

       coutfam = which(needfam == x)

       fwrite(genodtput,file=paste("hib",coutfam,".ped",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

       fwrite(mapdt, file = paste("hib",coutfam,".map",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

       fwrite(dtfam,file = paste("pheno",coutfam,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

       system(paste("plink --file hib",coutfam," --geno --mind --maf 0.01 --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --make-xrm --threads 32 --bfile hib",coutfam," --add --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 6 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[2,2],",",(varP(pop)[2,2]-varA(pop)[2,2]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_ST<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_ST)[2] = "hib.GA"

       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 7 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[1,1],",",(varP(pop)[1,1]-varA(pop)[1,1]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_BW<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_BW)[2] = "hib.GA"
       
        sebv_BW = hiebv_BW$hib.GA[match(pop_v@id,hiebv_BW$ID)]

        mebv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]
        
        febv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]

        sebv_ST = hiebv_ST$hib.GA[match(pop_v@id,hiebv_ST$ID)]

        mebv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]
        
        febv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]

        pop_v@ebv = cbind(sebv_BW + ((mebv_BW+febv_BW)/2),sebv_ST + ((mebv_ST+febv_ST)/2))
        
        return(pop_v)

    })
    }
     #=============four generation no qc=========================
    cand_list_four = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,.(id,suvtime)]
        pheno_t = pheno_a[id%in%pop_test@id,]
        
        pop_v = pop_cand[pop_cand@id%in%pheno_a$id]

       #genofam= genodt[match(c(pop_v@id,pheno_t$id),rownames(genodt)),]
       genofam= genodt[match(c(pop_v@id,pheno_t$id,searchParentsFam(pop = pop_v)),rownames(genodt)),]

       mapdt = makemap()

       genodtput = makeped(z = genofam)

       dtfam = alldt[id %in% rownames(genofam),]

       coutfam = which(needfam == x)

       fwrite(genodtput,file=paste("hib",coutfam,".ped",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

       fwrite(mapdt, file = paste("hib",coutfam,".map",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

       fwrite(dtfam,file = paste("pheno",coutfam,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

       system(paste("plink --file hib",coutfam," --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --make-xrm --threads 32 --bfile hib",coutfam," --add --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 6 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[2,2],",",(varP(pop)[2,2]-varA(pop)[2,2]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_ST<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_ST)[2] = "hib.GA"

       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 7 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[1,1],",",(varP(pop)[1,1]-varA(pop)[1,1]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_BW<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_BW)[2] = "hib.GA"
       
        sebv_BW = hiebv_BW$hib.GA[match(pop_v@id,hiebv_BW$ID)]

        mebv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]
        
        febv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]

        sebv_ST = hiebv_ST$hib.GA[match(pop_v@id,hiebv_ST$ID)]

        mebv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]
        
        febv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]

        pop_v@ebv = cbind(sebv_BW + ((mebv_BW+febv_BW)/2),sebv_ST + ((mebv_ST+febv_ST)/2))
        
        return(pop_v)

    })
    

    cand_list_rel = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,.(id,suvtime)]
        pheno_t = pheno_a[id%in%pop_test@id,]
        
        pop_v = pop_cand[pop_cand@id%in%pheno_a$id]

       #genofam= genodt[match(c(pop_v@id,pheno_t$id),rownames(genodt)),]
       genofam= genodt[match(c(pop_v@id,pheno_t$id,searchrealationsFam(genodt = genodt,pop = pop_v,pheno_t = pheno_t)),rownames(genodt)),]

       mapdt = makemap()

       genodtput = makeped(z = genofam)

       dtfam = alldt[id %in% rownames(genofam),]

       coutfam = which(needfam == x)

       fwrite(genodtput,file=paste("hib",coutfam,".ped",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

       fwrite(mapdt, file = paste("hib",coutfam,".map",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

       fwrite(dtfam,file = paste("pheno",coutfam,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

       system(paste("plink --file hib",coutfam," --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --make-xrm --threads 32 --bfile hib",coutfam," --add --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 6 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[2,2],",",(varP(pop)[2,2]-varA(pop)[2,2]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_ST<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_ST)[2] = "hib.GA"

       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 7 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[1,1],",",(varP(pop)[1,1]-varA(pop)[1,1]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_BW<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_BW)[2] = "hib.GA"
       
        sebv_BW = hiebv_BW$hib.GA[match(pop_v@id,hiebv_BW$ID)]

        mebv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]
        
        febv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]

        sebv_ST = hiebv_ST$hib.GA[match(pop_v@id,hiebv_ST$ID)]

        mebv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]
        
        febv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]

        pop_v@ebv = cbind(sebv_BW + ((mebv_BW+febv_BW)/2),sebv_ST + ((mebv_ST+febv_ST)/2))
        
        return(pop_v)

    })


    #===============one generation fam with extra data======================

    cand_list_extra_50 = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,.(id,suvtime)]
        pheno_t = pheno_a[id%in%pop_test@id,]
        
        pop_v = pop_cand[pop_cand@id%in%pheno_a$id]

       genofam= genodt[match(c(pop_v@id,sample(pheno_t$id,size = 50)),rownames(genodt)),]
       #genofam= genodt[match(c(pop_v@id,pheno_t$id,searchParentsFam(pop = pop_v)),rownames(genodt)),]
       
       #dtextra = alldt_extra[fam==x,]

       #idtemp= sample(dtextra$id,size = nref)

       #dtextra = dtextra[id %in% idtemp,]

       #genoextra = pullSnpGeno(pop_gaint[pop_gaint@id %in% dtextra$id])

       #genofam = rbind(genofam,genoextra)

       mapdt = makemap()

       genodtput = makeped(z = genofam)

       dtfam = alldt[id %in% rownames(genofam),]
       #dtfam = rbind(dtfam,dtextra)

       coutfam = which(needfam == x)

       fwrite(genodtput,file=paste("hib",coutfam,".ped",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

       fwrite(mapdt, file = paste("hib",coutfam,".map",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

       fwrite(dtfam,file = paste("pheno",coutfam,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

       system(paste("plink --file hib",coutfam," --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --make-xrm --threads 32 --bfile hib",coutfam," --add --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 6 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[2,2],",",(varP(pop)[2,2]-varA(pop)[2,2]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_ST<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_ST)[2] = "hib.GA"

       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 7 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[1,1],",",(varP(pop)[1,1]-varA(pop)[1,1]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_BW<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_BW)[2] = "hib.GA"
       
        sebv_BW = hiebv_BW$hib.GA[match(pop_v@id,hiebv_BW$ID)]

        mebv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]
        
        febv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]

        sebv_ST = hiebv_ST$hib.GA[match(pop_v@id,hiebv_ST$ID)]

        mebv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]
        
        febv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]

        pop_v@ebv = cbind(sebv_BW + ((mebv_BW+febv_BW)/2),sebv_ST + ((mebv_ST+febv_ST)/2))
        
        return(pop_v)

    })


    nref = 50

    cand_list_extra_150 = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,.(id,suvtime)]
        pheno_t = pheno_a[id%in%pop_test@id,]
        
        pop_v = pop_cand[pop_cand@id%in%pheno_a$id]

       genofam= genodt[match(c(pop_v@id,pheno_t$id),rownames(genodt)),]
       #genofam= genodt[match(c(pop_v@id,pheno_t$id,searchParentsFam(pop = pop_v)),rownames(genodt)),]
       
       dtextra = alldt_extra[fam==x,]

       idtemp= sample(dtextra$id,size = nref)

       dtextra = dtextra[id %in% idtemp,]

       genoextra = pullSnpGeno(pop_gaint[pop_gaint@id %in% dtextra$id])

       genofam = rbind(genofam,genoextra)

       mapdt = makemap()

       genodtput = makeped(z = genofam)

       dtfam = alldt[id %in% rownames(genofam),]
       dtfam = rbind(dtfam,dtextra)

       coutfam = which(needfam == x)

       fwrite(genodtput,file=paste("hib",coutfam,".ped",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

       fwrite(mapdt, file = paste("hib",coutfam,".map",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

       fwrite(dtfam,file = paste("pheno",coutfam,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

       system(paste("plink --file hib",coutfam," --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --make-xrm --threads 32 --bfile hib",coutfam," --add --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 6 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[2,2],",",(varP(pop)[2,2]-varA(pop)[2,2]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_ST<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_ST)[2] = "hib.GA"

       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 7 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[1,1],",",(varP(pop)[1,1]-varA(pop)[1,1]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_BW<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_BW)[2] = "hib.GA"
       
        sebv_BW = hiebv_BW$hib.GA[match(pop_v@id,hiebv_BW$ID)]

        mebv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]
        
        febv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]

        sebv_ST = hiebv_ST$hib.GA[match(pop_v@id,hiebv_ST$ID)]

        mebv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]
        
        febv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]

        pop_v@ebv = cbind(sebv_BW + ((mebv_BW+febv_BW)/2),sebv_ST + ((mebv_ST+febv_ST)/2))
        
        return(pop_v)

    })

    nref = 150

    cand_list_extra_250 = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,.(id,suvtime)]
        pheno_t = pheno_a[id%in%pop_test@id,]
        
        pop_v = pop_cand[pop_cand@id%in%pheno_a$id]

       genofam= genodt[match(c(pop_v@id,pheno_t$id),rownames(genodt)),]
       #genofam= genodt[match(c(pop_v@id,pheno_t$id,searchParentsFam(pop = pop_v)),rownames(genodt)),]
       
       dtextra = alldt_extra[fam==x,]

       idtemp= sample(dtextra$id,size = nref)

       dtextra = dtextra[id %in% idtemp,]

       genoextra = pullSnpGeno(pop_gaint[pop_gaint@id %in% dtextra$id])

       genofam = rbind(genofam,genoextra)

       mapdt = makemap()

       genodtput = makeped(z = genofam)

       dtfam = alldt[id %in% rownames(genofam),]
       dtfam = rbind(dtfam,dtextra)

       coutfam = which(needfam == x)

       fwrite(genodtput,file=paste("hib",coutfam,".ped",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

       fwrite(mapdt, file = paste("hib",coutfam,".map",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

       fwrite(dtfam,file = paste("pheno",coutfam,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

       system(paste("plink --file hib",coutfam," --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --make-xrm --threads 32 --bfile hib",coutfam," --add --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 6 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[2,2],",",(varP(pop)[2,2]-varA(pop)[2,2]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_ST<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_ST)[2] = "hib.GA"

       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 7 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[1,1],",",(varP(pop)[1,1]-varA(pop)[1,1]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_BW<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_BW)[2] = "hib.GA"
       
        sebv_BW = hiebv_BW$hib.GA[match(pop_v@id,hiebv_BW$ID)]

        mebv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]
        
        febv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]

        sebv_ST = hiebv_ST$hib.GA[match(pop_v@id,hiebv_ST$ID)]

        mebv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]
        
        febv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]

        pop_v@ebv = cbind(sebv_BW + ((mebv_BW+febv_BW)/2),sebv_ST + ((mebv_ST+febv_ST)/2))
        
        return(pop_v)

    })


     nref = 400

    cand_list_extra_500 = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,.(id,suvtime)]
        pheno_t = pheno_a[id%in%pop_test@id,]
        
        pop_v = pop_cand[pop_cand@id%in%pheno_a$id]

       genofam= genodt[match(c(pop_v@id,pheno_t$id),rownames(genodt)),]
       #genofam= genodt[match(c(pop_v@id,pheno_t$id,searchParentsFam(pop = pop_v)),rownames(genodt)),]
       
       dtextra = alldt_extra[fam==x,]

       idtemp= sample(dtextra$id,size = nref)

       dtextra = dtextra[id %in% idtemp,]

       genoextra = pullSnpGeno(pop_gaint[pop_gaint@id %in% dtextra$id])

       genofam = rbind(genofam,genoextra)

       mapdt = makemap()

       genodtput = makeped(z = genofam)

       dtfam = alldt[id %in% rownames(genofam),]
       dtfam = rbind(dtfam,dtextra)

       coutfam = which(needfam == x)

       fwrite(genodtput,file=paste("hib",coutfam,".ped",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

       fwrite(mapdt, file = paste("hib",coutfam,".map",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

       fwrite(dtfam,file = paste("pheno",coutfam,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

       system(paste("plink --file hib",coutfam," --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --make-xrm --threads 32 --bfile hib",coutfam," --add --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 6 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[2,2],",",(varP(pop)[2,2]-varA(pop)[2,2]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_ST<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_ST)[2] = "hib.GA"

       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 7 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[1,1],",",(varP(pop)[1,1]-varA(pop)[1,1]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_BW<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_BW)[2] = "hib.GA"
       
        sebv_BW = hiebv_BW$hib.GA[match(pop_v@id,hiebv_BW$ID)]

        mebv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]
        
        febv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]

        sebv_ST = hiebv_ST$hib.GA[match(pop_v@id,hiebv_ST$ID)]

        mebv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]
        
        febv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]

        pop_v@ebv = cbind(sebv_BW + ((mebv_BW+febv_BW)/2),sebv_ST + ((mebv_ST+febv_ST)/2))
        
        return(pop_v)

    })

    #extra ref
    nref = 900

    cand_list_extra_1000 = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,.(id,suvtime)]
        pheno_t = pheno_a[id%in%pop_test@id,]
        
        pop_v = pop_cand[pop_cand@id%in%pheno_a$id]

       genofam= genodt[match(c(pop_v@id,pheno_t$id),rownames(genodt)),]
       #genofam= genodt[match(c(pop_v@id,pheno_t$id,searchParentsFam(pop = pop_v)),rownames(genodt)),]
       
       dtextra = alldt_extra[fam==x,]

       idtemp= sample(dtextra$id,size = nref)

       dtextra = dtextra[id %in% idtemp,]

       genoextra = pullSnpGeno(pop_gaint[pop_gaint@id %in% dtextra$id])

       genofam = rbind(genofam,genoextra)

       mapdt = makemap()

       genodtput = makeped(z = genofam)

       dtfam = alldt[id %in% rownames(genofam),]
       dtfam = rbind(dtfam,dtextra)

       coutfam = which(needfam == x)

       fwrite(genodtput,file=paste("hib",coutfam,".ped",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

       fwrite(mapdt, file = paste("hib",coutfam,".map",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

       fwrite(dtfam,file = paste("pheno",coutfam,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

       system(paste("plink --file hib",coutfam," --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --make-xrm --threads 32 --bfile hib",coutfam," --add --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 6 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[2,2],",",(varP(pop)[2,2]-varA(pop)[2,2]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_ST<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_ST)[2] = "hib.GA"

       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 7 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[1,1],",",(varP(pop)[1,1]-varA(pop)[1,1]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_BW<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_BW)[2] = "hib.GA"
       
        sebv_BW = hiebv_BW$hib.GA[match(pop_v@id,hiebv_BW$ID)]

        mebv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]
        
        febv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]

        sebv_ST = hiebv_ST$hib.GA[match(pop_v@id,hiebv_ST$ID)]

        mebv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]
        
        febv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]

        pop_v@ebv = cbind(sebv_BW + ((mebv_BW+febv_BW)/2),sebv_ST + ((mebv_ST+febv_ST)/2))
        
        return(pop_v)

    })

    nref = 1900

    cand_list_extra_2000 = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,.(id,suvtime)]
        pheno_t = pheno_a[id%in%pop_test@id,]
        
        pop_v = pop_cand[pop_cand@id%in%pheno_a$id]

       genofam= genodt[match(c(pop_v@id,pheno_t$id),rownames(genodt)),]
       #genofam= genodt[match(c(pop_v@id,pheno_t$id,searchParentsFam(pop = pop_v)),rownames(genodt)),]
       
       dtextra = alldt_extra[fam==x,]

       idtemp= sample(dtextra$id,size = nref)

       dtextra = dtextra[id %in% idtemp,]

       genoextra = pullSnpGeno(pop_gaint[pop_gaint@id %in% dtextra$id])

       genofam = rbind(genofam,genoextra)

       mapdt = makemap()

       genodtput = makeped(z = genofam)

       dtfam = alldt[id %in% rownames(genofam),]
       dtfam = rbind(dtfam,dtextra)

       coutfam = which(needfam == x)

       fwrite(genodtput,file=paste("hib",coutfam,".ped",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

       fwrite(mapdt, file = paste("hib",coutfam,".map",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

       fwrite(dtfam,file = paste("pheno",coutfam,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

       system(paste("plink --file hib",coutfam," --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --make-xrm --threads 32 --bfile hib",coutfam," --add --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 6 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[2,2],",",(varP(pop)[2,2]-varA(pop)[2,2]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_ST<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_ST)[2] = "hib.GA"

       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 7 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[1,1],",",(varP(pop)[1,1]-varA(pop)[1,1]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_BW<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_BW)[2] = "hib.GA"
       
        sebv_BW = hiebv_BW$hib.GA[match(pop_v@id,hiebv_BW$ID)]

        mebv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]
        
        febv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]

        sebv_ST = hiebv_ST$hib.GA[match(pop_v@id,hiebv_ST$ID)]

        mebv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]
        
        febv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]

        pop_v@ebv = cbind(sebv_BW + ((mebv_BW+febv_BW)/2),sebv_ST + ((mebv_ST+febv_ST)/2))
        
        return(pop_v)

    })

    nref = 2900

    cand_list_extra_3000 = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,.(id,suvtime)]
        pheno_t = pheno_a[id%in%pop_test@id,]
        
        pop_v = pop_cand[pop_cand@id%in%pheno_a$id]

       genofam= genodt[match(c(pop_v@id,pheno_t$id),rownames(genodt)),]
       #genofam= genodt[match(c(pop_v@id,pheno_t$id,searchParentsFam(pop = pop_v)),rownames(genodt)),]
       
       dtextra = alldt_extra[fam==x,]

       idtemp= sample(dtextra$id,size = nref)

       dtextra = dtextra[id %in% idtemp,]

       genoextra = pullSnpGeno(pop_gaint[pop_gaint@id %in% dtextra$id])

       genofam = rbind(genofam,genoextra)

       mapdt = makemap()

       genodtput = makeped(z = genofam)

       dtfam = alldt[id %in% rownames(genofam),]
       dtfam = rbind(dtfam,dtextra)

       coutfam = which(needfam == x)

       fwrite(genodtput,file=paste("hib",coutfam,".ped",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

       fwrite(mapdt, file = paste("hib",coutfam,".map",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

       fwrite(dtfam,file = paste("pheno",coutfam,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

       system(paste("plink --file hib",coutfam," --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --make-xrm --threads 32 --bfile hib",coutfam," --add --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 6 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[2,2],",",(varP(pop)[2,2]-varA(pop)[2,2]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_ST<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_ST)[2] = "hib.GA"

       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 7 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[1,1],",",(varP(pop)[1,1]-varA(pop)[1,1]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_BW<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_BW)[2] = "hib.GA"
       
        sebv_BW = hiebv_BW$hib.GA[match(pop_v@id,hiebv_BW$ID)]

        mebv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]
        
        febv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]

        sebv_ST = hiebv_ST$hib.GA[match(pop_v@id,hiebv_ST$ID)]

        mebv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]
        
        febv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]

        pop_v@ebv = cbind(sebv_BW + ((mebv_BW+febv_BW)/2),sebv_ST + ((mebv_ST+febv_ST)/2))
        
        return(pop_v)

    })

    nref = 3900

    cand_list_extra_4000 = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,.(id,suvtime)]
        pheno_t = pheno_a[id%in%pop_test@id,]
        
        pop_v = pop_cand[pop_cand@id%in%pheno_a$id]

       genofam= genodt[match(c(pop_v@id,pheno_t$id),rownames(genodt)),]
       #genofam= genodt[match(c(pop_v@id,pheno_t$id,searchParentsFam(pop = pop_v)),rownames(genodt)),]
       
       dtextra = alldt_extra[fam==x,]

       idtemp= sample(dtextra$id,size = nref)

       dtextra = dtextra[id %in% idtemp,]

       genoextra = pullSnpGeno(pop_gaint[pop_gaint@id %in% dtextra$id])

       genofam = rbind(genofam,genoextra)

       mapdt = makemap()

       genodtput = makeped(z = genofam)

       dtfam = alldt[id %in% rownames(genofam),]
       dtfam = rbind(dtfam,dtextra)

       coutfam = which(needfam == x)

       fwrite(genodtput,file=paste("hib",coutfam,".ped",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

       fwrite(mapdt, file = paste("hib",coutfam,".map",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

       fwrite(dtfam,file = paste("pheno",coutfam,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

       system(paste("plink --file hib",coutfam," --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --make-xrm --threads 32 --bfile hib",coutfam," --add --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 6 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[2,2],",",(varP(pop)[2,2]-varA(pop)[2,2]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_ST<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_ST)[2] = "hib.GA"

       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 7 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[1,1],",",(varP(pop)[1,1]-varA(pop)[1,1]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_BW<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_BW)[2] = "hib.GA"
       
        sebv_BW = hiebv_BW$hib.GA[match(pop_v@id,hiebv_BW$ID)]

        mebv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]
        
        febv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]

        sebv_ST = hiebv_ST$hib.GA[match(pop_v@id,hiebv_ST$ID)]

        mebv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]
        
        febv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]

        pop_v@ebv = cbind(sebv_BW + ((mebv_BW+febv_BW)/2),sebv_ST + ((mebv_ST+febv_ST)/2))
        
        return(pop_v)

    })


    nref = 4900

    cand_list_extra_5000 = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,.(id,suvtime)]
        pheno_t = pheno_a[id%in%pop_test@id,]
        
        pop_v = pop_cand[pop_cand@id%in%pheno_a$id]

       genofam= genodt[match(c(pop_v@id,pheno_t$id),rownames(genodt)),]
       #genofam= genodt[match(c(pop_v@id,pheno_t$id,searchParentsFam(pop = pop_v)),rownames(genodt)),]
       
       dtextra = alldt_extra[fam==x,]

       idtemp= sample(dtextra$id,size = nref)

       dtextra = dtextra[id %in% idtemp,]

       genoextra = pullSnpGeno(pop_gaint[pop_gaint@id %in% dtextra$id])

       genofam = rbind(genofam,genoextra)

       mapdt = makemap()

       genodtput = makeped(z = genofam)

       dtfam = alldt[id %in% rownames(genofam),]
       dtfam = rbind(dtfam,dtextra)

       coutfam = which(needfam == x)

       fwrite(genodtput,file=paste("hib",coutfam,".ped",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

       fwrite(mapdt, file = paste("hib",coutfam,".map",sep = ""),col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

       fwrite(dtfam,file = paste("pheno",coutfam,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

       system(paste("plink --file hib",coutfam," --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --make-xrm --threads 32 --bfile hib",coutfam," --add --out hib",coutfam,sep = ""))
 
       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 6 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[2,2],",",(varP(pop)[2,2]-varA(pop)[2,2]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_ST<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_ST)[2] = "hib.GA"

       system(paste("hiblup --mme --pheno pheno",coutfam,".csv --pheno-pos 7 --xrm hib",coutfam,".GA --vc-priors ",varA(pop)[1,1],",",(varP(pop)[1,1]-varA(pop)[1,1]),
               " --pcg --threads 32 --out hib",coutfam,sep = ""))

       hiebv_BW<- fread(paste("hib",coutfam,".rand",sep = ""),sep = "\t")

       colnames(hiebv_BW)[2] = "hib.GA"
       
        sebv_BW = hiebv_BW$hib.GA[match(pop_v@id,hiebv_BW$ID)]

        mebv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]
        
        febv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]

        sebv_ST = hiebv_ST$hib.GA[match(pop_v@id,hiebv_ST$ID)]

        mebv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]
        
        febv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]

        pop_v@ebv = cbind(sebv_BW + ((mebv_BW+febv_BW)/2),sebv_ST + ((mebv_ST+febv_ST)/2))
        
        return(pop_v)

    })

#=======================================================================================================
    if(qc){

        onefam_ac_qc_BW = calacforeach(cand_list = cand_list_qc, pname = "onefamqc",ttype = "BW")

        onefam_ac_qc_ST = calacforeach(cand_list = cand_list_qc, pname = "onefamqc",ttype = "ST")

    }
    

    onefam_ac_BW = calacforeach(cand_list = cand_list, pname = "onefam",ttype = "BW")
    
    fourfam_ac_BW = calacforeach(cand_list = cand_list_four, pname = "fourfam",ttype = "BW")

    relfam_ac_BW = calacforeach(cand_list = cand_list_rel, pname = "relfam",ttype = "BW")

    extra_ac_BW_50 = calacforeach(cand_list = cand_list_extra_50, pname = "extra50",ttype = "BW")
    extra_ac_BW_150 = calacforeach(cand_list = cand_list_extra_150, pname = "extra150",ttype = "BW")
    extra_ac_BW_250 = calacforeach(cand_list = cand_list_extra_250, pname = "extra250",ttype = "BW")
    extra_ac_BW_500 = calacforeach(cand_list = cand_list_extra_500, pname = "extra500",ttype = "BW")
    extra_ac_BW_1000 = calacforeach(cand_list = cand_list_extra_1000, pname = "extra1000",ttype = "BW")
    extra_ac_BW_2000 = calacforeach(cand_list = cand_list_extra_2000, pname = "extra2000",ttype = "BW")
    extra_ac_BW_3000 = calacforeach(cand_list = cand_list_extra_3000, pname = "extra3000",ttype = "BW")
    extra_ac_BW_4000 = calacforeach(cand_list = cand_list_extra_4000, pname = "extra4000",ttype = "BW")
    extra_ac_BW_5000 = calacforeach(cand_list = cand_list_extra_5000, pname = "extra5000",ttype = "BW")

    onefam_ac_ST = calacforeach(cand_list = cand_list, pname = "onefam",ttype = "ST")
    
    fourfam_ac_ST = calacforeach(cand_list = cand_list_four, pname = "fourfam",ttype = "ST")

    relfam_ac_ST = calacforeach(cand_list = cand_list_rel, pname = "relfam",ttype = "ST")

    extra_ac_ST_50 = calacforeach(cand_list = cand_list_extra_50, pname = "extra50",ttype = "ST")
    extra_ac_ST_150 = calacforeach(cand_list = cand_list_extra_150, pname = "extra150",ttype = "ST")
    extra_ac_ST_250 = calacforeach(cand_list = cand_list_extra_250, pname = "extra250",ttype = "ST")
    extra_ac_ST_500 = calacforeach(cand_list = cand_list_extra_500, pname = "extra500",ttype = "ST")
    extra_ac_ST_1000 = calacforeach(cand_list = cand_list_extra_1000, pname = "extra1000",ttype = "ST")
    extra_ac_ST_2000 = calacforeach(cand_list = cand_list_extra_2000, pname = "extra2000",ttype = "ST")
    extra_ac_ST_3000 = calacforeach(cand_list = cand_list_extra_3000, pname = "extra3000",ttype = "ST")
    extra_ac_ST_4000 = calacforeach(cand_list = cand_list_extra_4000, pname = "extra4000",ttype = "ST")
    extra_ac_ST_5000 = calacforeach(cand_list = cand_list_extra_5000, pname = "extra5000",ttype = "ST")

}else{

     cand_list_four = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,.(id,suvtime)]
        pheno_t = pheno_a[id%in%pop_test@id,]
        
        pop_v = pop_cand[pop_cand@id%in%pheno_a$id]

       
       #genofam= genodt[match(c(pop_v@id,pheno_t$id,searchParentsFam(pop = pop_v)),rownames(genodt)),]

        meff = RRBLUP(pop = genodt_t[genodt_t@id %in% c(pheno_t$id,searchParentsFam(pop = pop_v))],traits = 2, snpChip = 1)
        
        pop_v = setEBV(pop = pop_v,meff,simParam = SP)
       
        sebv_ST = pop_v@ebv[,1]

        meff = RRBLUP(pop = genodt_t[genodt_t@id %in% c(pheno_t$id,searchParentsFam(pop = pop_v))],traits = 1, snpChip = 1)
        
        pop_v = setEBV(pop = pop_v,meff,simParam = SP)

        sebv_BW = pop_v@ebv[,1]

        mebv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]
        
        febv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]

        mebv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]
        
        febv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]

        pop_v@ebv = cbind(sebv_BW + ((mebv_BW+febv_BW)/2),sebv_ST + ((mebv_ST+febv_ST)/2))
        
        return(pop_v)

    })

     cand_list = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,.(id,suvtime)]
        pheno_t = pheno_a[id%in%pop_test@id,]
        
        pop_v = pop_cand[pop_cand@id%in%pheno_a$id]


        meff = RRBLUP(pop = genodt_t[genodt_t@id %in%pheno_t$id],traits = 1,snpChip = 1)
        
        pop_v = setEBV(pop = pop_v,meff,simParam = SP)
       
        sebv_BW = pop_v@ebv[,1]

        meff = RRBLUP(pop = genodt_t[genodt_t@id %in%pheno_t$id],traits = 2,snpChip = 1)
        
        pop_v = setEBV(pop = pop_v,meff,simParam = SP)
       
        sebv_ST = pop_v@ebv[,1]


        mebv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]
        
        febv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]

        mebv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]
        
        febv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]

        pop_v@ebv = cbind(sebv_BW + ((mebv_BW+febv_BW)/2),sebv_ST + ((mebv_ST+febv_ST)/2))
        
        return(pop_v)

    })

    cand_list_rel = lapply(needfam,function(x){
        pheno_a = alldt[fam==x,.(id,suvtime)]
        pheno_t = pheno_a[id%in%pop_test@id,]
        
        pop_v = pop_cand[pop_cand@id%in%pheno_a$id]

       
       #genofam= genodt[match(c(pop_v@id,pheno_t$id,searchParentsFam(pop = pop_v)),rownames(genodt)),]

        meff = RRBLUP(pop = genodt_t[genodt_t@id %in% c(pheno_t$id,searchrealationsFam(pop = pop_v))],traits = 2, snpChip = 1)
        
        pop_v = setEBV(pop = pop_v,meff,simParam = SP)
       
        sebv_ST = pop_v@ebv[,1]

        meff = RRBLUP(pop = genodt_t[genodt_t@id %in% c(pheno_t$id,searchrealationsFam(pop = pop_v))],traits = 1, snpChip = 1)
        
        pop_v = setEBV(pop = pop_v,meff,simParam = SP)

        sebv_BW = pop_v@ebv[,1]

        mebv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]
        
        febv_BW = alldt_pblup$ebv_BW[match(pop_v@mother,alldt_pblup$id)]

        mebv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]
        
        febv_ST = alldt_pblup$ebv_ST[match(pop_v@mother,alldt_pblup$id)]

        pop_v@ebv = cbind(sebv_BW + ((mebv_BW+febv_BW)/2),sebv_ST + ((mebv_ST+febv_ST)/2))
        
        return(pop_v)

    })

    onefam_ac = calacforeach(cand_list = cand_list, pname = "onefam")

    fourfam_ac = calacforeach(cand_list = cand_list_four, pname = "fourfam")
 
}
#====================================================================================================
if(qc){

    finalac= rbind(onepop_ac_BW,fourpop_ac_BW,onefam_ac_BW,onefam_ac_qc_BW,fourfam_ac_BW,relfam_ac_BW,extra_ac_BW_50,extra_ac_BW_150,extra_ac_BW_250,extra_ac_BW_500,extra_ac_BW_1000,extra_ac_BW_2000,extra_ac_BW_3000,extra_ac_BW_4000,extra_ac_BW_5000,
                   onepop_ac_ST,fourpop_ac_ST,onefam_ac_ST,onefam_ac_qc_ST,fourfam_ac_ST,relfam_ac_ST,extra_ac_ST_50,extra_ac_ST_150,extra_ac_ST_250,extra_ac_ST_500,extra_ac_ST_1000,extra_ac_ST_2000,extra_ac_ST_3000,extra_ac_ST_4000,extra_ac_ST_5000) 

}else{
    finalac= rbind(onepop_ac_BW,fourpop_ac_BW,onefam_ac_BW,fourfam_ac_BW,relfam_ac_BW,extra_ac_BW_50,extra_ac_BW_150,extra_ac_BW_250,extra_ac_BW_500,extra_ac_BW_1000,extra_ac_BW_2000,extra_ac_BW_3000,extra_ac_BW_4000,extra_ac_BW_5000,
                   onepop_ac_ST,fourpop_ac_ST,onefam_ac_ST,fourfam_ac_ST,relfam_ac_ST,extra_ac_ST_50,extra_ac_ST_150,extra_ac_ST_250,extra_ac_ST_500,extra_ac_ST_1000,extra_ac_ST_2000,extra_ac_ST_3000,extra_ac_ST_4000,extra_ac_ST_5000) 
}


#============================================calidenticalsnp============================================

 famsnpchek = calidenticalfamsnp(needfam)

#==============================================================================================================




#==============================output result===============================================

fwrite(finalac,paste("acck",r,".csv",sep = ""),sep = ",")


fwrite(famsnpchek,paste("snpck",r,".csv",sep = ""),sep = ",")

#==========================================================================================

