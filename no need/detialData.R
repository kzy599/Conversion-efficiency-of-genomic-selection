setwd("/home/GTDisk1/kangziyi/wssvcode/output/")
zt1 = fread("allzt.csv",sep = ",")
output = fread("alldt.csv",sep = ",")
output$GenerationInterval <- cut(output$nGeneration,
                                 breaks = seq(0, 20, by = 5),
                                 include.lowest = TRUE,
                                 labels = c("0-5", "6-10", "11-15", "16-20"))
output[denSNP >= 3 & denSNP <= 46,SnpInterval:="3-46"]
output[114<=denSNP,SnpInterval:="114-1250"]
output[nRef>=10&nRef<=70,RefInterval:="10-70"]
output[nRef>=100&nRef<=200,RefInterval:="100-200"]
unique(output$RefInterval)
unique(output$SnpInterval)
unique(output$GenerationInterval)

#only GS
outGS = output[denSNP!=0,]
outGS$SnpInterval = factor(outGS$SnpInterval,levels = c("3-46","114-1250"))
outGS$RefInterval = factor(outGS$RefInterval,levels = c("10-70","100-200"))
outGS$SnpInterval %>% unique()
outGS$RefInterval %>% unique()
outGS[nGeneration<=10,GenerationInterval:="1-10"]
outGS[nGeneration>=11,GenerationInterval:="11-20"]

#genetic gain
system("python3 3dplot.py --ITYPE F --VTYPE gg0 --VSE gg0se --plotname 3dplotFgeneticgain.pdf --ZTYPE 'Genetic gain' ")
system("python3 3dplot.py --ITYPE T --VTYPE gg0 --VSE gg0se --plotname 3dplotTgeneticgain.pdf --ZTYPE 'Genetic gain' ")
system("python3 maketable.py --VTYPE gg0 --VSE gg0se --TNAME geneticgaintable.docx")

#genetic diversity
system("python3 3dplot.py --ITYPE F --VTYPE genicVa --VSE genicVase --plotname 3dplotFgenicVa.pdf --ZTYPE 'Genetic diversity' ")
system("python3 3dplot.py --ITYPE T --VTYPE genicVa --VSE genicVase --plotname 3dplotTgenicVa.pdf --ZTYPE 'Genetic diversity' ")
system("python3 maketable.py --VTYPE genicVa --VSE genicVase --TNAME genicVatable.docx")

#effective population size
system("python3 3dplot.py --ITYPE F --VTYPE genicNe --VSE genicNese --plotname 3dplotFgenicNe.pdf --ZTYPE 'Effective population size' ")
system("python3 3dplot.py --ITYPE T --VTYPE genicNe --VSE genicNese --plotname 3dplotTgenicNe.pdf --ZTYPE 'Effective population size' ")
system("python3 maketable.py --VTYPE genicNe --VSE genicNese --TNAME genicNetable.docx")
system("python3 maketable.py --VTYPE Ne --VSE Nese --TNAME Netable.docx")
#conversion efficiency
system("python3 3dplot.py --ITYPE F --VTYPE genicEff --VSE genicEffse --plotname 3dplotFgenicEff.pdf --ZTYPE 'Conversion efficiency' ")
system("python3 3dplot.py --ITYPE T --VTYPE genicEff --VSE genicEffse --plotname 3dplotTgenicEff.pdf --ZTYPE 'Conversion efficiency' ")
system("python3 maketable.py --VTYPE genicEff --VSE genicEffse --TNAME genicEfftable.docx")

im = "F"
anova_model<- aov(genicVa~factor(nRef)+factor(denSNP)+factor(denSNP)*factor(nRef),data =outGS[nGeneration==20&imp==im,] )
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"factor(denSNP)",conf.level = 0.95)
tukey_results <- TukeyHSD(anova_model,"factor(nRef)",conf.level = 0.95)
tukey_results

allmodel = lm(log(genicVa)~denSNP+nRef+denSNP*nRef,data = outGS[nGeneration==20 & imp == "F"&SnpInterval=="114-1250",])
allmodel = lm(log(genicVa)~denSNP+nRef+denSNP*nRef,data = outGS[nGeneration==20 & imp == "F"&SnpInterval=="3-46",])
summary(allmodel)
(exp(allmodel$coefficients[3])-1)*100*10

g = 20
nR = 100
im = "F"
denS = 1250
chkt = t.test(x = output[nGeneration==g & imp == im & denSNP == 23,genicVa],y=output[nGeneration==g & imp == im & denSNP ==46,genicVa])
chkt$p.value
lh = zt1[nGeneration==g&imp==im&denSNP==114,gg0] %>% round(2) %>% mean()
rh = zt1[nGeneration==g&imp==im&denSNP==46,gg0] %>% round(2) %>% mean()
(((lh - rh)/rh)*100) %>% round(2)

#contrast with ps
#genetic gain
ps_gg = zt1[nGeneration==g&denSNP==0,gg0] %>% round(2) %>% mean()
percentg = c()
im = "F"
for(den in c(3,12,23,46,114,228,455,682,1250)){
  if(den == 1250) {
    temp_gg = zt1[nGeneration==g&imp=="F"&denSNP==den,gg0] %>% round(2) %>% mean()
  }else{
    temp_gg = zt1[nGeneration==g&imp==im&denSNP==den,gg0] %>% round(2) %>% mean()
  }
  
  percentg =c(percentg, (((temp_gg - ps_gg)/ps_gg)*100) %>% round(2))
}
percentg
percentg[-1] %>% max()
percentg[-1] %>% min()

#genetic diversity
ps_gg = zt1[nGeneration==g&denSNP==0,genicVa] %>% round(2) %>% mean()
percentg = c()
im = "F"
for(den in c(3,12,23,46,114,228,455,682,1250)){
  if(den == 1250) {
    temp_gg = zt1[nGeneration==g&imp=="F"&denSNP==den,genicVa] %>% round(2) %>% mean()
  }else{
    temp_gg = zt1[nGeneration==g&imp==im&denSNP==den,genicVa] %>% round(2) %>% mean()
    }
  
  percentg =c(percentg, (((temp_gg - ps_gg)/ps_gg)*100) %>% round(2))
}
percentg
percentg %>% max()
percentg %>% min()

#effective population size

ps_gg = zt1[nGeneration==g&denSNP==0,genicNe] %>% round(0) %>% mean() %>% round(0)
percentg = c()
im = "F"
for(den in c(3,12,23,46,114,228,455,682,1250)){
  if(den == 1250) {
    temp_gg = zt1[nGeneration==g&imp=="F"&denSNP==den,genicNe] %>% round(0) %>% mean()
  }else{
    temp_gg = zt1[nGeneration==g&imp==im&denSNP==den,genicNe] %>% round(0) %>% mean()
  }
  
  percentg =c(percentg, temp_gg)
}

ps_gg = zt1[nGeneration==g&denSNP==0,Ne] %>% round(0) %>% mean() %>% round(0)
percentg = c()
im = "F"
for(den in c(3,12,23,46,114,228,455,682,1250)){
  if(den == 1250) {
    temp_gg = zt1[nGeneration==g&imp=="F"&denSNP==den,Ne] %>% round(0) %>% mean()
  }else{
    temp_gg = zt1[nGeneration==g&imp==im&denSNP==den,Ne] %>% round(0) %>% mean()
  }
  
  percentg =c(percentg, temp_gg)
}
percentg %>% round(0) %>% mean()
(percentg %>% round(0)) - ps_gg
percentg %>% max()
percentg %>% min()

#conversion efficiency
ps_gg = zt1[nGeneration==g&denSNP==0,genicEff] %>% round(2) %>% mean()
percentg = c()
im = "T"
for(den in c(3,12,23,46,114,228,455,682,1250)){
  if(den == 1250) {
    temp_gg = zt1[nGeneration==g&imp=="F"&denSNP==den,genicEff] %>% round(2) %>% mean()
  }else{
    temp_gg = zt1[nGeneration==g&imp==im&denSNP==den,genicEff] %>% round(2) %>% mean()
  }
  
  percentg =c(percentg, (((temp_gg - ps_gg)/ps_gg)*100) %>% round(2))
}
percentg
percentg[-1] %>% max()
percentg[-1] %>% min()

den = 23
t.test(x = output[nGeneration==g & imp == "T" & denSNP == den,genicEff],y=output[nGeneration==g & denSNP ==0,genicEff])


#contrast with imputation

#genetic gain
percentg = c()
for(den in c(3,12,23,46,114)){
  temp_gg = zt1[nGeneration==g&imp=="T"&denSNP==den,gg0] %>% round(2) %>% mean()
  ps_gg = zt1[nGeneration==g&imp=="F"&denSNP==den,gg0] %>% round(2) %>% mean()
  percentg =c(percentg, (((temp_gg - ps_gg)/ps_gg)*100) %>% round(2))
}
percentg
percentg %>% max()
percentg %>% min()

den = c(3,12,23,46,114,228,455,682)
ref = c(10,30,50,70,100,150,200)
anova_model<- aov(gg0~imp,data =output[nGeneration==g & nRef %in% ref & denSNP %in% den,])
summary(anova_model)

den = 3
t.test(x = output[nGeneration==g & imp == "F" & denSNP == den,gg0],y=output[nGeneration==g & imp == "T" & denSNP ==den,gg0])

#genetic diversity
percentg = c()
for(den in c(3,12,23,46,114)){
  temp_gg = zt1[nGeneration==g&imp=="T"&denSNP==den,genicVa] %>% round(2) %>% mean()
  ps_gg = zt1[nGeneration==g&imp=="F"&denSNP==den,genicVa] %>% round(2) %>% mean()
  percentg =c(percentg, (((temp_gg - ps_gg)/ps_gg)*100) %>% round(2))
}
percentg
percentg %>% max()
percentg %>% min()

den = c(3,12,23,46,114,228,455,682)
ref = c(10,30,50,70,100,150,200)
anova_model<- aov(genicVa~imp,data =output[nGeneration==g & nRef %in% ref & denSNP %in% den,])
summary(anova_model)

den = 23
t.test(x = output[nGeneration==g & imp == "F" & denSNP == den,genicVa],y=output[nGeneration==g & imp == "T" & denSNP ==den,genicVa])

#effective populatuon size
percentg = c()
percentgF = c()
for(den in c(3,12,23,46,114)){
  temp_gg = zt1[nGeneration==g&imp=="T"&denSNP==den,genicNe] %>% round(0) %>% mean()
  ps_gg = zt1[nGeneration==g&imp=="F"&denSNP==den,genicNe] %>% round(0) %>% mean()
  percentg =c(percentg, temp_gg)
  percentgF = c(percentgF,ps_gg)
}
percentg
percentg %>% max()
percentg %>% min()

percentgF
percentgF %>% max()
percentgF %>% min()

den = c(3,12,23,46,114,228,455,682)
ref = c(10,30,50,70,100,150,200)
anova_model<- aov(genicNe~imp,data =output[nGeneration==g & nRef %in% ref & denSNP %in% den,])
summary(anova_model)

den = 12
t.test(x = output[nGeneration==g & imp == "F" & denSNP == den,genicNe],y=output[nGeneration==g & imp == "T" & denSNP ==den,genicNe])

#conversion efficiency
percentg = c()
dens = c(3,12,23,46,114)
dens = c(3,12,23)
dens = c(46,114,228)
g = 20
for(den in dens){
  temp_gg = zt1[nGeneration==g&imp=="T"&denSNP==den,genicEff] %>% round(2) %>% mean()
  ps_gg = zt1[nGeneration==g&imp=="F"&denSNP==den,genicEff] %>% round(2) %>% mean()
  percentg =c(percentg, (((temp_gg - ps_gg)/ps_gg)*100) %>% round(2))
}
percentg

percentg %>% max()

percentg %>% min()

percentg %>% mean()

den = c(46,114,228,455,682)
den = c(3,12,23)
ref = c(10,30,50,70,100,150,200)
anova_model<- aov(genicEff~imp,data =output[nGeneration==g & nRef %in% ref & denSNP %in% den,])
summary(anova_model)

den = 46
t.test(x = output[nGeneration==g & imp == "F" & denSNP == den,genicEff],y=output[nGeneration==g & imp == "T" & denSNP ==den,genicEff])

t.test(x = output[nGeneration==g & imp == "F" & denSNP == 0,genicEff],y=output[nGeneration==g & imp == "T" & denSNP ==46,genicEff])

#Impact of panel density 

#genetic gain
den = c(3,12,23,46,114)
den = c(114,228,455,682)
den = c(3,12,23,46,114,228,455,682)
modelrate = lm(log(gg0)~denSNP*nRef,data = output[nGeneration==g & nRef %in% c(10,30,50,70,100,150,200) & imp == "F" & denSNP %in% den,])
(((exp(modelrate$coefficients[2])-1))*100*10) %>% round(2)
summary(modelrate)

den = c(3,12,23,46,114,228,455,682,1250)
den = c(114,228,455,682,1250)
anova_model<- aov(gg0~factor(denSNP),data =output[nGeneration==g & nRef %in% c(10,30,50,70,100,150,200) & imp=="F" &  denSNP %in% den,])
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"factor(denSNP)",conf.level = 0.95)
tukey_results

#genetic diversity
den = c(3,12,23,46,114)
den = c(114,228,455,682)
den = c(3,12,23,46,114,228,455,682)
den = c(3,12,23,46,114,228,455,682,1250)
modelrate = lm(log(genicVa)~denSNP+nRef,data = output[nGeneration==g & nRef %in% c(10,30,50,70,100,150,200) & imp == "F" & denSNP %in% den,])
(((exp(modelrate$coefficients[2])-1))*100*10) %>% round(2)
summary(modelrate)

den = c(3,12,23,46,114,228,455,682,1250)
anova_model<- aov(genicVa~factor(denSNP),data =output[nGeneration==g & nRef %in% c(10,30,50,70,100,150,200) & imp=="F" &  denSNP %in% den,])
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"factor(denSNP)",conf.level = 0.95)
tukey_results

#effective population size


den = c(3,12,23,46,114,228,455,682,1250)
den = c(114,228,455,682,1250)
den = c(3,12,23,46,114)
anova_model<- aov(genicNe~factor(denSNP),data =output[nGeneration==g & nRef %in% c(10,30,50,70,100,150,200) & imp=="F" &  denSNP %in% den,])
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"factor(denSNP)",conf.level = 0.95)
tukey_results


#conversion efficiency
den = c(3,12,23,46)
den = c(46,114,228,455,682,1250)
den = c(3,12,23,46,114,228,455,682,1250)
modelrate = lm(log(genicEff)~denSNP*nRef,data = output[nGeneration==g & nRef %in% c(10,30,50,70,100,150,200) & imp == "T" & denSNP %in% den,])
(((exp(modelrate$coefficients[2])-1))*100*10) %>% round(2)
summary(modelrate)

den = c(3,12,23,46)
den = c(46,114,228,455,682,1250)
anova_model<- aov(genicEff~factor(denSNP),data =output[nGeneration==g & nRef %in% c(10,30,50,70,100,150,200) & imp=="F" &  denSNP %in% den,])
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"factor(denSNP)",conf.level = 0.95)
tukey_results

#within-family accuracy
#average on all generations and reference group size
den = c(12,23,46,114,228,455,682,1250)
famac_allgen = zt1[denSNP%in%den&imp == "F",.(denSNP,nGeneration,acforfam)]
famac_allgen[,lapply(.SD, mean), by = denSNP]
famac_allgen[,lapply(.SD, sd), by = denSNP]


famac_allgen[, lapply(.SD, mean), by = c("denSNP", "nGeneration")]
famac_allgen[, lapply(.SD, sd), by = c("denSNP", "nGeneration")]

#Impact of reference size
#genetic gain
den = c(3,12,23,46,114,228,455,682)
den = c(3,12,23,46,114,228,455,682,1250)
ref = c(10,30,50,70,100,150,200)
ref = c(10,30,50,70)
ref = c(70,100,150,200)
modelrate = lm(log(gg0)~denSNP*nRef,data = output[nGeneration==g & nRef %in% ref & imp=="F" &  denSNP %in% den,])
(((exp(modelrate$coefficients[3])-1))*100*10) %>% round(2)
summary(modelrate)
modelrate = lm(log(gg0)~nRef,data = output[nGeneration==g & nRef %in% ref & denSNP==1250,])
(((exp(modelrate$coefficients[2])-1))*100*10) %>% round(2)
summary(modelrate)
modelrate = lm(log(gg0)~nRef,data = output[nGeneration==g & nRef %in% ref & denSNP==0,])
(((exp(modelrate$coefficients[2])-1))*100*10) %>% round(2)
summary(modelrate)

anova_model<- aov(gg0~factor(nRef),data =output[nGeneration==g & nRef %in% ref & imp=="F" &  denSNP %in% den,])
summary(anova_model)
anova_model<- aov(gg0~factor(nRef),data =output[nGeneration==g & nRef %in% ref & imp=="T" &  denSNP %in% den,])
summary(anova_model)
anova_model<- aov(gg0~factor(nRef),data =output[nGeneration==g & nRef %in% ref & denSNP == 1250,])
summary(anova_model)
anova_model<- aov(gg0~factor(nRef),data =output[nGeneration==g & nRef %in% ref & denSNP== 0,])
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"factor(nRef)",conf.level = 0.95)
tukey_results


#genetic diversity
den = c(3,12,23,46,114,228,455,682)
den = c(3,12,23,46,114)
den = c(114,228,455,682,1250)
den = c(3,12,23,46,114,228,455,682,1250)
den = c(3,12)
den = c(23,46,114,228,455,682,1250)
ref = c(10,30,50,70,100,150,200)
ref = c(10,30,50,70)
ref = c(70,100,150,200)
modelrate = lm(log(genicVa)~SnpInterval*nRef,data = output[nGeneration==g & nRef %in% ref & imp=="F" &  denSNP %in% den,])
(((exp(modelrate$coefficients[3]+modelrate$coefficients[4])-1))*100*10) %>% round(2)
(((exp(modelrate$coefficients[3])-1))*100*10) %>% round(2)
summary(modelrate)

modelrate = lm(log(genicVa)~nRef+denSNP,data = output[nGeneration==g & nRef %in% ref & imp=="F" &  denSNP %in% den,])
(((exp(modelrate$coefficients[2])-1))*100*10) %>% round(2)
summary(modelrate)

modelrate = lm(log(genicVa)~nRef,data = output[nGeneration==g & nRef %in% ref & denSNP==1250,])
(((exp(modelrate$coefficients[2])-1))*100*10) %>% round(2)
summary(modelrate)
modelrate = lm(log(genicVa)~nRef,data = output[nGeneration==g & nRef %in% ref & denSNP==0,])
(((exp(modelrate$coefficients[2])-1))*100*10) %>% round(2)
summary(modelrate)

anova_model<- aov(genicVa~factor(nRef)*imp,data =output[nGeneration==g & nRef %in% ref &  denSNP %in% den,])
summary(anova_model)

anova_model<- aov(genicVa~factor(nRef),data =output[nGeneration==g & nRef %in% ref & imp=="F" &  denSNP %in% den,])
summary(anova_model)

anova_model<- aov(genicVa~factor(nRef),data =output[nGeneration==g & nRef %in% ref & denSNP == 1250,])
summary(anova_model)
anova_model<- aov(genicVa~factor(nRef),data =output[nGeneration==g & nRef %in% ref & denSNP== 0,])
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"factor(nRef)",conf.level = 0.95)
tukey_results


#effective population size

den = c(3,12)
den = c(23,46,114,228,455,682,1250)
ref = c(10,30,50,70,100,150,200)
anova_model<- aov(genicNe~factor(nRef),data =output[nGeneration==g & nRef %in% ref & imp=="F" &  denSNP %in% den,])
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"factor(nRef)",conf.level = 0.95)
tukey_results

percentg = data.table()
ref = zt1$nRef %>% unique()
for(den in c(23,46,114,228,455,682,1250)){
  ps_gg = zt1[nGeneration==g&imp=="F"&denSNP==den,genicNe] %>% round(0)
  temp_gg= data.table(nRef = ref, genicNe = ps_gg)
  percentg =rbind(percentg, temp_gg)
}
percentgF = c()
for (r in ref) {
  temp_gg = percentg[nRef == r,genicNe] %>% mean()
  percentgF = c(percentgF,temp_gg)
}
percentg = data.table(nRef = ref, genicNe  = round(percentgF))
percentg

ck = zt1[nGeneration==g&imp=="F"&denSNP%in%den,.(nRef,denSNP,genicNe,Ne)]

ck[,lapply(.SD, mean), by = denSNP]
ck[,lapply(.SD, sd), by = denSNP]

ck[,lapply(.SD, mean), by = nRef]
ck[,lapply(.SD, sd), by = nRef]

ck[, lapply(.SD, mean), by = c("denSNP", "nRef")]
ck[, lapply(.SD, sd), by = c("denSNP", "nRef")]

#coversion efficiency
den = c(3,12,23,46)
den = c(46,114,228,455,682,1250)
den = c(3,12,23,46,114,228,455,682,1250)
ref = c(10,30,50,70,100,150,200)
ref = c(30,50,70,100,150,200)
modelrate = lm(log(genicEff)~denSNP*nRef,data = output[nGeneration==g & nRef %in% ref & imp=="F" &  denSNP %in% den,])
(((exp(modelrate$coefficients[3])-1))*100*10) %>% round(2)
summary(modelrate)
modelrate = lm(log(genicEff)~nRef,data = output[nGeneration==g & nRef %in% ref & denSNP==1250,])
(((exp(modelrate$coefficients[2])-1))*100*10) %>% round(2)
summary(modelrate)
modelrate = lm(log(genicEff)~nRef,data = output[nGeneration==g & nRef %in% ref & denSNP==0,])
(((exp(modelrate$coefficients[2])-1))*100*10) %>% round(2)
summary(modelrate)

anova_model<- aov(genicEff~factor(nRef),data =output[nGeneration==g & nRef %in% ref & imp=="F" &  denSNP %in% den,])
summary(anova_model)
anova_model<- aov(genicEff~factor(nRef),data =output[nGeneration==g & nRef %in% ref & imp=="T" &  denSNP %in% den,])
summary(anova_model)
anova_model<- aov(genicEff~factor(nRef),data =output[nGeneration==g & nRef %in% ref & denSNP== 0,])
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"factor(nRef)",conf.level = 0.95)
tukey_results

percentg = data.table()
ref = zt1$nRef %>% unique()
for(den in c(3,12,23,46,114,228,455,682,1250)){
  ps_gg = zt1[nGeneration==g&imp=="F"&denSNP==den,genicEff]
  temp_gg= data.table(nRef = ref, genicEff = ps_gg)
  percentg =rbind(percentg, temp_gg)
}
percentgF = c()
for (r in ref) {
  temp_gg = percentg[nRef == r,genicEff]  %>% round(2) %>% mean()
  percentgF = c(percentgF,temp_gg)
}
percentg = data.table(nRef = ref, genicEff  = round(percentgF,2))
percentg = (percentg$genicEff - percentg[nRef==10,genicEff])/percentg[nRef==10,genicEff]
(percentg*100) %>% round(2)
t.test(x = output[nGeneration==g & nRef == 70 & imp == im & denSNP %in% denS,genicNe],y=output[nGeneration==g & nRef == 50 & imp == im & denSNP %in% denS,genicNe])
t.test(x = output[nGeneration==g & nRef == 70 & imp == im & denSNP %in% denS,genicNe],y=output[nGeneration==g & nRef == 50 & imp == im & denSNP %in% denS,genicNe])
backward_model <- step(allmodel, direction = "backward")
summary(backward_model)


den = c(3,12,23,46,114,228,455,682,1250)
ck = zt1[nGeneration==g&imp=="F"&denSNP%in%den,.(nRef,denSNP,genicEff)]

ck[,lapply(.SD, mean), by = denSNP]
ck[,lapply(.SD, sd), by = denSNP]

ck[,lapply(.SD, mean), by = nRef]
x= ck[,lapply(.SD, mean), by = nRef]$genicEff %>% round(2)
y = x[1]
(((x-y)/y)*100) %>% round(2)
ck[,lapply(.SD, sd), by = nRef]

ck[, lapply(.SD, mean), by = c("denSNP", "nRef")]
ck[, lapply(.SD, sd), by = c("denSNP", "nRef")]


#Additional file
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.line = element_line(colour = "black"),
          #axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.title.x=element_text(face = "bold",size = 14),
          axis.text = element_text(face = "bold",size = 12),#坐标轴刻度标签加粗
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          #legend.title=element_blank(),
          #legend.position=c(0.5, 0.95),#图例在绘图区域的位置
          #legend.position="none",
          legend.position="right",
          #legend.direction = "horizontal",
          legend.direction = "vertical",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=8)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
}
pv_c = c("gg0","genicVa","genicNe","genicEff")
pvse_c = c("gg0se","genicVase","genicNese","genicEffse")
yb_c = c("Genetic gain","Genetic diversity","Effective population size","Conversion efficiency")
fn_c = c("FigureSgg01.pdf","FigureSgenicVa1.pdf","FigureSgenicNe1.pdf","FigureSgenicEff1.pdf")
fn2_c = c("FigureSgg02.pdf","FigureSgenicVa2.pdf","FigureSgenicNe2.pdf","FigureSgenicEff2.pdf")
#Additional file 2, Figures S1 to S3
for(i in 1:length(pv_c)){
  pv = pv_c[i]
  pvse = pvse_c[i]
  yb = yb_c[i]
  fn = fn_c[i]
  fn2 = fn2_c[i]
  zt = zt1[nGeneration==20,]
  zt$denSNP = as.character(zt$denSNP)
  # zt$denSNP = paste(zt$denSNP," ","snp/chr",sep = "")
  levls = c("0","3","12","23","46","114","228","455","682","1250")
  # levls = paste(levls," ","snp/chr",sep = "")
  zt$denSNP = factor(zt$denSNP,levels = levls)
  zt[denSNP=="0",imp:="Ped"]
  zt[denSNP=="1250",imp:="HD"]
  zt[imp == "F",imp:="LD"]
  zt[imp == "T",imp:="Imputation"]
  zt$imp = factor(zt$imp,levels = c("Ped","LD","Imputation","HD"))
  zt$value = zt[,..pv]
  zt$valuese = zt[,..pvse]
  P<- ggplot(data = zt,aes(x=nRef,y=value,group=imp,color = imp))+
    geom_point()+
    geom_line()+
    #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
    xlab("Reference group size (per family)")+
    ylab(yb)+
    theme_zg()+theme(legend.position=c(0.85, 0.12),legend.title=element_blank())+
    #geom_ribbon(aes(ymin = pgvr -pgvrse, ymax = pgvr+pgvrse), alpha = 0.4)+
    #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
    geom_errorbar(aes(ymin=value-valuese,
                      ymax=value+valuese),
                  width=0.05,alpha = 0.5)+
    scale_x_continuous(limits = c(10,200),breaks = c(10,30,50,70,100,150,200))+labs(color="imputation")+scale_color_aaas()+facet_wrap(~ denSNP)
  # P+geom_hline(yintercept = zt[denSNP=="1250"&nRef==200,gg0])
  P
  ggsave(fn, P, width = 20, height = 10, dpi = 300)
  
  
  zt = zt1[nGeneration==20,]
  zt$nRef = as.character(zt$nRef)
  zt$nRef = factor(zt$nRef,levels = c("10","30","50","70","100","150","200"))
  zt[denSNP=="0",imp:="Ped"]
  zt[denSNP=="1250",imp:="HD"]
  zt[imp == "F",imp:="LD"]
  zt[imp == "T",imp:="Imputation"]
  zt$imp = factor(zt$imp,levels = c("Ped","LD","Imputation","HD"))
  zt$snpinterval = ifelse(zt$denSNP<=23,"3-46","46-1250")
  # zt$denSNP = as.character(zt$denSNP)
  # zt$denSNP = factor(zt$denSNP,levels = c("3","12","23","46","114","228","455","682","1250"))
  # pdf("Figure2.pdf", width = 20, height = 10)
  # showtext_auto(enable = TRUE)
  zt$value = zt[,..pv]
  zt$valuese = zt[,..pvse]
  P = ggplot(data = zt[imp%in%c("LD","HD","Ped","Imputation"),],aes(x=denSNP,y=value,group=imp,color = imp))+
    geom_point()+
    geom_line()+
    #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
    xlab("Panel density (per chromosome)")+
    ylab(yb)+
    theme_bw()+
    theme(panel.grid.major = element_line(colour = NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank(),
          text = element_text(family = "STXihei"),
          legend.position = "right",
          legend.background = element_rect(colour = "black"))+
    theme_zg()+theme(legend.position=c(0.85, 0.12),legend.title=element_blank())+
    #geom_ribbon(aes(ymin = pgvr -pgvrse, ymax = pgvr+pgvrse), alpha = 0.4)+
    #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
    geom_errorbar(aes(ymin=value-valuese,
                      ymax=value+valuese),
                  width=0.05,alpha = 0.5)+
    # scale_x_discrete()+
    scale_x_continuous(limits = c(0,1250),breaks = c(0,46,114,228,455,682,1250))+
    labs(color="Imputation")+scale_color_aaas()+facet_wrap(~ nRef)
  # dev.off()
  P
  # P+geom_hline(yintercept = zt[denSNP=="0",gg0])
  ggsave(fn2, P, width = 20, height = 10, dpi = 300)
  
}
