#============================================================

#detialData.R包含了文章中的数据
#============================================================
options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(ggplot2)
library(grid)
library(ggsci)
library(data.table)
library(stringr)
library(purrr)
library(dplyr)
library(devtools)
library(forestplot)
install.packages("kableExtra")
devtools::install_github("tidyverse/ggplot2")
setwd("/home/GTDisk1/kangziyi/wssvcode/output/")
rm(list = ls())
gc()
allnum = list.files()
allnum = allnum[allnum%flike%".csv"]
allnum = allnum[allnum%flike%"Pop"]
# allnum = allnum[!allnum%in%c("alldt.csv")]
allnum = str_extract_all(allnum, "\\d+")
allnum = sapply(c(1:length(allnum)), function(x){
  x = allnum[[x]][1]
  return(x)
})
allnum = as.numeric(unique(allnum))
allnum = allnum[order(allnum)]
pname = c()
for(gs in c("Pop")){
  for(ns in allnum){
    for(im in c("IF","IT")){
      for(rn in c(10,30,50,70,100,150,200)){
        if(ns%in%c(0,12500)&im=="IT") next
        pname = c(pname,paste(gs,ns,im,rn,sep = ""))
      }
    }
  }
}

#===========================================================================================

###
# Enfile = c()
# for(i in 1:10){
#   Enfile = c(Enfile,paste(pname,i,".csv",sep = ""))
# }
# allnum = list.files()
# snpdensity = Enfile[which(!Enfile%in%allnum)] %>% str_extract_all( "\\d+")
# snpdensity = lapply(snpdensity,function(x){
#   x = as.numeric(x[[1]])/10
#   return(x)
# })
# snpdensity = unique(unlist(snpdensity))
#===========================================================================================
####
rm(list=setdiff(ls(),"pname"))
gc()
rep = 10
for(i in 1:length(pname)){
  for(r in 1:rep){
    out= fread(paste(pname[i],r,".csv",sep = ""),sep = ",")
    proname<- pname[i]
    gg0<- out$mean_gv-out$mean_gv[1]
    #水平值（绝对值）的变化可以这样来分段，
    #如果看不同世代区间的增长率（相对值），需要确保基础一致，
    #要用与第0世代的差值作为进展，再用方程同时拟合所有数据再去分段
    #线性模型的系数需要相对一个基准来解释，所以分段比较时需要在一个统一基准下
    genicVariance = out$genicVa
    
    genicall <- sqrt(out$genicVa[1])
    
    genicper <- 1 - sqrt(out$genicVa)/sqrt(genicall)
    
    ggst <- gg0/genicall
    
    dt_cv = data.table(meanGZ = ggst , sdGenicZ = genicper )
    
    fitmodel = lm(meanGZ ~ sdGenicZ, data=dt_cv )
    
    genicEff = coef(fitmodel)[2]
    
    gg5<- out$mean_gv-out$mean_gv[6]
    gg10<- out$mean_gv-out$mean_gv[11]
    gg15<- out$mean_gv-out$mean_gv[16]
    denSNP = (unlist(str_extract_all(pname[i],"\\d+"))[1] %>% as.numeric()) / 10
    nRef = (unlist(str_extract_all(pname[i],"\\d+"))[2] %>% as.numeric())
    if(proname%flike%"IF"){
      imp = "F"
    }else{
      imp = "T"
    }
    inb <- out$inbreeding_plink
    inbp <-out$inbreeding
    acccand <- out$accuracy
    accallcand <- out$acall
    accfam <- out$acforeach
    acimp <- out$acimp
    #gxxe <- out$gxe
    variance <- out$Va
    
    geneticall <- sqrt(out$Va[1])
    
    geneticper <- 1 - sqrt(out$Va)/sqrt(geneticall)
    
    ggstg <- gg0/geneticall
    
    dt_cv = data.table(meanGZ = ggstg , sdGenicZ = geneticper )
    
    fitmodel = lm(meanGZ ~ sdGenicZ, data=dt_cv )
    
    geneticEff = coef(fitmodel)[2]
    
    fitmodel = glm(Va ~ generation, family=Gamma(link="log"), data=out)
    dC = 1 - exp(coef(fitmodel)[2])
    geneticNe = 1/(2*dC)
    
    fitmodel = glm(genicVa ~ generation, family=Gamma(link="log"), data=out)
    
    dC = 1 - exp(coef(fitmodel)[2])
    
    genicNe = 1/(2*dC)
    
    #gg0 = gg0/sqrt(variance[1])
    #gg10 = gg10/sqrt(variance[11])
    #gg20 = gg20/sqrt(variance[21])
    effsize <- out$Ne
    #g2 <- out$g2
    #i2 <- out$i2
    nGeneration<- out$generation
    
    if(i==1 & r==1){
      output<- as.data.table(cbind(proname,nGeneration,imp,denSNP,nRef,ggst,ggstg,genicEff,geneticEff,gg0,gg5,gg10,gg15,inb,acccand,accallcand,accfam,acimp,variance,genicVariance,effsize,genicNe,geneticNe))
    }else{
      output = rbind(output,as.data.table(cbind(proname,nGeneration,imp,denSNP,nRef,ggst,ggstg,genicEff,geneticEff,gg0,gg5,gg10,gg15,inb,acccand,accallcand,accfam,acimp,variance,genicVariance,effsize,genicNe,geneticNe)))
    }
  }
}
colnames(output) = c("pname","nGeneration","imp","denSNP","nRef","ggst","ggstg","genicEff","geneticEff","gg0","gg5","gg10","gg15","inb","acforcand","acforallcand","acforfam","acimp","Va","genicVa","Ne","genicNe","geneticNe")


# temp_dt = output[pname%flike%"Pop12500IF",]
# 
# temp_dt[,imp:="T"]
# 
# substr(temp_dt$pname,1,nchar("Pop12500IF")) <- "Pop12500IT"
# 
# output = rbind(output,temp_dt)

fwrite(output,file ="alldt.csv",sep = ",")


#==================================================================================================================================================================

rm(list = ls())
gc()
output = fread("alldt.csv",sep = ",")

getzt = function(x){
  proname = unique(x$pname)
  for(i in 1:length(proname)){
    for(g in 1:21){
      x[pname == proname[i] & nGeneration == (g-1),gg0se:=sd(x[pname == proname[i] & nGeneration == (g-1),gg0])]
      x[pname == proname[i] & nGeneration == (g-1),gg0:=mean(x[pname == proname[i] & nGeneration == (g-1),gg0])]
      
      x[pname == proname[i] & nGeneration == (g-1),gg5se:=sd(x[pname == proname[i] & nGeneration == (g-1),gg5])]
      x[pname == proname[i] & nGeneration == (g-1),gg5:=mean(x[pname == proname[i] & nGeneration == (g-1),gg5])]
      
      x[pname == proname[i] & nGeneration == (g-1),gg10se:=sd(x[pname == proname[i] & nGeneration == (g-1),gg10])]
      x[pname == proname[i] & nGeneration == (g-1),gg10:=mean(x[pname == proname[i] & nGeneration == (g-1),gg10])]
      
      x[pname == proname[i] & nGeneration == (g-1),gg15se:=sd(x[pname == proname[i] & nGeneration == (g-1),gg15])]
      x[pname == proname[i] & nGeneration == (g-1),gg15:=mean(x[pname == proname[i] & nGeneration == (g-1),gg15])]
      
      x[pname == proname[i] & nGeneration == (g-1),inbse:=sd(x[pname == proname[i] & nGeneration == (g-1),inb])]
      x[pname == proname[i] & nGeneration == (g-1),inb:=mean(x[pname == proname[i] & nGeneration == (g-1),inb])]
      
      x[pname == proname[i] & nGeneration == (g-1),acforcandse:=sd(x[pname == proname[i] & nGeneration == (g-1),acforcand])]
      x[pname == proname[i] & nGeneration == (g-1),acforcand:=mean(x[pname == proname[i] & nGeneration == (g-1),acforcand])]
      
      x[pname == proname[i] & nGeneration == (g-1),acforallcandse:=sd(x[pname == proname[i] & nGeneration == (g-1),acforallcand])]
      x[pname == proname[i] & nGeneration == (g-1),acforallcand:=mean(x[pname == proname[i] & nGeneration == (g-1),acforallcand])]
      
      x[pname == proname[i] & nGeneration == (g-1),acforfamse:=sd(x[pname == proname[i] & nGeneration == (g-1),acforfam])]
      x[pname == proname[i] & nGeneration == (g-1),acforfam:=mean(x[pname == proname[i] & nGeneration == (g-1),acforfam])]
      
      if(any(colnames(x)=="acimp")){
        x[pname == proname[i] & nGeneration == (g-1),acimpse:=sd(x[pname == proname[i] & nGeneration == (g-1),acimp])]
        x[pname == proname[i] & nGeneration == (g-1),acimp:=mean(x[pname == proname[i] & nGeneration == (g-1),acimp])]
      }
      
      x[pname == proname[i] & nGeneration == (g-1),Vase:=sd(x[pname == proname[i] & nGeneration == (g-1),Va])]
      x[pname == proname[i] & nGeneration == (g-1),Va:=mean(x[pname == proname[i] & nGeneration == (g-1),Va])]
      
      x[pname == proname[i] & nGeneration == (g-1),genicVase:=sd(x[pname == proname[i] & nGeneration == (g-1),genicVa])]
      x[pname == proname[i] & nGeneration == (g-1),genicVa:=mean(x[pname == proname[i] & nGeneration == (g-1),genicVa])]
      
      x[pname == proname[i] & nGeneration == (g-1),Nese:=sd(x[pname == proname[i] & nGeneration == (g-1),Ne])]
      x[pname == proname[i] & nGeneration == (g-1),Ne:=mean(x[pname == proname[i] & nGeneration == (g-1),Ne])]
      
      x[pname == proname[i] & nGeneration == (g-1),genicNese:=sd(x[pname == proname[i] & nGeneration == (g-1),genicNe])]
      x[pname == proname[i] & nGeneration == (g-1),genicNe:=mean(x[pname == proname[i] & nGeneration == (g-1),genicNe])]
      
      x[pname == proname[i] & nGeneration == (g-1),geneticNese:=sd(x[pname == proname[i] & nGeneration == (g-1),geneticNe])]
      x[pname == proname[i] & nGeneration == (g-1),geneticNe:=mean(x[pname == proname[i] & nGeneration == (g-1),geneticNe])]
      
      x[pname == proname[i] & nGeneration == (g-1),ggstse:=sd(x[pname == proname[i] & nGeneration == (g-1),ggst])]
      x[pname == proname[i] & nGeneration == (g-1),ggst:=mean(x[pname == proname[i] & nGeneration == (g-1),ggst])]
      
      x[pname == proname[i] & nGeneration == (g-1),ggstgse:=sd(x[pname == proname[i] & nGeneration == (g-1),ggstg])]
      x[pname == proname[i] & nGeneration == (g-1),ggstg:=mean(x[pname == proname[i] & nGeneration == (g-1),ggstg])]
      
      x[pname == proname[i] & nGeneration == (g-1),genicEffse:=sd(x[pname == proname[i] & nGeneration == (g-1),genicEff])]
      x[pname == proname[i] & nGeneration == (g-1),genicEff:=mean(x[pname == proname[i] & nGeneration == (g-1),genicEff])]
      
      x[pname == proname[i] & nGeneration == (g-1),geneticEffse:=sd(x[pname == proname[i] & nGeneration == (g-1),geneticEff])]
      x[pname == proname[i] & nGeneration == (g-1),geneticEff:=mean(x[pname == proname[i] & nGeneration == (g-1),geneticEff])]
    }
  }
  
  return(unique(x))
}

zt1 = getzt(x = output)

zt1$nGeneration = as.numeric(zt1$nGeneration)

fwrite(zt1,file ="allzt.csv",sep = ",")

#==================================================================================================================================================================
rm(list=ls())
gc()

zt1 = fread("allzt.csv",sep = ",")

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

library(showtext)
library(ggplot2)
# showtext_auto(enable = TRUE)


#======================================================
zt = zt1[nRef == 100&imp=="T",]
zt$denSNP = as.character(zt$denSNP)
# zt$denSNP = paste(zt$denSNP," ","snp/chr",sep = "")
levls = c("0","3","12","23","46","114","228","455","682","1250")

ggplot(zt,aes(x=sqrt(genicVa),y=gg0,colour=denSNP))+
  geom_path()+
  geom_point()+
  theme_bw()+
  ylab("Genetic mean")+
  xlab("Genic standard deviation")+
  scale_x_reverse(sec.axis=sec_axis(trans=~1-.,
                                    name="Converted/Lost genic standard deviation"))



#======================================================

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

P<- ggplot(data = zt,aes(x=nRef,y=gg0,group=imp,color = imp))+
  geom_point()+
  geom_line()+
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  xlab("The number of Reference Individuals within each family")+
  ylab("genetic gain")+
  theme_zg()+theme(legend.position=c(0.85, 0.12),legend.title=element_blank())+
  #geom_ribbon(aes(ymin = pgvr -pgvrse, ymax = pgvr+pgvrse), alpha = 0.4)+
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
  geom_errorbar(aes(ymin=gg0-gg0se,
                    ymax=gg0+gg0se),
                width=0.05,alpha = 0.5)+
  scale_x_continuous(limits = c(10,200),breaks = c(10,30,50,70,100,150,200))+labs(color="imputation")+scale_color_aaas()+facet_wrap(~ denSNP)
# P+geom_hline(yintercept = zt[denSNP=="1250"&nRef==200,gg0])
P
ggsave("Figure1.pdf", P, width = 20, height = 10, dpi = 300)


zt$refinterval <- ifelse(zt$nRef <= 50, "10-70", "70-200")
pv = "gg0" #ppt指的硕士答辩的ppt，用的gg0
pvse = "gg0se"
zt$value = zt[,..pv]
zt$valuese = zt[,..pvse]
P<- ggplot(data = zt[imp%in%c("Ped","HD")],aes(x=nRef,y=value,group=1,color = refinterval))+
  geom_point()+
  geom_line()+
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  xlab(NULL)+
  ylab("Genetic gain")+
  theme_zg()+theme(legend.title=element_blank())+
  #geom_ribbon(aes(ymin = pgvr -pgvrse, ymax = pgvr+pgvrse), alpha = 0.4)+
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
  geom_errorbar(aes(ymin=value-valuese,
                    ymax=value+valuese),
                width=0.05,alpha = 0.5)+
  scale_x_continuous(limits = c(10,200),breaks = c(10,30,50,70,100,150,200))+labs(color="imputation")+scale_color_aaas()+facet_wrap(~ denSNP)
# P+geom_hline(yintercept = zt[denSNP=="1250"&nRef==200,gg0])
P
ggsave("Figure1_ppt.pdf", P, width = 10, height = 5, dpi = 300)


pv = "gg0" #ppt指的是给栾老师新作的ppt中的图片
pvse = "gg0se"
zt$value = zt[,..pv]
zt$valuese = zt[,..pvse]
P<- ggplot(data = zt[imp%in%c("Ped","HD")],aes(x=nRef,y=value,group=imp,color = imp))+
  geom_point()+
  geom_line()+
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  xlab(NULL)+
  ylab("Genetic gain")+
  theme_zg()+theme(legend.title=element_blank())+
  #geom_ribbon(aes(ymin = pgvr -pgvrse, ymax = pgvr+pgvrse), alpha = 0.4)+
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
  geom_errorbar(aes(ymin=value-valuese,
                    ymax=value+valuese),
                width=0.05,alpha = 0.5)+
  scale_x_continuous(limits = c(10,200),breaks = c(10,30,50,70,100,150,200))+labs(color="imputation")+scale_color_aaas()#+facet_wrap(~ denSNP)
# P+geom_hline(yintercept = zt[denSNP=="1250"&nRef==200,gg0])
P
ggsave("FigurenewCE_ppt.pdf", P, width = 10, height = 5, dpi = 300)

ggsave("FigurenewGG_ppt.pdf", P, width = 10, height = 5, dpi = 300)



zt = zt1[denSNP==46,]
zt$denSNP = as.character(zt$denSNP)
zt[imp == "F",imp:="LD"]
zt[imp == "T",imp:="Imputation"]
zt$imp = factor(zt$imp,levels = c("LD","Imputation"))
P<- ggplot(data = zt,aes(x=nGeneration,y=acforfam,group=imp,color = imp))+
  geom_point()+
  geom_line()+
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  xlab("nGeneration")+
  ylab("Accuracy within family")+
  theme_zg()+theme(legend.position=c(0.85, 0.12),legend.title=element_blank())+
  #geom_ribbon(aes(ymin = pgvr -pgvrse, ymax = pgvr+pgvrse), alpha = 0.4)+
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
  geom_errorbar(aes(ymin=acforfam-acforfamse,
                    ymax=acforfam+acforfamse),
                width=0.05,alpha = 0.5)+
  scale_x_continuous(limits = c(0,20),breaks = seq(0,20))+labs(color="imputation")+scale_color_aaas()+facet_wrap(~ nRef)
# P+geom_hline(yintercept = zt[denSNP=="1250"&nRef==200,gg0])
P

zt = zt1[nGeneration%in%c(0:20),]
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
zt$GenerationInterval <- cut(zt$nGeneration,
                             breaks = seq(0, 20, by = 5),
                             include.lowest = TRUE,
                             labels = c("0-5", "6-10", "11-15", "16-20"))
P = ggplot(data = zt[imp%in%c("LD"),],aes(x=factor(denSNP),y=gg0,fill = GenerationInterval))+
  # geom_violin() +
  # geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  geom_boxplot() +
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  xlab(" panel size")+
  ylab("Genetic gain")+
  theme_bw()+
  theme_zg()+theme(legend.title=element_blank())+
  # geom_ribbon(aes(ymin = pgvr -pgvrse, ymax = pgvr+pgvrse), alpha = 0.4)+
  # scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
  # geom_errorbar(aes(ymin=gg0-gg0se,
  #                   ymax=gg0+gg0se),
  #               width=0.05,alpha = 0.5)+
  # scale_x_discrete()+
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,5))+
  labs(fill="Generation Interval")+scale_fill_aaas()+facet_wrap(~ GenerationInterval,scales = "free_y")
# dev.off()
P

#====================================================================================================
library(showtext)
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
P = ggplot(data = zt[imp%in%c("LD","HD","Ped","Imputation"),],aes(x=denSNP,y=Ne,group=imp,color = imp))+
  geom_point()+
  geom_line()+
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  xlab(NULL)+
  ylab("Genetic gain")+
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
  geom_errorbar(aes(ymin=Ne-Nese,
                    ymax=Ne+Nese),
                width=0.05,alpha = 0.5)+
  # scale_x_discrete()+
  scale_x_continuous(limits = c(0,1250),breaks = c(0,46,114,228,455,682,1250))+
  labs(color="imputation")+scale_color_aaas()+facet_wrap(~ nRef)
# dev.off()
P
P+geom_hline(yintercept = zt[denSNP=="0",gg0])
ggsave("Figure2.pdf", P, width = 20, height = 10, dpi = 300)

##
# zt = zt1[nGeneration==20&imp=="F",]
# panel_density = unique(zt$denSNP)
# reference_group_size = unique(zt$nRef)
# geneticgain = matrix(nrow = length(panel_density),ncol = length(reference_group_size))
# for (r in 1:length(panel_density)) {
#   for(c in 1:length(reference_group_size)){
#     geneticgain[r,c] = zt[denSNP==panel_density[r]&nRef == reference_group_size[c],gg0]
#   }
# }
# Chkplot = persp3d(panel_density, reference_group_size, geneticgain, col = "lightblue",
#         xlab = "Panel Density", ylab = "Reference Group Size", zlab = "Genetic Gain")
# rglwidget()
# rgl.postscript("3D_plot.pdf", fmt = "pdf")
# rgl.snapshot("3D_plot.png")
# writeWebGL("3D_plot.html")
##
P = ggplot(data = zt[nRef == "70" &imp%in%c("LD","HD"),],aes(x=denSNP,y=gg0,group=1,color = snpinterval))+
  geom_point()+
  geom_line()+
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  xlab(NULL)+
  ylab("Genetic gain")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text = element_text(family = "STXihei"),
        legend.position = "right",
        legend.background = element_rect(colour = "black"))+
  theme_zg()+theme(legend.title=element_blank())+
  #geom_ribbon(aes(ymin = pgvr -pgvrse, ymax = pgvr+pgvrse), alpha = 0.4)+
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
  geom_errorbar(aes(ymin=gg0-gg0se,
                    ymax=gg0+gg0se),
                width=0.05,alpha = 0.5)+
  # scale_x_discrete()+
  scale_x_continuous(limits = c(3,1250),breaks = c(3,46,1250))+
  labs(color="SNPinterval")+scale_color_aaas()#+facet_wrap(~ nRef)
# dev.off()
P
ggsave("Figure2_ppt.pdf", P, width = 10, height = 5, dpi = 300)

# zt[denSNP>0&denSNP<=46,snpinterval:="3-46"]
# zt[denSNP>46&denSNP<=1250,snpinterval:="114-682"]
levls = c("0","3","12","23","46","114","228","455","682","1250")
zt$denSNP = factor(zt$denSNP,levels = levls)
# zt$snpinterval = factor(zt$snpinterval,levels = c("3-46","116-482"))
P = ggplot(data = zt[nRef == "70" &imp%in%c("LD","Imputation"),],aes(x=denSNP,y=gg0,group=imp,color = snpinterval,shape = imp))+
  geom_point()+
  geom_line()+
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  xlab(NULL)+
  ylab("Genetic gain")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text = element_text(family = "STXihei"),
        legend.position = "right",
        legend.background = element_rect(colour = "black"))+
  theme_zg()+theme(legend.title=element_blank())+
  #geom_ribbon(aes(ymin = pgvr -pgvrse, ymax = pgvr+pgvrse), alpha = 0.4)+
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
  geom_errorbar(aes(ymin=gg0-gg0se,
                    ymax=gg0+gg0se),
                width=0.05,alpha = 0.5)+
  scale_x_discrete()+
  # scale_x_continuous(limits = c(3,682))+
  labs(color="SNPinterval")+scale_color_aaas()#+facet_wrap(~ snpinterval)
# dev.off()
P
ggsave("Figure2_imputation_ppt.pdf", P, width = 10, height = 5, dpi = 300)

zt = zt1[nGeneration%in%c(0:20),]
#zt$nRef = as.character(zt$nRef)
#zt$nRef = factor(zt$nRef,levels = c("10","30","50","70","100","150","200"))
zt[denSNP=="0",imp:="Ped"]
zt[denSNP=="1250",imp:="HD"]
zt[imp == "F",imp:="LD"]
zt[imp == "T",imp:="Imputation"]
zt$imp = factor(zt$imp,levels = c("Ped","LD","Imputation","HD"))
zt$GenerationInterval <- cut(zt$nGeneration,
                             breaks = seq(0, 20, by = 5),
                             include.lowest = TRUE,
                             labels = c("0-5", "6-10", "11-15", "16-20"))
# zt$denSNP = as.character(zt$denSNP)
# zt$denSNP = factor(zt$denSNP,levels = c("3","12","23","46","114","228","455","682","1250"))
# pdf("Figure2.pdf", width = 20, height = 10)
# showtext_auto(enable = TRUE)
P = ggplot(data = zt[imp%in%c("HD"),],aes(x=nRef,y=gg0,group=nGeneration,color = GenerationInterval))+
  geom_point()+
  geom_line()+
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  xlab("The number of SNP per chromosome")+
  ylab("The genetic gain")+
  theme_bw()+
  theme_zg()+theme(legend.position=c(0.85, 0.12),legend.title=element_blank())+
  #geom_ribbon(aes(ymin = pgvr -pgvrse, ymax = pgvr+pgvrse), alpha = 0.4)+
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
  geom_errorbar(aes(ymin=gg0-gg0se,
                    ymax=gg0+gg0se),
                width=0.05,alpha = 0.5)+
  # scale_x_discrete()+
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,5))+
  labs(color="imputation")+scale_color_aaas()+facet_wrap(~ GenerationInterval)
# dev.off()
P

P = ggplot(data = zt[imp%in%c("HD"),],aes(x=factor(nRef),y=gg0,fill = GenerationInterval))+
  geom_violin() +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  # geom_boxplot() +
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  xlab(" reference group size")+
  ylab("The genetic gain")+
  theme_bw()+
  theme_zg()+theme(legend.title=element_blank())+
  # geom_ribbon(aes(ymin = pgvr -pgvrse, ymax = pgvr+pgvrse), alpha = 0.4)+
  # scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
  # geom_errorbar(aes(ymin=gg0-gg0se,
  #                   ymax=gg0+gg0se),
  #               width=0.05,alpha = 0.5)+
  # scale_x_discrete()+
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,5))+
  labs(fill="Generation Interval")+scale_fill_aaas()+facet_wrap(~ GenerationInterval,scales = "free_y")
# dev.off()
P

zt %>%
  group_by(GenerationInterval, nRef) %>%
  summarize(MeanGeneticProgress = mean(gg0))

#====================================================================================
zt = zt1[nGeneration==20,]
zt$nRef = as.character(zt$nRef)
zt$denSNP = as.character(zt$denSNP)

zt[denSNP=="0",imp:="Ped"]

zt[denSNP=="1250",imp:="HD"]

zt[imp == "F",imp:="LD"]

zt[imp == "T",imp:="Imputation"]

zt$imp = factor(zt$imp,levels = c("Ped","LD","Imputation","HD"))

zt$denSNP = factor(zt$denSNP,levels = c("0","3","12","23","46","114","228","455","682","1250"))

zt$nRef = factor(zt$nRef,levels = c("10","30","50","70","100","150","200"))

P<- ggplot(data = zt,aes(x = denSNP, y = gg0, fill = imp))+
  geom_bar(stat = "identity",position = "dodge",)+
  ylab("Genetic gain")+xlab("The number of SNP per chromosome")+
  geom_errorbar(aes(ymax = gg0+gg0se, ymin = gg0-gg0se),
                position = position_dodge(0.9), width = 0.15)+
  theme_zg()+facet_wrap(~nRef)+ggsci::scale_fill_aaas()+theme(legend.position=c(0.85, 0.12),legend.title=element_blank())
P
ggsave("Figure3.pdf", P, width = 20, height = 10, dpi = 300)

P<- ggplot(data = zt[nRef=="50",],aes(x = denSNP, y = gg0, fill = imp))+
  geom_bar(stat = "identity",position = "dodge",)+
  ylab("Genetic gain")+xlab("The number of SNP per chromosome")+
  geom_errorbar(aes(ymax = gg0+gg0se, ymin = gg0-gg0se),
                position = position_dodge(0.9), width = 0.15)+
  theme_zg()+facet_wrap(~nRef)+ggsci::scale_fill_aaas()+theme(legend.title=element_blank())
P
# P<- ggplot(data = zt,aes(x = nRef, y = gg0, fill = imp))+
#   geom_bar(stat = "identity",position = "dodge",)+
#   ylab("Genetic gain")+
#   geom_errorbar(aes(ymax = gg0+gg0se, ymin = gg0-gg0se),
#                 position = position_dodge(0.9), width = 0.15)+
#   theme_zg()+facet_wrap(~denSNP)
# P
#====================================================================================
zt = zt1[nRef == 100&imp=="F",]
P1<- ggplot(data = zt,aes(x=nGeneration,y=Ne,group=pname,color = pname))+
  geom_point()+
  geom_line()+
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  xlab("Generation")+
  ylab("Genomic inbreedingr")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text = element_text(),
        legend.position = "right",
        legend.background = element_rect(colour = "black"))+
  #geom_ribbon(aes(ymin = pgvr -pgvrse, ymax = pgvr+pgvrse), alpha = 0.4)+
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
  geom_errorbar(aes(ymin=Ne-Nese,
                    ymax=Ne+Nese),
                width=0.05,alpha = 0.5)+
  scale_x_continuous(limits = c(1,20),breaks = seq(1,20,1))+labs(color="Breeding scheme")+scale_color_aaas()#+facet_wrap(~ ebv_calcualtion)
#scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(linetype="Breeding system",shape="Proportion of\nthe introduced\nMP individuals (%)")
P1

#=========================================================================================

#===========================create the table===========================
library(tidyverse)
library(knitr)
library(kableExtra)

# Prepare your data

createDT = function(zt1,Value,Valuese){
  SnpDensity = unique(zt1[nGeneration==20&imp=="F",denSNP])
  RefSize = unique(zt1[nGeneration==20&imp=="F",nRef])
  V1 = c()
  V1SE = c()
  V2 = c()
  V2SE = c()
  f1 = c()
  f2 = c()
  for (s in SnpDensity) {
    for (r in RefSize) {
      value1 = zt1[nGeneration==20&imp=="F"&denSNP==s&nRef==r,..Value]
      valueSE1 = zt1[nGeneration==20&imp=="F"&denSNP==s&nRef==r,..Valuese]
      value1 = value1 %>% as.numeric() %>% round(2)
      valueSE1 = valueSE1 %>% as.numeric() %>% round(2)
      
      V1 = c(V1,value1)
      V1SE = c(V1SE,valueSE1)
      
      value2 = zt1[nGeneration==20&imp=="T"&denSNP==s&nRef==r,..Value]
      valueSE2 = zt1[nGeneration==20&imp=="T"&denSNP==s&nRef==r,..Valuese]
      value2 = value2 %>% as.numeric() %>% round(2)
      valueSE2 = valueSE2 %>% as.numeric() %>% round(2)
      
      V2 = c(V2,value2)
      V2SE = c(V2SE,valueSE2)
      
      f1 = c(f1,s)
      f2 = c(f2,r)
    }
  }
  dt_table <- data.frame(
    Factor1 = f1,
    Factor2 = f2,
    Value1 = V1,
    SE1 = V1SE,
    Value2 = V2,
    SE2 = V2SE
  )
  return(dt_table)
}
Value = "gg0"
Valuese = "gg0se"

dt_table = createDT(zt1 = zt1, Value = Value,Valuese = Valuese)


# Combine Value1, SE1, Value2, and SE2 into a single column
dt_table <- dt_table %>%
  mutate(Values = paste(Value1, "±", SE1, "/", Value2, "±", SE2))

# Reshape the data into a matrix format
data_matrix <- dt_table %>%
  dplyr::select(Factor1, Factor2, Values) %>%
  tidyr::spread(key = Factor2, value = Values)

# Rename the rows to match Factor1
rownames(data_matrix) <- data_matrix$Factor1
data_matrix <- data_matrix[, -1]

table <- data_matrix %>%
  knitr::kable(
    caption = "Two Values and Standard Errors for Combinations of Factor1 and Factor2 in Matrix Format",
    format = "html"  # Change to "latex" if needed
  ) %>%
  kableExtra::kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    position = "center"
  )
table
# Save the formatted HTML table
save_kable(table, file = "table.html")

# Additionally, save as CSV and Excel
write.csv(data_matrix, "table.csv", row.names = TRUE)

install.packages("remarkdown")
install.packages("flextable")
install.packages("officer")
remotes::install_github("davidgohel/flextable")
remotes::install_github("davidgohel/officer")
library()
library(officer)
ft <- flextable(data_matrix)

# 创建一个新的 Word 文档
doc <- read_docx()

# 将 flextable 插入到 Word 文档
doc <- body_add_flextable(doc, value = ft)

# 保存 Word 文档
print(doc, target = "table.docx")

#
dt_table = as.data.table(dt_table)
dt_table[Factor1==0,Value1]
#
system("python3 3dplot.py --ITYPE F --VTYPE genicEff --VSE genicEffse --plotname 3dplotFgenicEff.pdf --ZTYPE 'Conversion Efficiency' ")

system("python3 3dplot.py --ITYPE T --VTYPE genicEff --VSE genicEffse --plotname 3dplotTgenicEff.pdf --ZTYPE 'Conversion Efficiency' ")

system("python3 3dplot.py --ITYPE F --VTYPE gg0 --VSE gg0se --plotname 3dplotFgeneticgain.pdf --ZTYPE 'Genetic gain' ")

system("python3 3dplot.py --ITYPE T --VTYPE genicVa --VSE genicVase --plotname 3dplotTgenicVa.pdf --ZTYPE 'Genetic diversity' ")
system("python3 3dplot.py --ITYPE F --VTYPE genicVa --VSE genicVase --plotname 3dplotFgenicVa.pdf --ZTYPE 'Genetic diversity' ")

system("python3 maketable.py --VTYPE gg0 --VSE gg0se --TNAME geneticgaintable.docx")

#=================================================================================

#====================model and contrast=====================================================================
g = 20
im = "F"
denS = c(3,12,23,46,114,228,455,682,1250)
Ref = c(10,30,50,70,100,150,200)
zt = zt1[nGeneration==g,]
calcontra = function(x,y){
  return((((x - y) / y)*100) %>% mean() %>% round(digits = 2))
}
library(magrittr)

#cal for all breeding schemes==
all10 = zt[nRef == 10,.(pname,nRef,gg0)]$gg0 %>% round(2)
all50 = zt[nRef == 50,.(pname,nRef,gg0)]$gg0 %>% round(2)
all100 = zt[nRef == 100,.(pname,nRef,gg0)]$gg0 %>% round(2)
calcontra(x = all50 , y = all10)
calcontra(x = all100 , y = all50)
#=============================
###=======calpara=================calpara==========
gfor10 = zt[nRef == 10 & imp == im & denSNP %in% denS,.(pname,nRef,gg0)]$gg0
gfor30 = zt[nRef == 30 & imp == im & denSNP %in% denS,.(pname,nRef,gg0)]$gg0
gfor50 = zt[nRef == 50 & imp == im & denSNP %in% denS,.(pname,nRef,gg0)]$gg0
gfor70 = zt[nRef == 70 & imp == im & denSNP %in% denS,.(pname,nRef,gg0)]$gg0
gfor100 = zt[nRef == 100 & imp == im & denSNP %in% denS,.(pname,nRef,gg0)]$gg0
gfor150 = zt[nRef == 150 & imp == im & denSNP %in% denS,.(pname,nRef,gg0)]$gg0
gfor200 = zt[nRef == 200 & imp == im & denSNP %in% denS,.(pname,nRef,gg0)]$gg0
calcontra(x = gfor30 , y = gfor10)
calcontra(x = gfor50 , y = gfor30)
calcontra(x = gfor70 , y = gfor50)
calcontra(x = gfor100 , y = gfor70)
calcontra(x = gfor100 , y = gfor70)
calcontra(x = gfor150 , y = gfor100)
calcontra(x = gfor200 , y = gfor100)

gfor10lf = zt[nRef == 10 & imp == "F" & denSNP %in% c(3,12,23,46),.(pname,nRef,gg0)]$gg0
gfor10lt = zt[nRef == 10 & imp == "T" & denSNP %in% c(3,12,23,46),.(pname,nRef,gg0)]$gg0
calcontra(x = gfor10lt , y = gfor10lf)
gfor200lf = zt[nRef == 200 & imp == "F" & denSNP %in% c(3,12,23,46),.(pname,nRef,gg0)]$gg0
gfor200lt = zt[nRef == 200 & imp == "T" & denSNP %in% c(3,12,23,46),.(pname,nRef,gg0)]$gg0
calcontra(x = gfor200lt , y = gfor200lf)

gfor70l = zt[nRef == 70 & imp == im & denSNP %in% c(3,12,23,46),.(pname,nRef,gg0)]$gg0
gfor150l = zt[nRef == 150 & imp == im & denSNP %in% c(3,12,23,46),.(pname,nRef,gg0)]$gg0
calcontra(x = gfor150l , y = gfor70l)

gfor70h = zt[nRef == 70 & imp == im & denSNP %in% c(114,228,455,682),.(pname,nRef,gg0)]$gg0
gfor150h = zt[nRef == 150 & imp == im & denSNP %in% c(114,228,455,682),.(pname,nRef,gg0)]$gg0
calcontra(x = gfor150h , y = gfor70h)
#===========================================================================================
library(interactions)
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

####
anova_model<- aov(ld~pname+genoratio+gev,data =dt_test )
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"pname",conf.level = 0.95)

##


#==========model of all factor for genomic selection=============================================================================================
outGS = output[denSNP!=0,]
outGS$SnpInterval = factor(outGS$SnpInterval,levels = c("3-46","114-1250"))
outGS$RefInterval = factor(outGS$RefInterval,levels = c("10-70","100-200"))
outGS$SnpInterval %>% unique()
outGS$RefInterval %>% unique()
outGS[nGeneration<=10,GenerationInterval:="1-10"]
outGS[nGeneration>=11,GenerationInterval:="11-20"]
# outGS = outGS[-which(outGS$gg0<0)]
allmodel = lm(log(gg0)~denSNP+nRef+nGeneration*GenerationInterval+denSNP*GenerationInterval+nRef*GenerationInterval+nGeneration*SnpInterval+nRef*SnpInterval+denSNP*SnpInterval+nGeneration*RefInterval+nRef*RefInterval+denSNP*RefInterval+nGeneration,data = outGS[nGeneration %in%c(1:20) & imp == "F",])
summary(allmodel)

anova_model<- aov(gg0~factor(nRef)+factor(denSNP),data =outGS[nGeneration==20&imp=="F",] )
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"factor(denSNP)",conf.level = 0.95)
tukey_results

backward_model <- step(allmodel, direction = "backward")
summary(backward_model)

library(car)
vif(allmodel)
coeValue = allmodel$coefficients["nGeneration"] #+ allmodel$coefficients["denSNP:SnpInterval114-1250"]
coeValue*100
(exp(coeValue)-1)*100

interact_plot(allmodel, pred = denSNP, modx = GenerationInterval, plot.points = FALSE) +
  labs(title = "交互效应图", x = "densnp", y = "预期 y 值")+
  theme_zg()+scale_x_continuous(limits = c(3,1250),breaks = c(3,46,1250))

chk = summary(allmodel) 
chk = as.data.frame(chk$coefficients)
coefficients <- data.frame(
  term = rownames(chk),
  estimate = chk$Estimate,
  std_error = chk$`Std. Error`
)

# 添加置信区间
coefficients <- coefficients %>%
  mutate(
    lower_ci = estimate - 1.96 * std_error,
    upper_ci = estimate + 1.96 * std_error
  )
coefficients = as.data.table(coefficients)
# 绘制棘状图
ggplot(coefficients[!term%in%c("SnpInterval114-1250","RefInterval100-200","nGeneration","GenerationInterval11-20","(Intercept)"),], aes(x = term, y = estimate)) +
  geom_point(size = 4) +
  geom_segment(aes(x = term, xend = term, y = lower_ci, yend = upper_ci), color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "回归系数及其95%置信区间", x = "变量", y = "估计值") +
  coord_flip() +
  theme_minimal()+theme_zg()

ggplot(coefficients[!term%in%c("SnpInterval114-1250","RefInterval100-200","nGeneration","GenerationInterval11-20","(Intercept)"),], aes(x = term, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "回归系数及其95%置信区间", x = "变量", y = "估计值") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


chk = summary(allmodel) 

coefficients <- as.data.frame(chk$coefficients)
coefficients <- coefficients %>%
  mutate(
    term = rownames(coefficients),
    lower_ci = Estimate - 1.96 * `Std. Error`,
    upper_ci = Estimate + 1.96 * `Std. Error`
  ) %>%
  dplyr::select(term, Estimate, `Std. Error`, lower_ci, upper_ci)

library(kableExtra)
# 制作表格
kable(coefficients, format = "html", col.names = c("Term", "Estimate", "Std. Error", "Lower CI", "Upper CI"), digits = 3) %>%
  kable_styling(full_width = F, position = "center") %>%
  column_spec(1, bold = TRUE) %>%
  column_spec(2:5, width = "2cm")
#=========================================================================


allmodel = lm(log(gg0)~denSNP+nRef+denSNP*GenerationInterval+nRef*GenerationInterval+nRef*SnpInterval+denSNP*SnpInterval+nRef*RefInterval+denSNP*RefInterval+nGeneration,data = outGS[nGeneration %in%c(1:20) & imp == "T",])
summary(allmodel)

allmodel = lm(log(Ne)~denSNP+nRef+denSNP*GenerationInterval+nRef*GenerationInterval+nRef*SnpInterval+denSNP*SnpInterval+nRef*RefInterval+denSNP*RefInterval+nGeneration,data = outGS[nGeneration %in%c(1:20) & imp == "F",])
summary(allmodel)

allmodel = lm(log(Ne)~denSNP+nRef+nGeneration+GenerationInterval+SnpInterval+RefInterval+SnpInterval*RefInterval+SnpInterval*GenerationInterval+RefInterval*GenerationInterval,data = outGS[nGeneration %in%c(1:20) & imp == "F",])
summary(allmodel)

allmodel = lm(Ne~imp + imp*SnpInterval+ imp*RefInterval+ denSNP+nRef+nGeneration*GenerationInterval+denSNP*GenerationInterval+nRef*GenerationInterval+nGeneration*SnpInterval+nRef*SnpInterval+denSNP*SnpInterval+nGeneration*RefInterval+nRef*RefInterval+denSNP*RefInterval+nGeneration,data = outGS[nGeneration %in%c(1:20),])
summary(allmodel)

backward_model <- step(allmodel, direction = "backward")
summary(backward_model)


allmodel = lm(log(Ne)~denSNP+nRef+denSNP*GenerationInterval+nRef*GenerationInterval+nRef*SnpInterval+denSNP*SnpInterval+nRef*RefInterval+denSNP*RefInterval+nGeneration,data = outGS[nGeneration %in%c(1:20) & imp == "T",])
summary(allmodel)

backward_model <- step(allmodel, direction = "backward")
summary(backward_model)

emm <- emmeans::emmeans(allmodel, ~ nRef | GenerationInterval)


#======================================================================================================


#==========model of all factor for pedigree selection=============================================================================================
outPS = output[denSNP==0,]
outPS$RefInterval = factor(outPS$RefInterval,levels = c("10-70","100-200"))
outGS$RefInterval %>% unique()
# outGS = outGS[-which(outGS$gg0<0)]
allmodel = lm(log(gg0)~nRef+nRef*GenerationInterval+nRef*RefInterval+nGeneration,data = outPS[nGeneration %in%c(1:20) & imp == "F",])
summary(allmodel)

backward_model <- step(allmodel, direction = "backward")
summary(backward_model)

emm <- emmeans::emmeans(allmodel, ~ nRef | GenerationInterval)


#======================================================================================================

#==========cal gg0 ～ ref=============================================================================================


g = 20
modellow = lm(log(gg0)~nRef,data = output[nGeneration==g & nRef %in% c(70,100,150,200) & imp == "F" & denSNP ==0,])
summary(modellow )
coeValue = modellow$coefficients["nRef"] #+ allmodel$coefficients["denSNP:SnpInterval114-1250"]
coeValue*100
(exp(coeValue)-1)*100

modelhig = lm(log(gg0)~nRef,data = output[nGeneration==g & nRef %in% c(10,30,50,70) & imp == "F" & denSNP ==1250,])
summary(modelhig)
coeValue = modelhig$coefficients["nRef"] #+ allmodel$coefficients["denSNP:SnpInterval114-1250"]
coeValue*100
(exp(coeValue)-1)*100

g = 20
modelhig = lm(log(gg0)~nRef,data = output[nGeneration%in%c((g-4):g) & nRef %in% c(10,30,50,70) & imp == "F" & denSNP ==1250,])
summary(modelhig)

modelallg = lm(log(gg0)~nRef,data = output[nGeneration==g & nRef %in% c(10,30,50,70,100,150,200) & denSNP ==0,])
summary(modelallg)

modelallg = lm(log(gg0)~nRef,data = output[nGeneration==g & nRef %in% c(10,30,50,70,100,150,200) & denSNP !=0 ,])
summary(modelallg)

coeValue = modelallg$coefficients["nRef"] #+ allmodel$coefficients["denSNP:SnpInterval114-1250"]
coeValue*100
(exp(coeValue)-1)*100


im = "F"
im = "T"
modellow = lm(log(gg0)~nRef,data = output[nGeneration==g & nRef %in% c(10,30,50,70,100,150,200) & imp == im & denSNP %in% c(3,12,23,46),])
summary(modellow )

modelhig = lm(log(gg0)~nRef,data = output[nGeneration==g & nRef %in% c(10,30,50,70,100,150,200) & imp == im & denSNP %in% c(46,114,228,455,682),])
summary(modelhig)


#======================================================================================================

#==========cal gg0 ～ dens=============================================================================================
g = 20
im = "F"
denS = c(3,12,23,46,114,228,455,682,1250)
Ref = c(10,30,50,70,100,150,200)
modellow = lm(log(gg0)~denSNP,data = output[nGeneration==g & nRef %in% c(10,30,50,70,100,150,200) & imp == im & denSNP !=0,])
summary(modellow )

modellow = lm(log(gg0)~denSNP,data = output[nGeneration==g & nRef %in% c(10,30,50,70,100,150,200) & imp == im & denSNP %in% c(3,12,23,46),])
summary(modellow )

modelhig = lm(log(gg0)~denSNP,data = output[nGeneration==g & nRef %in% c(10,30,50,70,100,150,200) & imp == im & denSNP %in% c(46,114,228,455,682,1250),])
summary(modelhig)


#======================================================================================================

mymodel<- aov(gg0~nRef,data= output[nGeneration==g & nRef %in% c(50,70,100) & imp == im & denSNP %in% denS,])

summary(mymodel)

t.test(x = output[nGeneration==g & nRef == 70 & imp == im & denSNP %in% denS,gg0],y=output[nGeneration==g & nRef == 50 & imp == im & denSNP %in% denS,gg0])

t.test(x = output[nGeneration==g & nRef == 100 & imp == im & denSNP %in% denS,gg0],y=output[nGeneration==g & nRef == 50 & imp == im & denSNP %in% denS,gg0])

t.test(x = output[nGeneration==g & nRef == 100 & imp == im & denSNP %in% denS,gg0],y=output[nGeneration==g & nRef == 70 & imp == im & denSNP %in% denS,gg0])

t.test(x = output[nGeneration==g & nRef == 100 & imp == im & denSNP %in% denS,gg0],y=output[nGeneration==g & nRef == 150 & imp == im & denSNP %in% denS,gg0])


#=======calpara=================calpara==========
gfor3 = zt[denSNP == 3 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor12 = zt[denSNP == 12 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor23 = zt[denSNP == 23 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor46 = zt[denSNP == 46 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor114 = zt[denSNP == 114 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor228 = zt[denSNP == 228 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor455 = zt[denSNP == 455 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor682 = zt[denSNP == 682 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
if(im == "F"){
  gfor1250 = zt[denSNP == 1250 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0 
}
calcontra(x = gfor12 , y = gfor3)
calcontra(x = gfor23 , y = gfor12)
calcontra(x = gfor46 , y = gfor23)
calcontra(x = gfor114 , y = gfor46)
calcontra(x = gfor114 , y = gfor23)
calcontra(x = gfor228 , y = gfor114)
calcontra(x = gfor455 , y = gfor114)
calcontra(x = gfor682 , y = gfor114)
calcontra(x = gfor1250 , y = gfor114)
calcontra(x = gfor455 , y = gfor228)
calcontra(x = gfor682 , y = gfor228)
calcontra(x = gfor1250, y = gfor228)
calcontra(x = gfor1250, y = gfor455)
calcontra(x = gfor1250, y = gfor682)

output = fread("alldt.csv",sep = ",")
t.test(x = output[nGeneration==g & denSNP == 0 & imp == im & nRef %in% Ref,inb],y=output[nGeneration==g & denSNP == 23 & imp == im & nRef %in% Ref,inb])


t.test(x = output[nGeneration==g & denSNP == 46 & imp == im & nRef %in% Ref,gg0],y=output[nGeneration==g & denSNP == 114 & imp == im & nRef %in% Ref,gg0])
t.test(x = output[nGeneration==g & denSNP == 114 & imp == im & nRef %in% Ref,gg0],y=output[nGeneration==g & denSNP == 228 & imp == im & nRef %in% Ref,gg0])
t.test(x = output[nGeneration==g & denSNP == 114 & imp == im & nRef %in% Ref,gg0],y=output[nGeneration==g & denSNP == 455 & imp == im & nRef %in% Ref,gg0])
t.test(x = output[nGeneration==g & denSNP == 228 & imp == im & nRef %in% Ref,gg0],y=output[nGeneration==g & denSNP == 455 & imp == im & nRef %in% Ref,gg0])
t.test(x = output[nGeneration==g & denSNP == 228 & imp == im & nRef %in% Ref,gg0],y=output[nGeneration==g & denSNP == 682 & imp == im & nRef %in% Ref,gg0])
t.test(x = output[nGeneration==g & denSNP == 228 & imp == im & nRef %in% Ref,gg0],y=output[nGeneration==g & denSNP == 1250 & imp == im & nRef %in% Ref,gg0])
t.test(x = output[nGeneration==g & denSNP == 455 & imp == im & nRef %in% Ref,gg0],y=output[nGeneration==g & denSNP == 682 & imp == im & nRef %in% Ref,gg0])
t.test(x = output[nGeneration==g & denSNP == 455 & imp == im & nRef %in% Ref,gg0],y=output[nGeneration==g & denSNP == 1250 & imp == im & nRef %in% Ref,gg0])
#=======calpara=================calpara==========
gfor3T = zt[denSNP == 3 & imp == "T" & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor3F = zt[denSNP == 3 & imp == "F" & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
calcontra(x = gfor3T, y = gfor3F)
(zt[denSNP == 3 & imp == "T" & nRef %in% Ref,acimp]*100) %>% mean() %>% round(digits = 2)

gfor12T = zt[denSNP == 12 & imp == "T" & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor12F = zt[denSNP == 12 & imp == "F" & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
calcontra(x = gfor12T, y = gfor12F)
(zt[denSNP == 12 & imp == "T" & nRef %in% Ref,acimp]*100) %>% mean() %>% round(digits = 2)

gfor23T = zt[denSNP == 23 & imp == "T" & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor23F = zt[denSNP == 23 & imp == "F" & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
calcontra(x = gfor23T, y = gfor23F)
(zt[denSNP == 23 & imp == "T" & nRef %in% Ref,acimp]*100) %>% mean() %>% round(digits = 2)

gfor46T = zt[denSNP == 46 & imp == "T" & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor46F = zt[denSNP == 46 & imp == "F" & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
calcontra(x = gfor46T, y = gfor46F)
(zt[denSNP == 46 & imp == "T" & nRef %in% Ref,acimp]*100) %>% mean() %>% round(digits = 2)


gfor228T = zt[denSNP == 228 & imp == "T" & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor228F = zt[denSNP == 228 & imp == "F" & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
calcontra(x = gfor46T, y = gfor46F)
(zt[denSNP == 228 & imp == "T" & nRef %in% Ref,acimp]*100) %>% mean() %>% round(digits = 2)

output = fread("alldt.csv",sep = ",")
t.test(x = output[nGeneration==g & denSNP == 46 & imp == "T" & nRef %in% Ref,gg0],y=output[nGeneration==g & denSNP == 46 & imp == "F" & nRef %in% Ref,gg0])
t.test(x = output[nGeneration==g & denSNP == 114 & imp == "T" & nRef %in% Ref,gg0],y=output[nGeneration==g & denSNP == 114 & imp == "F" & nRef %in% Ref,gg0])
t.test(x = output[nGeneration==g & denSNP == 228 & imp == "T" & nRef %in% Ref,gg0],y=output[nGeneration==g & denSNP == 228 & imp == "F" & nRef %in% Ref,gg0])
t.test(x = output[nGeneration==g & denSNP == 455 & imp == "T" & nRef %in% Ref,gg0],y=output[nGeneration==g & denSNP == 455 & imp == "F" & nRef %in% Ref,gg0])
t.test(x = output[nGeneration==g & denSNP == 682 & imp == "T" & nRef %in% Ref,gg0],y=output[nGeneration==g & denSNP == 682 & imp == "F" & nRef %in% Ref,gg0])
#=======calpara=================calpara==========
#compared with/within pblup/ssgblup
g=20
Ref = c(10,30,50,70,100,150,200)
Ref = 10
im = "T"
t.test(x = output[nGeneration==g & nRef == 70  & denSNP == 0,gg0],y=output[nGeneration==g & nRef == 100 & denSNP == 0,gg0])

gforp = zt[nRef %in% Ref  & denSNP == 0,gg0]
gfor3 = zt[denSNP == 3 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor12 = zt[denSNP == 12 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor23 = zt[denSNP == 23 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor46 = zt[denSNP == 46 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor114 = zt[denSNP == 114 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor228 = zt[denSNP == 228 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor455 = zt[denSNP == 455 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor682 = zt[denSNP == 682 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
gfor1250 = zt[denSNP == 1250 & imp == im & nRef %in% Ref,.(pname,denSNP,gg0)]$gg0
calcontra(x = gfor3,y = gforp)
calcontra(x = gfor12,y = gforp)
calcontra(x = gfor1250,y = gforp)
gforp10 = zt[nRef == 10  & denSNP == 0,gg0]
gforp30 = zt[nRef == 30  & denSNP == 0,gg0]
gforp50 = zt[nRef == 50  & denSNP == 0,gg0]
gforp70 = zt[nRef == 70  & denSNP == 0,gg0]
gforp150 = zt[nRef == 150  & denSNP == 0,gg0]
calcontra(x = gforp30,y = gforp10)
calcontra(x = gforp70,y = gforp30)
calcontra(x = gforp150,y = gforp70)


t.test(x = output[nGeneration==g & nRef == 70  & denSNP == 1250,gg0],y=output[nGeneration==g & nRef == 100 & denSNP == 1250,gg0])

gfors10 = zt[nRef == 10  & denSNP == 1250,gg0]
gfors30 = zt[nRef == 30  & denSNP == 1250,gg0]
gfors50 = zt[nRef == 50  & denSNP == 1250,gg0]
gfors70 = zt[nRef == 70  & denSNP == 1250,gg0]
gfors100 = zt[nRef == 100  & denSNP == 1250,gg0]
gfors150 = zt[nRef == 150  & denSNP == 1250,gg0]
calcontra(x = gfors30,y = gfors10)
calcontra(x = gfors70,y = gfors30)
calcontra(x = gfors150,y = gfors70)

gfors150 = zt[nRef == 150  & denSNP == 1250,gg0]
gfors7046F = zt[nRef == 70  & denSNP == 46 &imp == "F",gg0]
gfors7046T = zt[nRef == 70  & denSNP == 46 &imp == "T",gg0]
gfors5046T = zt[nRef == 50  & denSNP == 46 &imp == "T",gg0]
calcontra(x = gfors7046F,y = gfors70)
calcontra(x = gfors7046T,y = gfors70)
calcontra(x = gfors7046T,y = gfors150)
calcontra(x = gfors5046T,y = gfors70)
