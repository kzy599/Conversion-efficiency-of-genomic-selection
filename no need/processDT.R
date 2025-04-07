# install.packages("kableExtra")
# devtools::install_github("tidyverse/ggplot2")
setwd("/home/GTDisk1/kangziyi/ProDTc/")
rm(list = ls())
gc()
# pname = paste("pop",c(0:6),sep = "")
# pname_s = paste("spop",c(1:6),sep = "")
# pname = c(pname,pname_s)
allnum = list.files()
allnum = allnum[allnum%flike%"10.csv"]
pname = substr(allnum,1,nchar(allnum)-nchar("**.csv"))
pname = unique(pname)
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

    accallcand <- out$acall
    accfam <- out$acforeach

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

    nGeneration<- out$generation
    
    ac = unlist(str_extract_all(pname[i],"\\d+")) %>% as.numeric()
    ac <- c(0, 1, 0.1, 0.3, 0.5, 0.7, 0.9)[match(ac, c(0, 1, 2, 3, 4, 5, 6))]
    
    
    if(pname[i]%flike%"spop"){
      if(ac == 0){
        relM = "G"
      }else{
        relM = "A"
      }
    }else{
      if(ac == 0){
        relM = "A"
      }else{
        relM = "G"
      }
    }
    
    tempchar = substr(pname[i],1,1)
    
    if(tempchar=="l"){
      TD = 25
    }else if(tempchar=="m"){
      TD = 65
    }else if(tempchar == "h"){
      TD = 90
    }else if(tempchar == "n"){
      TD = "rand"
      relM = "N"
    }else{
      TD = 45
    }
    
    if(tempchar == "e"){
      BVtype = "original"
    }else if(tempchar == "w"){
      BVtype ="weighted"
    }else if(tempchar == "c"){
      BVtype = "std"
    }else if(tempchar == "o"){
      BVtype = "oweighted"
    }else{
      BVtype = "old"
    }
    
    if(i==1 & r==1){
      output<- as.data.table(cbind(proname,nGeneration,BVtype,ac,relM,TD,ggst,ggstg,genicEff,geneticEff,gg0,gg5,gg10,gg15,accallcand,accfam,variance,genicVariance,genicNe,geneticNe))
    }else{
      output = rbind(output,as.data.table(cbind(proname,nGeneration,BVtype,ac,relM,TD,ggst,ggstg,genicEff,geneticEff,gg0,gg5,gg10,gg15,accallcand,accfam,variance,genicVariance,genicNe,geneticNe)))
    }
  }
}
colnames(output) = c("pname","nGeneration","BV","tac","relMat","TargetD","ggst","ggstg","genicEff","geneticEff","gg0","gg5","gg10","gg15","acforallcand","acforfam","Va","genicVa","genicNe","geneticNe")


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
      
      # x[pname == proname[i] & nGeneration == (g-1),inbse:=sd(x[pname == proname[i] & nGeneration == (g-1),inb])]
      # x[pname == proname[i] & nGeneration == (g-1),inb:=mean(x[pname == proname[i] & nGeneration == (g-1),inb])]
      # 
      # x[pname == proname[i] & nGeneration == (g-1),acforcandse:=sd(x[pname == proname[i] & nGeneration == (g-1),acforcand])]
      # x[pname == proname[i] & nGeneration == (g-1),acforcand:=mean(x[pname == proname[i] & nGeneration == (g-1),acforcand])]
      
      x[pname == proname[i] & nGeneration == (g-1),acforallcandse:=sd(x[pname == proname[i] & nGeneration == (g-1),acforallcand])]
      x[pname == proname[i] & nGeneration == (g-1),acforallcand:=mean(x[pname == proname[i] & nGeneration == (g-1),acforallcand])]
      
      x[pname == proname[i] & nGeneration == (g-1),acforfamse:=sd(x[pname == proname[i] & nGeneration == (g-1),acforfam])]
      x[pname == proname[i] & nGeneration == (g-1),acforfam:=mean(x[pname == proname[i] & nGeneration == (g-1),acforfam])]
      
      # if(any(colnames(x)=="acimp")){
      #   x[pname == proname[i] & nGeneration == (g-1),acimpse:=sd(x[pname == proname[i] & nGeneration == (g-1),acimp])]
      #   x[pname == proname[i] & nGeneration == (g-1),acimp:=mean(x[pname == proname[i] & nGeneration == (g-1),acimp])]
      # }
      
      x[pname == proname[i] & nGeneration == (g-1),Vase:=sd(x[pname == proname[i] & nGeneration == (g-1),Va])]
      x[pname == proname[i] & nGeneration == (g-1),Va:=mean(x[pname == proname[i] & nGeneration == (g-1),Va])]
      
      x[pname == proname[i] & nGeneration == (g-1),genicVase:=sd(x[pname == proname[i] & nGeneration == (g-1),genicVa])]
      x[pname == proname[i] & nGeneration == (g-1),genicVa:=mean(x[pname == proname[i] & nGeneration == (g-1),genicVa])]
      
      # x[pname == proname[i] & nGeneration == (g-1),Nese:=sd(x[pname == proname[i] & nGeneration == (g-1),Ne])]
      # x[pname == proname[i] & nGeneration == (g-1),Ne:=mean(x[pname == proname[i] & nGeneration == (g-1),Ne])]
      
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

# library(showtext)
# showtext_auto(enable = TRUE)


#======================================================
pv = "genicEff"
pvse = "genicEffse"
zt = zt1[nGeneration==20&TargetD==45&BV == "old",]
# zt$tac <- c(0, 1, 0.1, 0.3, 0.5, 0.7, 0.9)[match(zt$tac, c(0, 1, 2, 3, 4, 5, 6))]
zt$TargetD = as.character(zt$TargetD)
zt$value = zt[,..pv]
zt$valuese = zt[,..pvse]
pedValue = zt1[nGeneration==20&relMat=="A"&tac==0&TargetD==45,..pv] %>% as.numeric()
P<- ggplot(data = zt,aes(x=tac,y=value,group = relMat, colour = relMat))+
  geom_point()+
  geom_line()+
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  xlab("The accuracy within family")+
  ylab("Conversion efficiency")+
  theme_zg()+theme(legend.position.inside=c(0.85, 0.12),legend.title=element_blank())+
  #geom_ribbon(aes(ymin = pgvr -pgvrse, ymax = pgvr+pgvrse), alpha = 0.4)+
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
  geom_errorbar(aes(ymin=value-valuese,
                    ymax=value+valuese),
                width=0.05,alpha = 0.5)+labs(color="Mat")+scale_color_aaas()
# P+geom_hline(yintercept = zt[denSNP=="1250"&nRef==200,gg0])
P
P+geom_hline(yintercept = pedValue)
#======================================================


pv = "genicEff"
pvse = "genicEffse"
zt = rbind(zt1[nGeneration==20&relMat=="G"&tac!=0&BV == "old",],zt1[nGeneration==20&relMat=="A"&tac==0&BV == "old",],zt1[nGeneration==20&relMat=="N"&BV == "old",])
# zt$tac <- c(0, 1, 0.1, 0.3, 0.5, 0.7, 0.9)[match(zt$tac, c(0, 1, 2, 3, 4, 5, 6))]
zt$TargetD = as.character(zt$TargetD)
zt$value = zt[,..pv]
zt$valuese = zt[,..pvse]
P<- ggplot(data = zt,aes(x=tac,y=value,group = TargetD, colour = TargetD))+
  geom_point()+
  geom_line()+
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  xlab("The accuracy within family")+
  ylab("Conversion efficiency")+
  theme_zg()+theme(legend.position.inside=c(0.85, 0.12),legend.title=element_blank())+
  #geom_ribbon(aes(ymin = pgvr -pgvrse, ymax = pgvr+pgvrse), alpha = 0.4)+
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
  geom_errorbar(aes(ymin=value-valuese,
                    ymax=value+valuese),
                width=0.05,alpha = 0.5)+labs(color="Mat")+scale_color_aaas()
# P+geom_hline(yintercept = zt[denSNP=="1250"&nRef==200,gg0])
P


#======================================================


pv = "genicEff"
pvse = "genicEffse"

pedValue = zt1[nGeneration==20&relMat=="A"&tac==0&TargetD==45,..pv] %>% as.numeric()

zt = rbind(zt1[nGeneration==20&BV%in%c("original","weighted","std","oweighted"),],zt1[nGeneration==20&BV=="old"&relMat == "G"&TargetD==45&tac!=0,])

# zt$tac <- c(0, 1, 0.1, 0.3, 0.5, 0.7, 0.9)[match(zt$tac, c(0, 1, 2, 3, 4, 5, 6))]
zt$TargetD = as.character(zt$TargetD)
zt$value = zt[,..pv]
zt$valuese = zt[,..pvse]

zt = zt[BV%in%c("std","old","oweighted"),]

zt[BV=="old",BV:="original"]
zt[BV == "oweighted",BV:="weighted"]
P<- ggplot(data = zt,aes(x=tac,y=value,group = BV, colour = BV))+
  geom_point()+
  geom_line()+
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  xlab("The accuracy within family")+
  ylab("Conversion efficiency")+
  theme_zg()+theme(legend.position.inside=c(0.85, 0.12),legend.title=element_blank())+
  #geom_ribbon(aes(ymin = pgvr -pgvrse, ymax = pgvr+pgvrse), alpha = 0.4)+
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
  geom_errorbar(aes(ymin=value-valuese,
                    ymax=value+valuese),
                width=0.05,alpha = 0.5)+labs(color="GEBV")+scale_color_aaas()
P

P+geom_hline(yintercept = pedValue)



#Additional file 6, Figure S17
zt = zt[BV=="original",]

P<- ggplot(data = zt,aes(x=tac,y=value,group = BV, colour = BV))+
  geom_point()+
  geom_line()+
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  xlab("Within-family accuracy")+
  ylab("Conversion efficiency")+
  theme_zg()+theme(legend.position.inside=c(0.85, 0.12),legend.title=element_blank())+
  #geom_ribbon(aes(ymin = pgvr -pgvrse, ymax = pgvr+pgvrse), alpha = 0.4)+
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
  geom_errorbar(aes(ymin=value-valuese,
                    ymax=value+valuese),
                width=0.05,alpha = 0.5)+labs(color="GEBV")+scale_color_aaas()
P

P = P+theme(legend.position = "none")+geom_hline(yintercept = pedValue)

P
ggsave("FigrueS17.pdf", P , width = 15, height = 5, dpi = 300)
#
output = fread("alldt.csv",sep = ",")
zt = rbind(output[nGeneration==20&relMat=="G"&tac!=0,],output[nGeneration==20&relMat=="A"&tac==0,],output[nGeneration==20&relMat=="N",])
chkt = t.test(x = zt[TargetD ==45 & tac==0,genicEff],y = zt[TargetD ==45 & tac==0.7,genicEff])
chkt
chkt$p.value
