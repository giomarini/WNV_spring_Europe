
rm(list = ls())
graphics.off()

### required libraries

library(tidyverse)
library(gridExtra)
library(grid)
library(MASS)
library(pscl)
library(sf)
library(RColorBrewer)

to_plot=T #if TRUE, plots are saved in the specified output file

#### reading data #####
### the dataset was made merging the ECDC datasets
### temperature data was averaged from MODIS daily time series
### Temperature covariates:
### AVG_2003_2010: spring average for 2003-2010
### STD_ANOMALY: standardized spring anomaly wrt 2003-2010

WNV_data=read.table(file="dataset_for_analysis.txt",header=T)

WNV_data$WNV_before=as.factor(WNV_data$WNV_before)

nrow(WNV_data)
mean(WNV_data$CASES)
length(which(WNV_data$CASES==0))
length(which(WNV_data$CASES!=0))

round(mean(WNV_data$CASES),1)
round(mean(WNV_data$CASES[which(WNV_data$CASES>0)]),1)
max(WNV_data$CASES)


### Panels for figure 1 are saved in separated files which can be combined 
### with any graphic program (for instance Inkscape)
if(to_plot){
  pdf("Figure 1a.pdf",width=9,height=5)
  par(cex.lab=1.5,cex.axis=1.5,mar=c(5,7,2,2))
  
  hist(WNV_data$CASES,main="",xlab=expression("H"["y,i"]),
       ylab="",axes=F,breaks=300,freq=T,ylim=c(0,450))
  axis(1,at=seq(0,200,50))
  axis(2,las=2,at=seq(0,450,50))
  mtext(text="Count",side=2,line=4,cex=1.5)
  mtext(side=3,"a",adj=0,cex=2)
  par(fig = c(0.3,1, 0.3, 1), new = T)
  hist(WNV_data$CASES[which(WNV_data$CASES>=20)],main="",xlab="",
       ylab="",axes=F,breaks=300,freq=T,ylim=c(0,3))
  axis(1,at=seq(20,200,30))
  axis(2,las=2,at=seq(0,3,1))
  mtext(text="Count",side=2,line=2,cex=1.5)
  mtext(text=expression("H"["y,i"]),side=1,line=3,cex=1.5)
  
  dev.off()
}

if(to_plot){
  pdf("Figure 1b.pdf",width=9,height=5)
  par(cex.lab=1.5,cex.axis=1.5,mar=c(5,6,2,1))
  
  for_plot=WNV_data %>% 
    group_by(YEAR) %>%
    summarize(total_cases=sum(CASES))
  
  a=barplot(for_plot$total_cases,axes=F,ylim=c(0,1600),ylab="")
  axis(1,at=a,labels=2011:2019)
  axis(2,at=seq(0,1600,400),las=2)
  mtext(side=3,"b",adj=0,cex=2)
  mtext(text="Total WNV cases",side=2,line=4,cex=1.5)
  
  dev.off()
}

### Fig1c: Europe map with WNV total cases
### the shapefiles are publicly available from Eurostat
### https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/administrative-units-statistical-units/nuts
europe_nuts=st_read("Europe_nuts/NUTS_RG_20M_2016_3857.shp")

total_by_nuts=WNV_data %>%
  group_by(NUTS_ID) %>%
  summarise(total=sum(CASES,na.rm=T))

breaks_colors=c(1,2,10,20,50,100)

colors_map=rev(brewer.pal(length(breaks_colors),"RdYlBu"))

#removing UK, Ireland, Cyprus, Turkey
map_plot=europe_nuts[-which(europe_nuts$CNTR_CODE=="UK" | 
                              europe_nuts$CNTR_CODE=="TR" |
                              europe_nuts$CNTR_CODE=="IE" |
                              europe_nuts$CNTR_CODE=="CY"),]
map_plot=map_plot[,1]
map_plot$COLOR="white"
for(i in 1:nrow(total_by_nuts)){
  row_map=which(as.character(map_plot$NUTS_ID)==as.character(total_by_nuts$NUTS_ID[i]))
  which_color=findInterval(total_by_nuts$total[i],breaks_colors)
  map_plot$COLOR[row_map]=colors_map[which_color]
}

tmp=map_plot[,3]
xmin=-300000
xmax=2500000
ymin=4100000
ymax=7500000

if(to_plot)
  pdf("Figure 1c.pdf",width=10,height=5)
plot(tmp,xlim=c(xmin,xmax),ylim=c(ymin,ymax),main="",col=tmp$COLOR)
legenda=c(1,"2-9","10-19","20-49","50-99",">100")
legend(x="topleft",pch=19,bty="n",col=colors_map,
       legend=legenda,title="Total WNV cases",ncol=2,cex=1.45)
mtext(side=3,"c",adj=0,cex=2)
if(to_plot)
  dev.off()

#### all models computation ####

covariates <- c("","+AVG_2003_2010","+STD_ANOMALY","+WNV_before")

#formulas for count model
formulas=c("CASES~1")
for(l in 1:length(covariates)){
  ciao=combn(covariates, l,simplify=F)
  for(i in 1:length(ciao)){
    a="CASES~1"
    for(j in 1:length(ciao[[i]]))
      a=paste0(a,ciao[[i]][j])
    formulas=rbind(formulas,a)
  }
}
formulas=unique(formulas)

#formulas for binomial model
formulas_zero_model=c("1")
for(l in 1:length(covariates)){
  ciao=combn(covariates, l,simplify=F)
  for(i in 1:length(ciao)){
    a="1"
    for(j in 1:length(ciao[[i]]))
      a=paste0(a,ciao[[i]][j])
    formulas_zero_model=rbind(formulas_zero_model,a)
  }
}
formulas_zero_model=unique(formulas_zero_model)

#formulas for ZI model
formulas_ZI=c()
for(i in 1:nrow(formulas))
  for(j in 1:nrow(formulas_zero_model))
    formulas_ZI=rbind(formulas_ZI,paste0(formulas[i],"|",formulas_zero_model[j]))
formulas_ZI=unique(formulas_ZI)

AICtab <- matrix(nrow=nrow(formulas_ZI),ncol=3)
index_tab=1

#ZINB models
for(i in 1:nrow(formulas_ZI)){
  mod <- try(zeroinfl(formula(formulas_ZI[i,]),data=WNV_data,dist="negbin"),silent=T)
  if(is.list(mod)){
    dummy <- c("ZINB", formulas_ZI[i,],AIC(mod))
    AICtab[index_tab,]=dummy
    index_tab=index_tab+1
  }
}

AICdf=as.data.frame(AICtab)
names(AICdf)=c("model","formula","AIC")
AICdf$AIC=as.numeric(as.character(AICdf$AIC))
AICdf=AICdf %>% arrange(AIC) %>% 
  mutate(deltaAIC=AIC-min(AIC,na.rm=T))
AICdf$AIC=round(AICdf$AIC,2)
AICdf$deltaAIC=round(AICdf$deltaAIC,2)


#### finding best model ####

frm=formula(as.character(AICdf$formula[1]))
best_model1=zeroinfl(frm,data=WNV_data,dist="negbin")
summary(best_model1)

frm=formula(as.character(AICdf$formula[2]))
best_model2=zeroinfl(frm,data=WNV_data,dist="negbin")
summary(best_model2)

AIC(best_model1,best_model2)

# so we select best_model2 as best model, as it has all significant coefficients and a slighlty higher AIC
# We now perform some model diagnostics:
# residuals VS explanatory variables
res=residuals(best_model2,type="pearson")
par(mfrow=c(1,3))
plot(WNV_data$AVG_2003_2010,res,xlab="AVG_2003_2010",ylab="Pearson residuals")
plot(WNV_data$STD_ANOMALY,res,xlab="STD_ANOMALY",ylab="Pearson residuals")
boxplot(res~WNV_data$WNV_before,xlab="WNV_before",ylab="Pearson residuals")

# fitted values
fits=fitted(best_model2)
plot(fits,res,xlab="Fitted values",ylab="Pearson residuals")
plot(WNV_data$CASES,fits,xlab=expression("H"["i"]*"(y)"),ylab="Fitted values")
plot(WNV_data$CASES,res,xlab=expression("H"["i"]*"(y)"),ylab="Pearson residuals")


#### Conditional predictions #####

# plotting conditional to T_2003_2010, so STD_ANOMALY=0
tmin=12
tmax=22
newdb <- expand.grid(AVG_2003_2010=
                       seq(tmin,tmax,length=100),
                     STD_ANOMALY = 0,
                     WNV_before = factor(c(0,1,"NR")))

betazero <- best_model2$coefficients$zero
Xzero    <- model.matrix(~1 +STD_ANOMALY+WNV_before,data = newdb)

newdb$pi_zero   <- Xzero %*% betazero
newdb$prob_zero <- exp(newdb$pi_zero)/(1+exp(newdb$pi_zero)) #prob of having zero in the binomial process

betacount <- best_model2$coefficients$count
Xcount    <- model.matrix(~1+AVG_2003_2010+STD_ANOMALY+WNV_before,data = newdb)
newdb$mu   <- exp(Xcount %*% betacount)
newdb$mean <- (1-newdb$prob_zero)*newdb$mu

newdb$p0count <- dnbinom(0,size = best_model2$theta, mu = newdb$mu) #probability of zero in the count process

newdb$prob_zero_both=newdb$prob_zero+(1-newdb$prob_zero)*newdb$p0count

p1 <- newdb %>% transform(WNV_before = case_when(WNV_before== "1"~"Cases in previous year",
                                                 WNV_before == "0"~"No Cases in previous year",
                                                 WNV_before== "NR"~"Cases in previous year not reported" ))%>%
  ggplot(.,aes(x=AVG_2003_2010,y=mean,col=WNV_before))+geom_line(lwd=2)+theme_bw()+
  ylab(expression("E(H"["y,i"]*")")) + xlab(expression(hat("T")))+ labs(col="WNV_BEFORE")+
  scale_x_continuous(breaks = seq(tmin,tmax,2))+
  labs(tag = "a")+ylim(0,25)+
  theme(legend.position = c(0.35, 0.85),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15,face="bold"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.tag=element_text(size=20))
  

p3 <- newdb %>% transform(WNV_before = case_when(WNV_before== "1"~"Cases in previous year",
                                                 WNV_before == "0"~"No Cases in previous year",
                                                 WNV_before== "NR"~"Cases in previous year not reported" ))%>%
  ggplot(.,aes(x=AVG_2003_2010,y=prob_zero_both,col=WNV_before))+geom_line(lwd=2)+theme_bw()+
  ylab(expression("P(H"["y,i"]*"=0)")) + xlab(expression(hat("T")))+
  scale_x_continuous(breaks = seq(tmin,tmax,2))+
  theme(legend.position = "none")+ylim(0,1)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.tag=element_text(size=20))+
  labs(tag = "c")


# plotting conditional to STD_ANOMALY, so T_2003_2010=mean(T_2003_2010)
anom_min=-5
anom_max=5
newdb <- expand.grid(AVG_2003_2010=mean(WNV_data$AVG_2003_2010),
                     STD_ANOMALY = seq(anom_min,anom_max,length=100),
                     WNV_before = factor(c(0,1,"NR")))

betazero <- best_model2$coefficients$zero
Xzero    <- model.matrix(~1 +STD_ANOMALY+WNV_before,data = newdb)

newdb$pi_zero   <- Xzero %*% betazero
newdb$prob_zero <- exp(newdb$pi_zero)/(1+exp(newdb$pi_zero))

betacount <- best_model2$coefficients$count
Xcount    <- model.matrix(~1+AVG_2003_2010+STD_ANOMALY+WNV_before,data = newdb)
newdb$mu   <- exp(Xcount %*% betacount)
newdb$mean <- (1-newdb$prob_zero)*newdb$mu

p2 <- newdb %>% transform(WNV_before = case_when(WNV_before== "1"~"Cases in previous year",
                                                 WNV_before == "0"~"No Cases in previous year",
                                                 WNV_before== "NR"~"Cases in previous year not reported" ))%>%
  ggplot(.,aes(x=STD_ANOMALY,y=mean,col=WNV_before))+geom_line(lwd=2)+theme_bw()+
  ylab(expression("E(H"["y,i"]*")")) + 
  scale_x_continuous(breaks = seq(anom_min,anom_max,2))+
  theme(legend.position = "none")+ylim(0,25)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.tag=element_text(size=20))+
  labs(tag = "b")

newdb$p0count <- dnbinom(0,size = best_model2$theta, mu = newdb$mu) #probability of zero in the count process

newdb$prob_zero_both=newdb$prob_zero+(1-newdb$prob_zero)*newdb$p0count

p4 <- newdb %>% transform(WNV_before = case_when(WNV_before== "1"~"Cases in previous year",
                                                 WNV_before == "0"~"No Cases in previous year",
                                                 WNV_before== "NR"~"Cases in previous year not reported" ))%>%
  ggplot(.,aes(x=STD_ANOMALY,y=prob_zero_both,col=WNV_before))+geom_line(lwd=2)+theme_bw()+
  ylab(expression("P(H"["y,i"]*"=0)")) + 
  scale_x_continuous(breaks = seq(anom_min,anom_max,2))+
  theme(legend.position = "none")+ylim(0,1)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.tag=element_text(size=20))+
  labs(tag = "d")

grid.arrange(p1,p2,p3,p4,ncol=2)

if(to_plot){
  pdf("Figure 2.pdf",width=14,height=10)
  grid.arrange(p1,p2,p3,p4,ncol=2)
  dev.off()
}

#### using positive values only ####
### here we use the original ECDC combined datasets (so no zeros) to check whether our results hold
WNV_data_nozeros=WNV_data[-which(WNV_data$CASES==0),]

par(mfrow=c(1,1))
hist(log(WNV_data_nozeros$CASES))

# we will compute linear models with log-transformation of Y
covariates <- c("","+AVG_2003_2010","+STD_ANOMALY")
formulas=c("log(CASES)~1")
for(l in 1:length(covariates)){
  ciao=combn(covariates, l,simplify=F)
  for(i in 1:length(ciao)){
    a="log(CASES)~1"
    for(j in 1:length(ciao[[i]]))
      a=paste0(a,ciao[[i]][j])
    formulas=rbind(formulas,a)
  }
}
formulas=unique(formulas)

AICtab <- matrix(nrow=nrow(formulas),ncol=3)
index_tab=1

#LM
for(i in 1:nrow(formulas)){
  mod <- lm(formula(formulas[i,]),data=WNV_data_nozeros)
  dummy <- c("LM", formulas[i,],AIC(mod))
  # AICtab <- bind_rows(AICtab,dummy)
  AICtab[index_tab,]=dummy
  index_tab=index_tab+1
}

AICdf_nozeros=as.data.frame(AICtab)
names(AICdf_nozeros)=c("model","formula","AIC")
AICdf_nozeros$AIC=as.numeric(as.character(AICdf_nozeros$AIC))
AICdf_nozeros=AICdf_nozeros %>% arrange(AIC) %>% 
  mutate(deltaAIC=AIC-min(AIC,na.rm=T))
AICdf_nozeros$AIC=round(AICdf_nozeros$AIC,2)
AICdf_nozeros$deltaAIC=round(AICdf_nozeros$deltaAIC,2)

frm=formula(as.character(AICdf_nozeros$formula[1]))
best_model_nozeros1=lm(frm,data=WNV_data_nozeros)
summary(best_model_nozeros1)

# so the positive relationship with spring temperature (both covariates) is confirmed
# we can check the residuals
res=residuals(best_model_nozeros1,type="pearson")
par(mfrow=c(2,2))
hist(res)
qqnorm(res)
plot(WNV_data_nozeros$AVG_2003_2010,res,xlab="AVG_2003_2010",ylab="Pearson residuals")
plot(WNV_data_nozeros$STD_ANOMALY,res,xlab="STD_ANOMALY",ylab="Pearson residuals")
