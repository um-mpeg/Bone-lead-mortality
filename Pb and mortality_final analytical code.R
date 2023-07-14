library(survival)
library(epiDisplay)
library(Hmisc)
#### load the survey package
library(survey)
library(dplyr)
library(ggplot2)
library("survminer")
options(survey.lonely.psu="adjust")

setwd("C:\\Users\\envie\\Dropbox (University of Michigan)\\My Documents\\NHANES\\III\\updated mortality\\reanalysis")
setwd("C:\\Users\\sungkyun\\Dropbox (University of Michigan)\\My Documents\\NHANES\\III\\updated mortality\\reanalysis")

###########################################################################################
###################### DO NOT RUN #################################

library(haven)
pbmort2<-read_sas("pbmort2.sas7bdat")
names(pbmort2)
names(pbmort2)<-tolower(names(pbmort2))

### 11/23/2022, revision ###

### add age in month
agemon<-read_sas("C:\\Users\\envie\\Dropbox (University of Michigan)\\My Documents\\NHANES\\III\\agemon.sas7bdat")
names(agemon)
names(agemon)<-tolower(names(agemon))

pbmort2<-merge(pbmort2, agemon, by="seqn", all.x=T)
summ(pbmort2)

pbmort2[1:10,c("age","hsaitmor")]
summary(pbmort2$age)
summary(pbmort2$hsaitmor)
# for top-coded age, (90+3)*12=1116 were recoded 

# Also, recode top-coded age to 93
pbmort2$age<-ifelse(pbmort2$age==90, 93, pbmort2$age)
summ(pbmort2$age)

pbmort2$age.die<-pbmort2$hsaitmor+pbmort2$pmon_mec

summ(pbmort2$age)
summ(pbmort2$hsaitmor)

## create tertiles
pbmort2$bpb3<-cut2(pbmort2$bpb, g=3)
pbmort2$tibm3<-cut2(pbmort2$tibm, g=3)
pbmort2$patm3<-cut2(pbmort2$patm, g=3)

## create quartiles
pbmort2$bpb4<-cut2(pbmort2$bpb, g=4)
pbmort2$tibm4<-cut2(pbmort2$tibm, g=4)
pbmort2$patm4<-cut2(pbmort2$patm, g=4)
names(pbmort2)

### Standardize: comparing 10th vs 90th
quantile(pbmort2$bpb, c(0.1,0.9))
quantile(pbmort2$logbpb, c(0.1,0.9))
pbmort2$bpb.q91<-pbmort2$bpb/(quantile(pbmort2$bpb, 0.9)-quantile(pbmort2$bpb, 0.1))
pbmort2$logbpb.q91<-pbmort2$logbpb/(quantile(pbmort2$logbpb, 0.9)-quantile(pbmort2$logbpb, 0.1))
pbmort2$tibm.q91<-pbmort2$tibm/(quantile(pbmort2$tibm, 0.9)-quantile(pbmort2$tibm, 0.1))
pbmort2$patm.q91<-pbmort2$patm/(quantile(pbmort2$patm, 0.9)-quantile(pbmort2$patm, 0.1))

### log-transform tibm and patm
summ(tibm)
summ(patm)

tibm.flag<-ifelse(tibm<=0, 1, 0)
table(tibm.flag)

pbmort2$logtib<-log(tibm+8.815)
pbmort2$logpat<-log(patm)
summ(pbmort2$logtib)
summ(pbmort2$logpat)
hist(pbmort2$logtib,n=20)
hist(pbmort2$logpat,n=20)

pbmort2$logtib.q91<-pbmort2$logtib/(quantile(pbmort2$logtib, 0.9)-quantile(pbmort2$logtib, 0.1))
pbmort2$logpat.q91<-pbmort2$logpat/(quantile(pbmort2$logpat, 0.9)-quantile(pbmort2$logpat, 0.1))



#### 2/28/2023
### Add updated bone lead
bonepb<-read_sas("C:\\Users\\sungkyun\\Dropbox (University of Michigan)\\My Documents\\NHANES\\bone lead\\nhanes_iii_pred_bonepb.sas7bdat")
names(bonepb)
names(bonepb)<-tolower(names(bonepb))

pbmort3<-merge(pbmort2, bonepb, by="seqn", all.x=T)
summ(pbmort3$tibm)
summ(pbmort3$pred_tibx)
summ(pbmort3$patm)
summ(pbmort3$pred_patx)

tibx.flag<-ifelse(pbmort3$pred_tibx<=0, 1, 0)
table(tibx.flag)
tab1(tibx.flag)

pbmort3$logtib<-log(pbmort3$pred_tibx+9.607)
pbmort3$logpat<-log(pbmort3$pred_patx)
summ(pbmort3$logtib)
summ(pbmort3$logpat)
hist(pbmort3$logtib,n=20)
hist(pbmort3$logpat,n=20)

pbmort3$logtib.q91<-pbmort3$logtib/(quantile(pbmort3$logtib, 0.9)-quantile(pbmort3$logtib, 0.1))
pbmort3$logpat.q91<-pbmort3$logpat/(quantile(pbmort3$logpat, 0.9)-quantile(pbmort3$logpat, 0.1))

### check the proportion of pred_tibx<0
pbmort3$tibm.neg<-ifelse(pbmort3$pred_tibx<=0, 1, 0)
tab1(pbmort3$tibm.neg)
#n=1207, 10.4%

### bpb<1
tab1(pbmort3$bpb.lod)
#n=970, 8.3%

tabpct(pbmort3$bpb.lod, pbmort3$tibm.neg)

pbmort3$tibm.imp<-ifelse(pbmort3$pred_tibx<=0, 0.001, pbmort3$pred_tibx)
summ(pbmort3$tibm.imp)

pbmort3$logtib.imp<-log(pbmort3$tibm.imp)
pbmort3$logtib.imp.q91<-pbmort3$logtib.imp/(quantile(pbmort3$logtib.imp, 0.9)-quantile(pbmort3$logtib.imp, 0.1))

summ(pbmort3$logtib.q91)
summ(pbmort3$logtib.imp.q91)


### compare tertiles

## create tertiles
pbmort3$tibm3<-cut2(pbmort3$pred_tibx, g=3)
pbmort3$patm3<-cut2(pbmort3$pred_patx, g=3)

## create quartiles
pbmort3$tibm4<-cut2(pbmort3$pred_tibx, g=4)
pbmort3$patm4<-cut2(pbmort3$pred_patx, g=4)

save(pbmort3, file="pbmort3.rda")

# export csv
#write.csv(pbmort3, "pbmort3.csv")

#### DO NOT RUN

###########################################################################################

#############################################################################################################################################
############################ START FROM HERE #################################################
#############################################################################################################################################

load("pbmort3.rda")
names(pbmort3)

attach(pbmort3)

bpdsn<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=pbmort3, nest=T)

### incident rate
summ(pmon_mec)
summ(pyr_mec)
summary(pmon_mec)
summary(pmon_mec/12)
sum(pmon_mec)
sum(pyr_mec)
sum(wt_mh)
svymean(~pmon_mec, bpdsn)
confint(svymean(~pmon_mec, bpdsn))

table(m_total)
table(m_cvd)
table(m_heart)
table(m_cancer)

table(m_total,raceth)
table(m_cvd,raceth)
table(m_heart,raceth)

sum(m_total)/sum(pmon_mec)*100000
sum(m_cvd)/sum(pmon_mec)*100000
sum(m_heart)/sum(pmon_mec)*100000
sum(m_cancer)/sum(pmon_mec)*100000

ci.poisson(sum(m_total), sum(pmon_mec)/(100000*12), alpha=.05)
ci.poisson(sum(m_total*wt_mh), sum(pmon_mec*wt_mh)/(100000*12), alpha=.05)

ci.poisson(sum(m_cvd), sum(pmon_mec)/(100000*12), alpha=.05)
ci.poisson(sum(m_cvd*wt_mh), sum(pmon_mec*wt_mh)/(100000*12), alpha=.05)

ci.poisson(sum(m_heart), sum(pmon_mec)/(100000*12), alpha=.05)
ci.poisson(sum(m_heart*wt_mh), sum(pmon_mec*wt_mh)/(100000*12), alpha=.05)

ci.poisson(sum(m_cancer), sum(pmon_mec)/(100000*12), alpha=.05)
ci.poisson(sum(m_cancer*wt_mh), sum(pmon_mec*wt_mh)/(100000*12), alpha=.05)


### For Table 1: Pop Characteristics
### done in SAS

## check mortality cases at top-coded age
pbmort.age90<-subset(pbmort3, age==93)
table(pbmort.age90$m_total)
table(pbmort.age90$m_cvd)
table(pbmort.age90$m_heart)
table(pbmort.age90$m_cancer)


### correlation among lead measures
library(corrplot)

## calculate p-value matrix
names(pbmort3)
x <- pbmort3[,c("bpb","pred_tibx","pred_patx")]
logx <- pbmort3[,c("logbpb","logtib","logpat")]

cor(x)
cor(logx)

p<- cor.mtest(x, conf.level = .95)
p
p$p

corrx = cor(x) #get the correlations
class(cor(x))
corrplot(corrx, method = "color")

corrplot(corrx, method = "color", 
         type = "upper", order = "hclust")

library(jtools)
svycor(~bpb+logbpb+tibm+logtib+patm+logpat, bpdsn)
svycor(~bpb+tibm+patm, bpdsn)
svycor(~logbpb+logtib+logpat, bpdsn)

corrx = svycor(~bpb+tibm+patm, bpdsn)
attributes(corrx)
class(corrx)
corrx<-data.frame(corrx$cors)

corrlogx = svycor(~logbpb+logtib+logpat, bpdsn)

corrplot(corrx, method = "color")

corrplot(corrx, method = "color", 
         type = "upper", order = "hclust")

corrplot(corrlogx, method = "color", 
         type = "upper", order = "hclust")


### BLL by race and sex

svyby(~bpb, ~raceth, bpdsn, svymean)
svyby(~pred_tibx, ~raceth, bpdsn, svymean)
svyby(~pred_patx, ~raceth, bpdsn, svymean)

svyby(~bpb, ~raceth, bpdsn, svyquantile, quantiles=c(0.25,0.5,0.75))
svyby(~pred_tibx, ~raceth, bpdsn, svyquantile, quantiles=c(0.25,0.5,0.75))
svyby(~pred_patx, ~raceth, bpdsn, svyquantile, quantiles=c(0.25,0.5,0.75))

svyby(~bpb, ~sex, bpdsn, svyquantile, quantiles=c(0.25,0.5,0.75))

bpb.sex<-svyglm(logbpb~factor(sex), bpdsn)
summary(bpb.sex)
tib.sex<-svyglm(logtib~factor(sex), bpdsn)
summary(tib.sex)
pat.sex<-svyglm(logpat~factor(sex), bpdsn)
summary(pat.sex)

bpb.race<-svyglm(logbpb~factor(raceth), bpdsn)
summary(bpb.race)
tib.race<-svyglm(logtib~factor(raceth), bpdsn)
summary(tib.race)
pat.race<-svyglm(logpat~factor(raceth), bpdsn)
summary(pat.race)


######################################################################
### TABLE 2 Survival Analysis ###
######################################################################

#########################################
### use age.die instead of PERMTH_EXM ###
#########################################

########################### ALL-CAUSE MORTALITY ################################
## unadjusted
svycox.bpb.un<-svycoxph(Surv(hsaitmor, age.die, m_total)~bpb.q91, bpdsn)
summary(svycox.bpb.un)
svycox.logbpb.un<-svycoxph(Surv(hsaitmor, age.die, m_total)~logbpb.q91, bpdsn)
summary(svycox.logbpb.un)
svycox.tibm.un<-svycoxph(Surv(hsaitmor, age.die, m_total)~tibm.q91, bpdsn)
summary(svycox.tibm.un)
svycox.patm.un<-svycoxph(Surv(hsaitmor, age.die, m_total)~patm.q91, bpdsn)
summary(svycox.patm.un)
## Unadjusted are too high, so do not report


## fully adjusted
svycox.logbpb.q91<-svycoxph(Surv(hsaitmor, age.die, m_total)~logbpb.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                            +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logbpb.q91)

svycox.logtib.q91<-svycoxph(Surv(hsaitmor, age.die, m_total)~logtib.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                          +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logtib.q91)

svycox.logpat.q91<-svycoxph(Surv(hsaitmor, age.die, m_total)~logpat.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                          +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logpat.q91)

test.prop.bpb<-cox.zph(svycox.logbpb.q91)
test.prop.bpb
test.prop.logtib<-cox.zph(svycox.logtib.q91)
test.prop.logtib
test.prop.logpat<-cox.zph(svycox.logpat.q91)
test.prop.logpat

########################### CVD MORTALITY ################################
## fully adjusted
svycox.logbpb.q91.cvd<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logbpb.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logbpb.q91.cvd)
attributes(summary(svycox.logbpb.q91.cvd))
attributes(svycox.logbpb.q91.cvd)
summary(svycox.logbpb.q91.cvd)$conf.int[1,1]
summary(svycox.logbpb.q91.cvd)$conf.int[1,3]
summary(svycox.logbpb.q91.cvd)$conf.int[1,4]

#svycox.tibm.q91.cvd<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~tibm.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
#                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
#summary(svycox.tibm.q91.cvd)

#svycox.patm.q91.cvd<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~patm.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
#                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
#summary(svycox.patm.q91.cvd)

svycox.logtib.q91.cvd<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logtib.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logtib.q91.cvd)

svycox.logpat.q91.cvd<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logpat.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logpat.q91.cvd)

test.prop.bpb.cvd<-cox.zph(svycox.logbpb.q91.cvd)
test.prop.bpb.cvd
test.prop.logtib.cvd<-cox.zph(svycox.logtib.q91.cvd)
test.prop.logtib.cvd
test.prop.logpat.cvd<-cox.zph(svycox.logpat.q91.cvd)
test.prop.logpat.cvd

########################### HEART MORTALITY ################################
## fully adjusted
svycox.logbpb.q91.heart<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logbpb.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                  +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logbpb.q91.heart)

#svycox.tibm.q91.heart<-svycoxph(Surv(hsaitmor, age.die, m_heart)~tibm.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
#                                +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
#summary(svycox.tibm.q91.heart)

#svycox.patm.q91.heart<-svycoxph(Surv(hsaitmor, age.die, m_heart)~patm.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
#                                +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
#summary(svycox.patm.q91.heart)

svycox.logtib.q91.heart<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logtib.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logtib.q91.heart)

svycox.logpat.q91.heart<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logpat.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logpat.q91.heart)

test.prop.bpb.heart<-cox.zph(svycox.logbpb.q91.heart)
test.prop.bpb.heart
test.prop.tibm.heart<-cox.zph(svycox.tibm.q91.heart)
test.prop.tibm.heart
test.prop.patm.heart<-cox.zph(svycox.patm.q91.heart)
test.prop.patm.heart

########################### CANCER MORTALITY ################################
## fully adjusted
svycox.logbpb.q91.cancer<-svycoxph(Surv(hsaitmor, age.die, m_cancer)~logbpb.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                   +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logbpb.q91.cancer)

#svycox.tibm.q91.cancer<-svycoxph(Surv(hsaitmor, age.die, m_cancer)~tibm.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
#                                 +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
#summary(svycox.tibm.q91.cancer)

#svycox.patm.q91.cancer<-svycoxph(Surv(hsaitmor, age.die, m_cancer)~patm.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
#                                 +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
#summary(svycox.patm.q91.cancer)

svycox.logtib.q91.cancer<-svycoxph(Surv(hsaitmor, age.die, m_cancer)~logtib.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                 +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logtib.q91.cancer)

svycox.logpat.q91.cancer<-svycoxph(Surv(hsaitmor, age.die, m_cancer)~logpat.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                 +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logpat.q91.cancer)

test.prop.bpb.cancer<-cox.zph(svycox.logbpb.q91.cancer)
test.prop.bpb.cancer
test.prop.tibm.cancer<-cox.zph(svycox.tibm.q91.cancer)
test.prop.tibm.cancer
test.prop.patm.cancer<-cox.zph(svycox.patm.q91.cancer)
test.prop.patm.cancer

################################
########## Tertiles ############
################################

tabpct(bpb3, m_total, percent="row")
tabpct(tibm3, m_total, percent="row")
tabpct(patm3, m_total, percent="row")

names(pbmort2)
pbmort2$agegroup<-cut2(age, c(40,60,70,80))
table(pbmort2$agegroup)


### All cause
svycox.bpb3<-svycoxph(Surv(hsaitmor, age.die, m_total)~bpb3+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                      +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.bpb3)
svycox.bpb3.t<-svycoxph(Surv(hsaitmor, age.die, m_total)~as.numeric(bpb3)+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                      +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.bpb3.t)

svycox.tibm3<-svycoxph(Surv(hsaitmor, age.die, m_total)~tibm3+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                       +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.tibm3)
svycox.tibm3.t<-svycoxph(Surv(hsaitmor, age.die, m_total)~as.numeric(tibm3)+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                       +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.tibm3.t)

svycox.patm3<-svycoxph(Surv(hsaitmor, age.die, m_total)~patm3+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                       +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.patm3)
svycox.patm3.t<-svycoxph(Surv(hsaitmor, age.die, m_total)~as.numeric(patm3)+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                       +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.patm3.t)

test.prop.bpb3<-cox.zph(svycox.bpb3)
test.prop.bpb3
test.prop.tibm3<-cox.zph(svycox.tibm3)
test.prop.tibm3
test.prop.patm3<-cox.zph(svycox.patm3)
test.prop.patm3

### cvd
svycox.bpb3<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~bpb3+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                      +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.bpb3)
svycox.bpb3.t<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~as.numeric(bpb3)+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                        +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.bpb3.t)

svycox.tibm3<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~tibm3+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                       +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.tibm3)
svycox.tibm3.t<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~as.numeric(tibm3)+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                         +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.tibm3.t)

svycox.patm3<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~patm3+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                       +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.patm3)
svycox.patm3.t<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~as.numeric(patm3)+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                         +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.patm3.t)

### heart
svycox.bpb3<-svycoxph(Surv(hsaitmor, age.die, m_heart)~bpb3+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                      +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.bpb3)
svycox.bpb3.t<-svycoxph(Surv(hsaitmor, age.die, m_heart)~as.numeric(bpb3)+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                        +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.bpb3.t)

svycox.tibm3<-svycoxph(Surv(hsaitmor, age.die, m_heart)~tibm3+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                       +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.tibm3)
svycox.tibm3.t<-svycoxph(Surv(hsaitmor, age.die, m_heart)~as.numeric(tibm3)+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                         +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.tibm3.t)

svycox.patm3<-svycoxph(Surv(hsaitmor, age.die, m_heart)~patm3+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                       +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.patm3)
svycox.patm3.t<-svycoxph(Surv(hsaitmor, age.die, m_heart)~as.numeric(patm3)+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                         +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.patm3.t)

### cancer
svycox.bpb3<-svycoxph(Surv(hsaitmor, age.die, m_cancer)~bpb3+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                      +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.bpb3)
svycox.bpb3.t<-svycoxph(Surv(hsaitmor, age.die, m_cancer)~as.numeric(bpb3)+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                        +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.bpb3.t)

svycox.tibm3<-svycoxph(Surv(hsaitmor, age.die, m_cancer)~tibm3+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                       +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.tibm3)
svycox.tibm3.t<-svycoxph(Surv(hsaitmor, age.die, m_cancer)~as.numeric(tibm3)+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                         +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.tibm3.t)

svycox.patm3<-svycoxph(Surv(hsaitmor, age.die, m_cancer)~patm3+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                       +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.patm3)
svycox.patm3.t<-svycoxph(Surv(hsaitmor, age.die, m_cancer)~as.numeric(patm3)+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                         +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.patm3.t)

#######################################
########## smoothing plots ############
#######################################
library(splines)

## center pred_tibx and pred_patx at their 10th
quantile(bpb, c(0.1))
quantile(pred_tibx, c(0.1))
quantile(pred_patx, c(0.1))
pbmort3$bpb10<-bpb-1
pbmort3$pred_tibx10<-pred_tibx+0.1239909
pbmort3$pred_patx10<-pred_patx-9.800195

save(pbmort3, file="pbmort3.rda")
bpdsn<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=pbmort3, nest=T)

svycox.bpb.sp.cvd<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~ns(bpb10,knots=c(quantile(bpb10, c(0.05, 0.5, 0.95))))+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                            +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.bpb.sp.cvd)

svycox.tibm.sp.cvd<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~ns(pred_tibx10, df=4)+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.tibm.sp.cvd)




svycox.tibm.sp.cvd<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~ns(pred_tibx,knots=c(quantile(pred_tibx, c(0.05, 0.5, 0.95))))+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.tibm.sp.cvd)

svycox.patm.sp.cvd<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~ns(pred_patx,knots=c(quantile(pred_patx, c(0.05, 0.5, 0.95))))+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.patm.sp.cvd)

### first calculate the reference value for OR -- predicted value for 10th percentile ###
quantile(pred_tibx, c(0.1))

#find index number for pred_tibx==-0.127205521, which is the 10th percentile
which(pred_tibx==-0.127205521) #1116
#order logtas
o<-order(pred_tibx)
#find index number within o for the reference value
which(o==1116) #173
o[173]



# Tibia lead
fit.tibm.cvd<-predict(svycox.tibm.sp.cvd, type="terms", se.fit=T)
fit.tibm.cvd$fit[1:10,]

p.tibm.cvd<-exp(fit.tibm.cvd$fit[,1])
l95.tibm.cvd<-exp(fit.tibm.cvd$fit[,1]-1.96*fit.tibm.cvd$se.fit[,1])
u95.tibm.cvd<-exp(fit.tibm.cvd$fit[,1]+1.96*fit.tibm.cvd$se.fit[,1])

o<-order(pred_tibx)
plot(pred_tibx[o], p.tibm.cvd[o], log="y", type="l", ylim=c(min(summary(l95.tibm.cvd)),max(summary(u95.tibm.cvd))),
     ylab="HR", xlab="Tibia lead (?g/g)")
lines(pred_tibx[o], l95.tibm.cvd[o], lty=2, col=2)
lines(pred_tibx[o], u95.tibm.cvd[o], lty=2, col=2)

## using ggplot
d.tibx<-data.frame(pred_tibx10, p.tibm.cvd, l95.tibm.cvd, u95.tibm.cvd)
summ(d.tibx)

ggplot(d.tibx,aes(x=pred_tibx10,y=p.tibm.cvd, group=1))+geom_line(col="black")+
  coord_cartesian(xlim = c(-8.5, 60), ylim = c(0, 4))+
  geom_hline(yintercept=1, color="red", linetype="dashed")+
  geom_ribbon(aes(ymin = l95.tibm.cvd, ymax = u95.tibm.cvd), alpha = 0.1)+
  ylab('HR')+xlab('Predicted tibia lead (?g/g)')




d.tibx<-data.frame(pred_tibx, p.tibm.cvd, l95.tibm.cvd, u95.tibm.cvd)
summ(d.tibx)


ggplot(d.tibx,aes(x=pred_tibx,y=p.tibm.cvd, group=1))+geom_line(col="black")+
  ylim(0, 4)+
  geom_hline(yintercept=1, color="red", linetype="dotted",size=1.5)+
  geom_ribbon(aes(ymin = l95.tibm.cvd, ymax = u95.tibm.cvd), alpha = 0.1)+
  ylab('HR')+xlab('Predicted tibia lead (?g/g)')

ggplot(d.tibx,aes(x=pred_tibx,y=p.tibm.cvd, group=1))+geom_line(col="black")+
  coord_cartesian(ylim = c(0, 4))+
  geom_hline(yintercept=1, color="red", linetype="dashed")+
  geom_ribbon(aes(ymin = l95.tibm.cvd, ymax = u95.tibm.cvd), alpha = 0.1)+
  ylab('HR')+xlab('Predicted tibia lead (?g/g)')



# Patella lead

quantile(pred_patx, c(0.1))


#find index number for pred_patx==9.800282, which is the 10th percentile
which(pred_patx==9.800282) #1116
#order logtas
o<-order(pred_tibx)
#find index number within o for the reference value
which(o==1116) #173
o[173]


fit.patm.cvd<-predict(svycox.patm.sp.cvd, type="terms", se.fit=T)
fit.patm.cvd$fit[1:10,]

p.patm.cvd<-exp(fit.patm.cvd$fit[,1])
l95.patm.cvd<-exp(fit.patm.cvd$fit[,1]-1.96*fit.patm.cvd$se.fit[,1])
u95.patm.cvd<-exp(fit.patm.cvd$fit[,1]+1.96*fit.patm.cvd$se.fit[,1])

o<-order(patm)
plot(patm[o], p.patm.cvd[o], log="y", type="l", ylim=c(min(summary(l95.patm.cvd)),max(summary(u95.patm.cvd))),
     ylab="HR", xlab="Patella lead (µg/g)")
lines(patm[o], l95.patm.cvd[o], lty=2, col=2)
lines(patm[o], u95.patm.cvd[o], lty=2, col=2)



#####################################################################
###### Figure 2: Effect modification by Sex and Race/ethnicity ######
#####################################################################

###########################
########## Sex ############
###########################

pbmort3.male<-subset(pbmort3, sex==1)
bpdsn.male<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=pbmort3.male, nest=T)

pbmort3.female<-subset(pbmort3, sex==2)
bpdsn.female<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=pbmort3.female, nest=T)

######### total #########
## male
total.logbpb.q91.m<-svycoxph(Surv(hsaitmor, age.die, m_total)~logbpb.q91+factor(income)+factor(raceth)+factor(obesity)
                            +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.male)
summary(total.logbpb.q91.m)

total.logtib.q91.m<-svycoxph(Surv(hsaitmor, age.die, m_total)~logtib.q91+factor(income)+factor(raceth)+factor(obesity)
                          +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.male)
summary(total.logtib.q91.m)

total.logpat.q91.m<-svycoxph(Surv(hsaitmor, age.die, m_total)~logpat.q91+factor(income)+factor(raceth)+factor(obesity)
                          +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.male)
summary(total.logpat.q91.m)

## female
total.logbpb.q91.f<-svycoxph(Surv(hsaitmor, age.die, m_total)~logbpb.q91+factor(income)+factor(raceth)+factor(obesity)
                            +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.female)
summary(total.logbpb.q91.f)

total.logtib.q91.f<-svycoxph(Surv(hsaitmor, age.die, m_total)~logtib.q91+factor(income)+factor(raceth)+factor(obesity)
                            +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.female)
summary(total.logtib.q91.f)

total.logpat.q91.f<-svycoxph(Surv(hsaitmor, age.die, m_total)~logpat.q91+factor(income)+factor(raceth)+factor(obesity)
                            +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.female)
summary(total.logpat.q91.f)

### interaction
total.logbpb.q91<-svycoxph(Surv(hsaitmor, age.die, m_total)~factor(sex)*(logbpb.q91+factor(income)+factor(raceth)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c), bpdsn)
summary(total.logbpb.q91)

total.logtib.q91<-svycoxph(Surv(hsaitmor, age.die, m_total)~factor(sex)*(logtib.q91+factor(income)+factor(raceth)+factor(obesity)
                                                                          +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c), bpdsn)
summary(total.logtib.q91)

total.logpat.q91<-svycoxph(Surv(hsaitmor, age.die, m_total)~factor(sex)*(logpat.q91+factor(income)+factor(raceth)+factor(obesity)
                                                                          +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c), bpdsn)
summary(total.logpat.q91)


########################### CVD MORTALITY ################################
## male
cvd.logbpb.q91.m<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logbpb.q91+factor(income)+factor(raceth)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.male)
summary(cvd.logbpb.q91.m)

cvd.logtib.q91.m<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logtib.q91+factor(income)+factor(raceth)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.male)
summary(cvd.logtib.q91.m)

cvd.logpat.q91.m<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logpat.q91+factor(income)+factor(raceth)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.male)
summary(cvd.logpat.q91.m)

## female
cvd.logbpb.q91.f<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logbpb.q91+factor(income)+factor(raceth)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.female)
summary(cvd.logbpb.q91.f)

cvd.logtib.q91.f<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logtib.q91+factor(income)+factor(raceth)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.female)
summary(cvd.logtib.q91.f)

cvd.logpat.q91.f<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logpat.q91+factor(income)+factor(raceth)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.female)
summary(cvd.logpat.q91.f)

### interaction
cvd.logbpb.q91<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~factor(sex)*(logbpb.q91+factor(income)+factor(raceth)+factor(obesity)
                                                                          +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c), bpdsn)
summary(cvd.logbpb.q91)

cvd.logtib.q91<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~factor(sex)*(logtib.q91+factor(income)+factor(raceth)+factor(obesity)
                                                                          +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c), bpdsn)
summary(cvd.logtib.q91)

cvd.logpat.q91<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~factor(sex)*(logpat.q91+factor(income)+factor(raceth)+factor(obesity)
                                                                          +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c), bpdsn)
summary(cvd.logpat.q91)


########################### HEART MORTALITY ################################
## male
heart.logbpb.q91.m<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logbpb.q91+factor(income)+factor(raceth)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.male)
summary(heart.logbpb.q91.m)

heart.logtib.q91.m<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logtib.q91+factor(income)+factor(raceth)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.male)
summary(heart.logtib.q91.m)

heart.logpat.q91.m<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logpat.q91+factor(income)+factor(raceth)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.male)
summary(heart.logpat.q91.m)

## female
heart.logbpb.q91.f<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logbpb.q91+factor(income)+factor(raceth)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.female)
summary(heart.logbpb.q91.f)

heart.logtib.q91.f<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logtib.q91+factor(income)+factor(raceth)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.female)
summary(heart.logtib.q91.f)

heart.logpat.q91.f<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logpat.q91+factor(income)+factor(raceth)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.female)
summary(heart.logpat.q91.f)

### interaction
heart.logbpb.q91<-svycoxph(Surv(hsaitmor, age.die, m_heart)~factor(sex)*(logbpb.q91+factor(income)+factor(raceth)+factor(obesity)
                                                                        +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c), bpdsn)
summary(heart.logbpb.q91)

heart.logtib.q91<-svycoxph(Surv(hsaitmor, age.die, m_heart)~factor(sex)*(logtib.q91+factor(income)+factor(raceth)+factor(obesity)
                                                                        +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c), bpdsn)
summary(heart.logtib.q91)

heart.logpat.q91<-svycoxph(Surv(hsaitmor, age.die, m_heart)~factor(sex)*(logpat.q91+factor(income)+factor(raceth)+factor(obesity)
                                                                        +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c), bpdsn)
summary(heart.logpat.q91)


############################
########## race ############
############################
table(raceth)
svytable(~raceth, bpdsn)

pbmort3.wh<-subset(pbmort3, raceth==1)
bpdsn.wh<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=pbmort3.wh, nest=T)

pbmort3.bl<-subset(pbmort3, raceth==2)
bpdsn.bl<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=pbmort3.bl, nest=T)

pbmort3.mx<-subset(pbmort3, raceth==3)
bpdsn.mx<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=pbmort3.mx, nest=T)


######### total #########
## White
total.logbpb.q91.wh<-svycoxph(Surv(hsaitmor, age.die, m_total)~logbpb.q91+factor(income)+factor(sex)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.wh)
summary(total.logbpb.q91.wh)

total.logtib.q91.wh<-svycoxph(Surv(hsaitmor, age.die, m_total)~logtib.q91+factor(income)+factor(sex)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.wh)
summary(total.logtib.q91.wh)

total.logpat.q91.wh<-svycoxph(Surv(hsaitmor, age.die, m_total)~logpat.q91+factor(income)+factor(sex)+factor(obesity)
                              +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.wh)
summary(total.logpat.q91.wh)

## Black
total.logbpb.q91.bl<-svycoxph(Surv(hsaitmor, age.die, m_total)~logbpb.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.bl)
summary(total.logbpb.q91.bl)

total.logtib.q91.bl<-svycoxph(Surv(hsaitmor, age.die, m_total)~logtib.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.bl)
summary(total.logtib.q91.bl)

total.logpat.q91.bl<-svycoxph(Surv(hsaitmor, age.die, m_total)~logpat.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.bl)
summary(total.logpat.q91.bl)

## Mexican
total.logbpb.q91.mx<-svycoxph(Surv(hsaitmor, age.die, m_total)~logbpb.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.mx)
summary(total.logbpb.q91.mx)

total.logtib.q91.mx<-svycoxph(Surv(hsaitmor, age.die, m_total)~logtib.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.mx)
summary(total.logtib.q91.mx)

total.logpat.q91.mx<-svycoxph(Surv(hsaitmor, age.die, m_total)~logpat.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.mx)
summary(total.logpat.q91.mx)

### interaction
#total.logbpb.q91<-svycoxph(Surv(hsaitmor, age.die, m_total)~factor(raceth)*(logbpb.q91+factor(income)+factor(sex)+factor(obesity)
#                                                                          +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c), bpdsn)
#summary(total.logbpb.q91)
#Error in solve.default(g$var, coef(g)) : 
#system is computationally singular: reciprocal condition number = 6.27254e-22

total.logbpb.q91<-svycoxph(Surv(hsaitmor, age.die, m_total)~factor(raceth)*logbpb.q91+factor(raceth)*factor(income)+factor(raceth)*factor(sex)+factor(raceth)*factor(obesity)
                   +factor(raceth)*factor(smk)+factor(raceth)*factor(htn_bp)+factor(raceth)*factor(ucd_cr3)+factor(raceth)*factor(alc4)+factor(raceth)*factor(phyact1)+factor(raceth)*factor(hei3)+chol+a1c, bpdsn)
summary(total.logbpb.q91)

total.logtib.q91<-svycoxph(Surv(hsaitmor, age.die, m_total)~factor(raceth)*logtib.q91+factor(raceth)*factor(income)+factor(raceth)*factor(sex)+factor(raceth)*factor(obesity)
                            +factor(raceth)*factor(smk)+factor(raceth)*factor(htn_bp)+factor(raceth)*factor(ucd_cr3)+factor(raceth)*factor(alc4)+factor(raceth)*factor(phyact1)+factor(raceth)*factor(hei3)+chol+a1c, bpdsn)
summary(total.logtib.q91)

total.logpat.q91<-svycoxph(Surv(hsaitmor, age.die, m_total)~factor(raceth)*logpat.q91+factor(raceth)*factor(income)+factor(raceth)*factor(sex)+factor(raceth)*factor(obesity)
                            +factor(raceth)*factor(smk)+factor(raceth)*factor(htn_bp)+factor(raceth)*factor(ucd_cr3)+factor(raceth)*factor(alc4)+factor(raceth)*factor(phyact1)+factor(raceth)*factor(hei3)+chol+a1c, bpdsn)
summary(total.logpat.q91)

total.logpat.q91<-svycoxph(Surv(hsaitmor, age.die, m_total)~factor(raceth)*(logpat.q91+factor(income)+factor(sex)+factor(obesity)
                                                                          +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c), bpdsn)
summary(total.logpat.q91)

######### CVD #########
## White
cvd.logbpb.q91.wh<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logbpb.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.wh)
summary(cvd.logbpb.q91.wh)

cvd.logtib.q91.wh<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logtib.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.wh)
summary(cvd.logtib.q91.wh)

cvd.logpat.q91.wh<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logpat.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.wh)
summary(cvd.logpat.q91.wh)

## Black
cvd.logbpb.q91.bl<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logbpb.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.bl)
summary(cvd.logbpb.q91.bl)

cvd.logtib.q91.bl<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logtib.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.bl)
summary(cvd.logtib.q91.bl)

cvd.logpat.q91.bl<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logpat.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.bl)
summary(cvd.logpat.q91.bl)

## Mexican
cvd.logbpb.q91.mx<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logbpb.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.mx)
summary(cvd.logbpb.q91.mx)

cvd.logtib.q91.mx<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logtib.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.mx)
summary(cvd.logtib.q91.mx)

cvd.logpat.q91.mx<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logpat.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.mx)
summary(cvd.logpat.q91.mx)

### interaction
cvd.logbpb.q91<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~factor(raceth)*logbpb.q91+factor(raceth)*factor(income)+factor(raceth)*factor(sex)+factor(raceth)*factor(obesity)
                            +factor(raceth)*factor(smk)+factor(raceth)*factor(htn_bp)+factor(raceth)*factor(ucd_cr3)+factor(raceth)*factor(alc4)+factor(raceth)*factor(phyact1)+factor(raceth)*factor(hei3)+chol+a1c, bpdsn)
summary(cvd.logbpb.q91)

cvd.logtib.q91<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~factor(raceth)*logtib.q91+factor(raceth)*factor(income)+factor(raceth)*factor(sex)+factor(raceth)*factor(obesity)
                            +factor(raceth)*factor(smk)+factor(raceth)*factor(htn_bp)+factor(raceth)*factor(ucd_cr3)+factor(raceth)*factor(alc4)+factor(raceth)*factor(phyact1)+factor(raceth)*factor(hei3)+chol+a1c, bpdsn)
summary(cvd.logtib.q91)

cvd.logpat.q91<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~factor(raceth)*logpat.q91+factor(raceth)*factor(income)+factor(raceth)*factor(sex)+factor(raceth)*factor(obesity)
                            +factor(raceth)*factor(smk)+factor(raceth)*factor(htn_bp)+factor(raceth)*factor(ucd_cr3)+factor(raceth)*factor(alc4)+factor(raceth)*factor(phyact1)+factor(raceth)*factor(hei3)+chol+a1c, bpdsn)
summary(cvd.logpat.q91)

#cvd.logpat.q91<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~factor(raceth)*(logpat.q91+factor(income)+factor(sex)+factor(obesity)
#                                                                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c), bpdsn)
#summary(cvd.logpat.q91)

######### Heart #########
## White
heart.logbpb.q91.wh<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logbpb.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.wh)
summary(heart.logbpb.q91.wh)

heart.logtib.q91.wh<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logtib.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.wh)
summary(heart.logtib.q91.wh)

heart.logpat.q91.wh<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logpat.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.wh)
summary(heart.logpat.q91.wh)

## Black
heart.logbpb.q91.bl<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logbpb.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.bl)
summary(heart.logbpb.q91.bl)

heart.logtib.q91.bl<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logtib.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.bl)
summary(heart.logtib.q91.bl)

heart.logpat.q91.bl<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logpat.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.bl)
summary(heart.logpat.q91.bl)

## Mexican
heart.logbpb.q91.mx<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logbpb.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.mx)
summary(heart.logbpb.q91.mx)

heart.logtib.q91.mx<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logtib.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.mx)
summary(heart.logtib.q91.mx)

heart.logpat.q91.mx<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logpat.q91+factor(income)+factor(sex)+factor(obesity)
                               +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn.mx)
summary(heart.logpat.q91.mx)

### interaction
heart.logbpb.q91<-svycoxph(Surv(hsaitmor, age.die, m_heart)~factor(raceth)*logbpb.q91+factor(raceth)*factor(income)+factor(raceth)*factor(sex)+factor(raceth)*factor(obesity)
                            +factor(raceth)*factor(smk)+factor(raceth)*factor(htn_bp)+factor(raceth)*factor(ucd_cr3)+factor(raceth)*factor(alc4)+factor(raceth)*factor(phyact1)+factor(raceth)*factor(hei3)+chol+a1c, bpdsn)
summary(heart.logbpb.q91)

heart.logtib.q91<-svycoxph(Surv(hsaitmor, age.die, m_heart)~factor(raceth)*logtib.q91+factor(raceth)*factor(income)+factor(raceth)*factor(sex)+factor(raceth)*factor(obesity)
                            +factor(raceth)*factor(smk)+factor(raceth)*factor(htn_bp)+factor(raceth)*factor(ucd_cr3)+factor(raceth)*factor(alc4)+factor(raceth)*factor(phyact1)+factor(raceth)*factor(hei3)+chol+a1c, bpdsn)
summary(heart.logtib.q91)

heart.logpat.q91<-svycoxph(Surv(hsaitmor, age.die, m_heart)~factor(raceth)*logpat.q91+factor(raceth)*factor(income)+factor(raceth)*factor(sex)+factor(raceth)*factor(obesity)
                            +factor(raceth)*factor(smk)+factor(raceth)*factor(htn_bp)+factor(raceth)*factor(ucd_cr3)+factor(raceth)*factor(alc4)+factor(raceth)*factor(phyact1)+factor(raceth)*factor(hei3)+chol+a1c, bpdsn)
summary(heart.logpat.q91)

#heart.logpat.q91<-svycoxph(Surv(hsaitmor, age.die, m_heart)~factor(raceth)*(logpat.q91+factor(income)+factor(sex)+factor(obesity)
#                                                                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c), bpdsn)
#summary(heart.logpat.q91)

### forest plots
library(ggplot2)

modif<-c("Blood Lead","Blood Lead","Blood Lead","Blood Lead","Blood Lead","Blood Lead","Blood Lead","Blood Lead",
         "Tibia Lead","Tibia Lead","Tibia Lead","Tibia Lead","Tibia Lead","Tibia Lead","Tibia Lead","Tibia Lead",
         "Patella Lead","Patella Lead","Patella Lead","Patella Lead","Patella Lead","Patella Lead","Patella Lead","Patella Lead")
group <- c("All","","Male","Female","","NHW","NHB","MA",
           "All","","Male","Female","","NHW","NHB","MA",
           "All","","Male","Female","","NHW","NHB","MA")
modif1<-factor(modif, levels=unique(modif))
group1<-factor(group, levels=unique(group))

### Total
HR.total  <- c(summary(svycox.logbpb.q91)$conf.int[1,1],NA,summary(total.logbpb.q91.m)$conf.int[1,1],summary(total.logbpb.q91.f)$conf.int[1,1],NA,summary(total.logbpb.q91.wh)$conf.int[1,1],summary(total.logbpb.q91.bl)$conf.int[1,1],summary(total.logbpb.q91.mx)$conf.int[1,1],
               summary(svycox.logtib.q91)$conf.int[1,1],NA,summary(total.logtib.q91.m)$conf.int[1,1],summary(total.logtib.q91.f)$conf.int[1,1],NA,summary(total.logtib.q91.wh)$conf.int[1,1],summary(total.logtib.q91.bl)$conf.int[1,1],summary(total.logtib.q91.mx)$conf.int[1,1],
               summary(svycox.logpat.q91)$conf.int[1,1],NA,summary(total.logpat.q91.m)$conf.int[1,1],summary(total.logpat.q91.f)$conf.int[1,1],NA,summary(total.logpat.q91.wh)$conf.int[1,1],summary(total.logpat.q91.bl)$conf.int[1,1],summary(total.logpat.q91.mx)$conf.int[1,1]) 
lower.total <- c(summary(svycox.logbpb.q91)$conf.int[1,3],NA,summary(total.logbpb.q91.m)$conf.int[1,3],summary(total.logbpb.q91.f)$conf.int[1,3],NA,summary(total.logbpb.q91.wh)$conf.int[1,3],summary(total.logbpb.q91.bl)$conf.int[1,3],summary(total.logbpb.q91.mx)$conf.int[1,3],
                 summary(svycox.logtib.q91)$conf.int[1,3],NA,summary(total.logtib.q91.m)$conf.int[1,3],summary(total.logtib.q91.f)$conf.int[1,3],NA,summary(total.logtib.q91.wh)$conf.int[1,3],summary(total.logtib.q91.bl)$conf.int[1,3],summary(total.logtib.q91.mx)$conf.int[1,3],
                 summary(svycox.logpat.q91)$conf.int[1,3],NA,summary(total.logpat.q91.m)$conf.int[1,3],summary(total.logpat.q91.f)$conf.int[1,3],NA,summary(total.logpat.q91.wh)$conf.int[1,3],summary(total.logpat.q91.bl)$conf.int[1,3],summary(total.logpat.q91.mx)$conf.int[1,3]) 
upper.total <- c(summary(svycox.logbpb.q91)$conf.int[1,4],NA,summary(total.logbpb.q91.m)$conf.int[1,4],summary(total.logbpb.q91.f)$conf.int[1,4],NA,summary(total.logbpb.q91.wh)$conf.int[1,4],summary(total.logbpb.q91.bl)$conf.int[1,4],summary(total.logbpb.q91.mx)$conf.int[1,4],
                 summary(svycox.logtib.q91)$conf.int[1,4],NA,summary(total.logtib.q91.m)$conf.int[1,4],summary(total.logtib.q91.f)$conf.int[1,4],NA,summary(total.logtib.q91.wh)$conf.int[1,4],summary(total.logtib.q91.bl)$conf.int[1,4],summary(total.logtib.q91.mx)$conf.int[1,4],
                 summary(svycox.logpat.q91)$conf.int[1,4],NA,summary(total.logpat.q91.m)$conf.int[1,4],summary(total.logpat.q91.f)$conf.int[1,4],NA,summary(total.logpat.q91.wh)$conf.int[1,4],summary(total.logpat.q91.bl)$conf.int[1,4],summary(total.logpat.q91.mx)$conf.int[1,4]) 

RR_total <- data.frame(modif1, group1, HR.total, lower.total, upper.total)

p1.total = ggplot(data=RR_total,
                aes(x = group1,y = HR.total, ymin = lower.total, ymax = upper.total))+
  geom_pointrange(aes(col=group1))+
  geom_hline(aes(),yintercept =1, linetype=2)+
  xlab('')+ ylab("Hazard Ratio (95% CI)")+
  geom_errorbar(aes(ymin=lower.total, ymax=upper.total, col=group1),width=0.3,cex=1)+ 
  facet_wrap(~modif1,strip.position="top",nrow=1,scales = "free_x") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(angle = 45),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12,face="bold"),
        legend.position="none")+
  scale_y_log10(breaks=c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.0,3.0))+guides(fill=FALSE)+
  ggtitle("All-cause mortality")
p1.total

### add significance stars (failed)
p1.total = ggplot(data=RR_total,
                  aes(x = group1,y = HR.total, ymin = lower.total, ymax = upper.total))+
  geom_pointrange(aes(col=group1))+
  geom_hline(aes(),yintercept =1, linetype=2)+
  xlab('')+ ylab("Hazard Ratio (95% CI)")+
  geom_errorbar(aes(ymin=lower.total, ymax=upper.total, col=group1),width=0.3,cex=1)+ 
  facet_wrap(~modif1,strip.position="top",nrow=1,scales = "free_x") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(angle = 45),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12,face="bold"),
        legend.position="none")+
  scale_y_log10(breaks=c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.0,3.0))+guides(fill=FALSE)+
  ggtitle("All-cause mortality")
p1.total




### CVD
HR.cvd  <- c(summary(svycox.logbpb.q91.cvd)$conf.int[1,1],NA,summary(cvd.logbpb.q91.m)$conf.int[1,1],summary(cvd.logbpb.q91.f)$conf.int[1,1],NA,summary(cvd.logbpb.q91.wh)$conf.int[1,1],summary(cvd.logbpb.q91.bl)$conf.int[1,1],summary(cvd.logbpb.q91.mx)$conf.int[1,1],
               summary(svycox.logtib.q91.cvd)$conf.int[1,1],NA,summary(cvd.logtib.q91.m)$conf.int[1,1],summary(cvd.logtib.q91.f)$conf.int[1,1],NA,summary(cvd.logtib.q91.wh)$conf.int[1,1],summary(cvd.logtib.q91.bl)$conf.int[1,1],summary(cvd.logtib.q91.mx)$conf.int[1,1],
               summary(svycox.logpat.q91.cvd)$conf.int[1,1],NA,summary(cvd.logpat.q91.m)$conf.int[1,1],summary(cvd.logpat.q91.f)$conf.int[1,1],NA,summary(cvd.logpat.q91.wh)$conf.int[1,1],summary(cvd.logpat.q91.bl)$conf.int[1,1],summary(cvd.logpat.q91.mx)$conf.int[1,1]) 
lower.cvd <- c(summary(svycox.logbpb.q91.cvd)$conf.int[1,3],NA,summary(cvd.logbpb.q91.m)$conf.int[1,3],summary(cvd.logbpb.q91.f)$conf.int[1,3],NA,summary(cvd.logbpb.q91.wh)$conf.int[1,3],summary(cvd.logbpb.q91.bl)$conf.int[1,3],summary(cvd.logbpb.q91.mx)$conf.int[1,3],
                 summary(svycox.logtib.q91.cvd)$conf.int[1,3],NA,summary(cvd.logtib.q91.m)$conf.int[1,3],summary(cvd.logtib.q91.f)$conf.int[1,3],NA,summary(cvd.logtib.q91.wh)$conf.int[1,3],summary(cvd.logtib.q91.bl)$conf.int[1,3],summary(cvd.logtib.q91.mx)$conf.int[1,3],
                 summary(svycox.logpat.q91.cvd)$conf.int[1,3],NA,summary(cvd.logpat.q91.m)$conf.int[1,3],summary(cvd.logpat.q91.f)$conf.int[1,3],NA,summary(cvd.logpat.q91.wh)$conf.int[1,3],summary(cvd.logpat.q91.bl)$conf.int[1,3],summary(cvd.logpat.q91.mx)$conf.int[1,3]) 
upper.cvd <- c(summary(svycox.logbpb.q91.cvd)$conf.int[1,4],NA,summary(cvd.logbpb.q91.m)$conf.int[1,4],summary(cvd.logbpb.q91.f)$conf.int[1,4],NA,summary(cvd.logbpb.q91.wh)$conf.int[1,4],summary(cvd.logbpb.q91.bl)$conf.int[1,4],summary(cvd.logbpb.q91.mx)$conf.int[1,4],
                 summary(svycox.logtib.q91.cvd)$conf.int[1,4],NA,summary(cvd.logtib.q91.m)$conf.int[1,4],summary(cvd.logtib.q91.f)$conf.int[1,4],NA,summary(cvd.logtib.q91.wh)$conf.int[1,4],summary(cvd.logtib.q91.bl)$conf.int[1,4],summary(cvd.logtib.q91.mx)$conf.int[1,4],
                 summary(svycox.logpat.q91.cvd)$conf.int[1,4],NA,summary(cvd.logpat.q91.m)$conf.int[1,4],summary(cvd.logpat.q91.f)$conf.int[1,4],NA,summary(cvd.logpat.q91.wh)$conf.int[1,4],summary(cvd.logpat.q91.bl)$conf.int[1,4],summary(cvd.logpat.q91.mx)$conf.int[1,4]) 

#################################################################
#HR.cvd  <- c(1.6339,NA,1.5596,1.6918,NA,1.6385,1.4757,1.0404,
#             3.1633,NA,2.9232,3.0488,NA,3.1555,2.3163,1.0061,
#             2.3055,NA,2.2662,2.1232,NA,2.1614,2.5371,1.1285) 
#lower.cvd <- c(1.2467,NA,1.021,1.2792,NA,1.2254,1.0159,0.6275,
#               1.8666,NA,1.4921,1.6328,NA,1.7613,1.0828,0.3884,
#               1.4920,NA,1.3623,1.2374,NA,1.3261,1.513,0.5728)
#upper.cvd <- c(2.1414,NA,2.3822,2.2377,NA,2.1908,2.1436,1.725,
#               5.3608,NA,5.7266,5.6926,NA,5.6535,4.955,2.6065,
#               3.5624,NA,3.77,3.6432,NA,3.5228,4.2543,2.223)
#################################################################

RR_cvd <- data.frame(modif1, group1, HR.cvd, lower.cvd, upper.cvd)

p1.cvd = ggplot(data=RR_cvd,
            aes(x = group1,y = HR.cvd, ymin = lower.cvd, ymax = upper.cvd))+
  geom_pointrange(aes(col=group1))+
  geom_hline(aes(),yintercept =1, linetype=2)+
  xlab('')+ ylab("Hazard Ratio (95% CI)")+
  geom_errorbar(aes(ymin=lower.cvd, ymax=upper.cvd, col=group1),width=0.3,cex=1)+ 
  facet_wrap(~modif1,strip.position="top",nrow=1,scales = "free_x") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(angle = 45),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12,face="bold"),
        legend.position="none")+
  scale_y_log10(breaks=c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.0,3.0))+guides(fill=FALSE)+
  ggtitle("CVD mortality")
p1.cvd

par(mfrow=c(2,1))
p1.total
p1.cvd

grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1)))
## and then use the print function to print the objects into the viewport.
print(p1.total, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(p1.cvd, vp=viewport(layout.pos.row=2, layout.pos.col=1))


# same colors within each lead
p2 = ggplot(data=RR_data,
            aes(x = group1,y = HR, ymin = lower, ymax = upper))+
  geom_pointrange(aes(col=modif1))+
  geom_hline(aes(),yintercept =1, linetype=2)+
  xlab('')+ ylab("Hazard Ratio (95% CI)")+
  geom_errorbar(aes(ymin=lower, ymax=upper, col=modif1),width=0.2,cex=1)+ 
  facet_wrap(~modif1,strip.position="top",nrow=1,scales = "free_x") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12,face="bold"),
        legend.position="none")+
  scale_y_log10(breaks=c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.0,3.0))+guides(fill=FALSE)+
  ggtitle("CVD mortality")
p2

#####################################################################################
### DO NOT USE

###########################################################################
## Using a single interaction (i.e., interaction only b/w race and lead) ##
###########################################################################

names(pbmort3)
table(sex)
table(raceth)
pbmort3$fem<-ifelse(pbmort3$sex==2,1,0)
pbmort3$nhw<-ifelse(pbmort3$raceth==1,1,0)
pbmort3$nhb<-ifelse(pbmort3$raceth==2,1,0)
pbmort3$ma<-ifelse(pbmort3$raceth==3,1,0)

bpdsn<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=pbmort3, nest=T)

###########################
### interaction w/ race ###
###########################
total.logbpb.int<-svycoxph(Surv(hsaitmor, age.die, m_total)~factor(raceth)*logbpb.q91+factor(income)+factor(sex)+factor(obesity)
                           +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(total.logbpb.int)

total.logbpb.int2<-svycoxph(Surv(hsaitmor, age.die, m_total)~nhw*logbpb.q91+ma*logbpb.q91+factor(income)+factor(sex)+factor(obesity)
                           +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(total.logbpb.int2)

total.logbpb.int3<-svycoxph(Surv(hsaitmor, age.die, m_total)~nhw*logbpb.q91+nhb*logbpb.q91+factor(income)+factor(sex)+factor(obesity)
                            +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(total.logbpb.int3)


total.logtib.int<-svycoxph(Surv(hsaitmor, age.die, m_total)~factor(raceth)*logtib.q91+factor(income)+factor(sex)+factor(obesity)
                           +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(total.logtib.int)

total.logtib.int2<-svycoxph(Surv(hsaitmor, age.die, m_total)~nhw*logtib.q91+ma*logtib.q91+factor(income)+factor(sex)+factor(obesity)
                           +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(total.logtib.int2)

total.logtib.int3<-svycoxph(Surv(hsaitmor, age.die, m_total)~nhw*logtib.q91+nhb*logtib.q91+factor(income)+factor(sex)+factor(obesity)
                            +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(total.logtib.int3)

### DO NOT USE
#####################################################################################


##############################################################
######## Comparison between bpb and bone lead ########
##############################################################

### June 16, 2023
## check bpb distribution
quantile(bpb, c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8,0.9,0.95))
names(pbmort3)
tab1(bpb2)
median(bpb)
tapply(bpb, bpb2, summ)

tabpct(bpb2, tibx4, percent="column")
tabpct(bpb2, patx4, percent="column")

ggplot(pbmort3, aes(bpb, pred_tibx))+geom_point()+geom_smooth()
ggplot(pbmort3, aes(logbpb, logtib))+geom_point()+geom_smooth()

ggplot(pbmort3, aes(bpb, pred_patx))+geom_point()+geom_smooth()
ggplot(pbmort3, aes(logbpb, logpat))+geom_point()+geom_smooth()

bpb.5<-cut2(bpb, 5)
tibx5<-cut2(pred_tibx, g=5)
patx5<-cut2(pred_patx, g=5)

table(bpb.5)
table(tibx5)
table(patx5)

tabpct(bpb.5, tibx5)
tabpct(bpb.5, patx5)

### survey weighted

svytable(~bpb.5, Ntotal=100, design=bpdsn)

svytable(~bpb.5+tibx5, design=bpdsn)
svytable(~bpb.5+patx5, design=bpdsn)

svytable(~bpb.5+tibx5, Ntotal=100, design=bpdsn)
svytable(~bpb.5+patx5, Ntotal=100, design=bpdsn)

svyby(~pred_tibx, ~bpb.5, bpdsn, svymean)
svyby(~pred_patx, ~bpb.5, bpdsn, svymean)

svyby(~pred_tibx, ~bpb.5, bpdsn, svyquantile, quantiles=c(0,0.1,0.25,0.5,0.75,0.9,1))
svyby(~pred_patx, ~bpb.5, bpdsn, svyquantile, quantiles=c(0,0.1,0.25,0.5,0.75,0.9,1))

### age cut at 50
pbmort3.age49<-subset(pbmort3, age<50)
bpdsn.age49<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=pbmort3.age49, nest=T)

pbmort3.age50<-subset(pbmort3, age>=50)
bpdsn.age50<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=pbmort3.age50, nest=T)

svyby(~pred_tibx, ~bpb.5, bpdsn.age49, svyquantile, quantiles=c(0,0.1,0.25,0.5,0.75,0.9,1))
svyby(~pred_patx, ~bpb.5, bpdsn.age49, svyquantile, quantiles=c(0,0.1,0.25,0.5,0.75,0.9,1))

svyby(~pred_tibx, ~bpb.5, bpdsn.age50, svyquantile, quantiles=c(0,0.1,0.25,0.5,0.75,0.9,1))
svyby(~pred_patx, ~bpb.5, bpdsn.age50, svyquantile, quantiles=c(0,0.1,0.25,0.5,0.75,0.9,1))




##############################################################
######## Sensitivity Analysis: 1. time-on-study ########
##############################################################

names(pbmort3)

# Total
svycox.logbpb.pmon<-svycoxph(Surv(pmon_mec, m_total)~logbpb.q91+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                            +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logbpb.pmon)

svycox.logtib.pmon<-svycoxph(Surv(pmon_mec, m_total)~logtib.q91+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                            +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logtib.pmon)

svycox.logpat.pmon<-svycoxph(Surv(pmon_mec, m_total)~logpat.q91+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                            +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logpat.pmon)

# CVD
svycox.logbpb.pmon<-svycoxph(Surv(pmon_mec, m_cvd)~logbpb.q91+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logbpb.pmon)

svycox.logtib.pmon<-svycoxph(Surv(pmon_mec, m_cvd)~logtib.q91+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logtib.pmon)

svycox.logpat.pmon<-svycoxph(Surv(pmon_mec, m_cvd)~logpat.q91+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logpat.pmon)

# Heart
svycox.logbpb.pmon<-svycoxph(Surv(pmon_mec, m_heart)~logbpb.q91+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logbpb.pmon)

svycox.logtib.pmon<-svycoxph(Surv(pmon_mec, m_heart)~logtib.q91+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logtib.pmon)

svycox.logpat.pmon<-svycoxph(Surv(pmon_mec, m_heart)~logpat.q91+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logpat.pmon)

### tertile
# Total
svycox.logbpb.pmon<-svycoxph(Surv(pmon_mec, m_total)~bpb3+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logbpb.pmon)
svycox.logbpb.pmon<-svycoxph(Surv(pmon_mec, m_total)~as.numeric(bpb3)+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logbpb.pmon)

svycox.logtib.pmon<-svycoxph(Surv(pmon_mec, m_total)~tibm3+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logtib.pmon)
svycox.logtib.pmon<-svycoxph(Surv(pmon_mec, m_total)~as.numeric(tibm3)+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logtib.pmon)

svycox.logpat.pmon<-svycoxph(Surv(pmon_mec, m_total)~patm3+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logpat.pmon)
svycox.logpat.pmon<-svycoxph(Surv(pmon_mec, m_total)~as.numeric(patm3)+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logpat.pmon)

# CVD
svycox.logbpb.pmon<-svycoxph(Surv(pmon_mec, m_cvd)~bpb3+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logbpb.pmon)
svycox.logbpb.pmon<-svycoxph(Surv(pmon_mec, m_cvd)~as.numeric(bpb3)+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logbpb.pmon)

svycox.logtib.pmon<-svycoxph(Surv(pmon_mec, m_cvd)~tibm3+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logtib.pmon)
svycox.logtib.pmon<-svycoxph(Surv(pmon_mec, m_cvd)~as.numeric(tibm3)+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logtib.pmon)

svycox.logpat.pmon<-svycoxph(Surv(pmon_mec, m_cvd)~patm3+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logpat.pmon)
svycox.logpat.pmon<-svycoxph(Surv(pmon_mec, m_cvd)~as.numeric(patm3)+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logpat.pmon)

# Heart
svycox.logbpb.pmon<-svycoxph(Surv(pmon_mec, m_heart)~bpb3+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logbpb.pmon)
svycox.logbpb.pmon<-svycoxph(Surv(pmon_mec, m_heart)~as.numeric(bpb3)+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logbpb.pmon)

svycox.logtib.pmon<-svycoxph(Surv(pmon_mec, m_heart)~tibm3+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logtib.pmon)
svycox.logtib.pmon<-svycoxph(Surv(pmon_mec, m_heart)~as.numeric(tibm3)+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logtib.pmon)

svycox.logpat.pmon<-svycoxph(Surv(pmon_mec, m_heart)~patm3+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logpat.pmon)
svycox.logpat.pmon<-svycoxph(Surv(pmon_mec, m_heart)~as.numeric(patm3)+age+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                             +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c, bpdsn)
summary(svycox.logpat.pmon)


########## 6/20/2023 ###############################################################
####################################################################################
### Sensitivity analysis 2. Test the added value of including the pred bone lead ###
####################################################################################

###########################################################################################################################
#### DO NOT RUN

library(survMisc)

summ(bmi)
IQR(bmi)

svycox.cvd.pred<-svycoxph(Surv(pmon_mec, m_cvd)~I(age/10)+factor(educ)+factor(white_l)+I(bmi/5)+factor(smk)+I(packyrs/10), bpdsn)
summary(svycox.cvd.pred)

svycox.cvd.pred.logbpb<-svycoxph(Surv(pmon_mec, m_cvd)~logbpb.q91+I(age/10)+factor(educ)+factor(white_l)+I(bmi/5)+factor(smk)+I(packyrs/10), bpdsn)
summary(svycox.cvd.pred.logbpb)

svycox.cvd.pred.logtib<-svycoxph(Surv(pmon_mec, m_cvd)~logtib.q91+I(age/10)+factor(educ)+factor(white_l)+I(bmi/5)+factor(smk)+I(packyrs/10), bpdsn)
summary(svycox.cvd.pred.logtib)

svycox.cvd.pred.logpat<-svycoxph(Surv(pmon_mec, m_cvd)~logpat.q91+I(age/10)+factor(educ)+factor(white_l)+I(bmi/5)+factor(smk)+I(packyrs/10), bpdsn)
summary(svycox.cvd.pred.logpat)

regTermTest(svycox.cvd.pred.logbpb, ~logbpb.q91)
regTermTest(svycox.cvd.pred.logbpb, ~logbpb.q91, method="LRT")
regTermTest(svycox.cvd.pred.logbpb, ~logbpb.q91, method="WorkingWald")

regTermTest(svycox.cvd.pred.logtib, ~logtib.q91)
regTermTest(svycox.cvd.pred.logpat, ~logpat.q91)

svycox.cvd.pred.logtib1<-svycoxph(Surv(pmon_mec, m_cvd)~logtib.q91+logbpb.q91+I(age/10)+factor(educ)+factor(white_l)+I(bmi/5)+factor(smk)+I(packyrs/10), bpdsn)
summary(svycox.cvd.pred.logtib1)
regTermTest(svycox.cvd.pred.logtib1, ~logtib.q91)

svycox.cvd.pred.logpat1<-svycoxph(Surv(pmon_mec, m_cvd)~logpat.q91+logbpb.q91+I(age/10)+factor(educ)+factor(white_l)+I(bmi/5)+factor(smk)+I(packyrs/10), bpdsn)
summary(svycox.cvd.pred.logpat1)
regTermTest(svycox.cvd.pred.logpat1, ~logpat.q91)

### using coxph
cox.cvd.pred<-coxph(Surv(pmon_mec, m_cvd)~I(age/10)+factor(educ)+factor(white_l)+I(bmi/5)+factor(smk)+I(packyrs/10), data=pbmort3)
summary(cox.cvd.pred)

cox.cvd.pred.logbpb<-coxph(Surv(pmon_mec, m_cvd)~logbpb.q91+I(age/10)+factor(educ)+factor(white_l)+I(bmi/5)+factor(smk)+I(packyrs/10), data=pbmort3)
summary(cox.cvd.pred.logbpb)

anova(cox.cvd.pred, cox.cvd.pred.logbpb)

cox.cvd.pred.logtib<-coxph(Surv(pmon_mec, m_cvd)~logtib.q91+I(age/10)+factor(educ)+factor(white_l)+I(bmi/5)+factor(smk)+I(packyrs/10), data=pbmort3)
summary(cox.cvd.pred.logtib)
anova(cox.cvd.pred, cox.cvd.pred.logtib)

cox.cvd.pred.logpat<-coxph(Surv(pmon_mec, m_cvd)~logpat.q91+I(age/10)+factor(educ)+factor(white_l)+I(bmi/5)+factor(smk)+I(packyrs/10), data=pbmort3)
summary(cox.cvd.pred.logpat)
anova(cox.cvd.pred, cox.cvd.pred.logpat)

anova(cox.cvd.pred.logbpb, cox.cvd.pred.logtib, cox.cvd.pred.logpat)
AIC(cox.cvd.pred.logbpb, cox.cvd.pred.logtib, cox.cvd.pred.logpat)

cox.cvd.pred.logtib1<-coxph(Surv(pmon_mec, m_cvd)~logtib.q91+logbpb.q91+I(age/10)+factor(educ)+factor(white_l)+I(bmi/5)+factor(smk)+I(packyrs/10), data=pbmort3)
summary(cox.cvd.pred.logtib1)
anova(cox.cvd.pred.logtib, cox.cvd.pred.logtib1)

cox.cvd.pred.logpat1<-coxph(Surv(pmon_mec, m_cvd)~logpat.q91+logbpb.q91+I(age/10)+factor(educ)+factor(white_l)+I(bmi/5)+factor(smk)+I(packyrs/10), data=pbmort3)
summary(cox.cvd.pred.logpat1)
anova(cox.cvd.pred.logpat, cox.cvd.pred.logpat1)

#### DO NOT RUN
##########################################################################################################

### START FROM HERE

### using the bone lead predictors (educ, bmi, packyrs)

### Total
svycox.total.pred<-svycoxph(Surv(pmon_mec, m_total)~I(age/10)+I(bmi/IQR(bmi))+factor(educ)+factor(white_l)+factor(smk)+I(packyrs/IQR(packyrs)), bpdsn)
summary(svycox.total.pred)

svycox.total.pred.logbpb<-svycoxph(Surv(pmon_mec, m_total)~logbpb.q91+I(age/10)+I(bmi/IQR(bmi))+factor(educ)+factor(white_l)+factor(smk)+I(packyrs/IQR(packyrs)), bpdsn)
summary(svycox.total.pred.logbpb)

svycox.total.pred.logtib<-svycoxph(Surv(pmon_mec, m_total)~logtib.q91+I(age/10)+I(bmi/IQR(bmi))+factor(educ)+factor(white_l)+factor(smk)+I(packyrs/IQR(packyrs)), bpdsn)
summary(svycox.total.pred.logtib)

svycox.total.pred.logpat<-svycoxph(Surv(pmon_mec, m_total)~logpat.q91+I(age/10)+I(bmi/IQR(bmi))+factor(educ)+factor(white_l)+factor(smk)+I(packyrs/IQR(packyrs)), bpdsn)
summary(svycox.total.pred.logpat)

### CVD
svycox.cvd.pred<-svycoxph(Surv(pmon_mec, m_cvd)~I(age/10)+I(bmi/IQR(bmi))+factor(educ)+factor(white_l)+factor(smk)+I(packyrs/IQR(packyrs)), bpdsn)
summary(svycox.cvd.pred)

svycox.cvd.pred.logbpb<-svycoxph(Surv(pmon_mec, m_cvd)~logbpb.q91+I(age/10)+I(bmi/IQR(bmi))+factor(educ)+factor(white_l)+factor(smk)+I(packyrs/IQR(packyrs)), bpdsn)
summary(svycox.cvd.pred.logbpb)

svycox.cvd.pred.logtib<-svycoxph(Surv(pmon_mec, m_cvd)~logtib.q91+I(age/10)+I(bmi/IQR(bmi))+factor(educ)+factor(white_l)+factor(smk)+I(packyrs/IQR(packyrs)), bpdsn)
summary(svycox.cvd.pred.logtib)

svycox.cvd.pred.logpat<-svycoxph(Surv(pmon_mec, m_cvd)~logpat.q91+I(age/10)+I(bmi/IQR(bmi))+factor(educ)+factor(white_l)+factor(smk)+I(packyrs/IQR(packyrs)), bpdsn)
summary(svycox.cvd.pred.logpat)

### using the same covariates + predictors for bone lead (educ, bmi, packyrs)

### Total
svycox.logbpb.q91.total1<-svycoxph(Surv(hsaitmor, age.die, m_total)~logbpb.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                 +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+I(chol/IQR(chol))+I(a1c/IQR(a1c))+factor(educ)+factor(white_l)+I(bmi/IQR(bmi))+I(packyrs/IQR(packyrs)), bpdsn)
summary(svycox.logbpb.q91.total1)

svycox.logtib.q91.total1<-svycoxph(Surv(hsaitmor, age.die, m_total)~logtib.q91+logbpb.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                 +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+I(chol/IQR(chol))+I(a1c/IQR(a1c))+factor(educ)+factor(white_l)+I(bmi/IQR(bmi))+I(packyrs/IQR(packyrs)), bpdsn)
summary(svycox.logtib.q91.total1)
regTermTest(svycox.logtib.q91.total1, ~logtib.q91)
regTermTest(svycox.logtib.q91.total1, ~logbpb.q91)

svycox.logpat.q91.total1<-svycoxph(Surv(hsaitmor, age.die, m_total)~logpat.q91+logbpb.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                 +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+I(chol/IQR(chol))+I(a1c/IQR(a1c))+factor(educ)+factor(white_l)+I(bmi/IQR(bmi))+I(packyrs/IQR(packyrs)), bpdsn)
summary(svycox.logpat.q91.total1)
regTermTest(svycox.logpat.q91.total1, ~logpat.q91)
regTermTest(svycox.logpat.q91.total1, ~logbpb.q91)

### CVD
svycox.logbpb.q91.cvd1<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logbpb.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                 +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+I(chol/IQR(chol))+I(a1c/IQR(a1c))+factor(educ)+factor(white_l)+I(bmi/IQR(bmi))+I(packyrs/IQR(packyrs)), bpdsn)
summary(svycox.logbpb.q91.cvd1)

svycox.logtib.q91.cvd1<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logtib.q91+logbpb.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                 +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+I(chol/IQR(chol))+I(a1c/IQR(a1c))+factor(educ)+factor(white_l)+I(bmi/IQR(bmi))+I(packyrs/IQR(packyrs)), bpdsn)
summary(svycox.logtib.q91.cvd1)
regTermTest(svycox.logtib.q91.cvd1, ~logtib.q91, method="LRT")
regTermTest(svycox.logtib.q91.cvd1, ~logbpb.q91)

svycox.logpat.q91.cvd1<-svycoxph(Surv(hsaitmor, age.die, m_cvd)~logpat.q91+logbpb.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                 +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+I(chol/IQR(chol))+I(a1c/IQR(a1c))+factor(educ)+factor(white_l)+I(bmi/IQR(bmi))+I(packyrs/IQR(packyrs)), bpdsn)
summary(svycox.logpat.q91.cvd1)
regTermTest(svycox.logpat.q91.cvd1, ~logpat.q91, method="LRT")
regTermTest(svycox.logpat.q91.cvd1, ~logbpb.q91)

### heart
svycox.logbpb.q91.heart1<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logbpb.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                 +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+I(chol/IQR(chol))+I(a1c/IQR(a1c))+factor(educ)+factor(white_l)+I(bmi/IQR(bmi))+I(packyrs/IQR(packyrs)), bpdsn)
summary(svycox.logbpb.q91.heart1)

svycox.logtib.q91.heart1<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logtib.q91+logbpb.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                 +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+I(chol/IQR(chol))+I(a1c/IQR(a1c))+factor(educ)+factor(white_l)+I(bmi/IQR(bmi))+I(packyrs/IQR(packyrs)), bpdsn)
summary(svycox.logtib.q91.heart1)
regTermTest(svycox.logtib.q91.heart1, ~logtib.q91)
regTermTest(svycox.logtib.q91.heart1, ~logbpb.q91)

svycox.logpat.q91.heart1<-svycoxph(Surv(hsaitmor, age.die, m_heart)~logpat.q91+logbpb.q91+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                 +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+I(chol/IQR(chol))+I(a1c/IQR(a1c))+factor(educ)+factor(white_l)+I(bmi/IQR(bmi))+I(packyrs/IQR(packyrs)), bpdsn)
summary(svycox.logpat.q91.heart1)
regTermTest(svycox.logpat.q91.heart1, ~logpat.q91)
regTermTest(svycox.logpat.q91.heart1, ~logbpb.q91)



## using time-on-study
svycox.logtib.q91.cvd2<-svycoxph(Surv(pmon_mec, m_cvd)~logtib.q91+logbpb.q91+I(age/10)+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                 +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c+factor(educ)+factor(white_l)+I(bmi/5)+I(packyrs/10), bpdsn)
summary(svycox.logtib.q91.cvd2)
svycox.logpat.q91.cvd2<-svycoxph(Surv(pmon_mec, m_cvd)~logpat.q91+logbpb.q91+I(age/10)+factor(sex)+factor(income)+factor(raceth)+factor(obesity)
                                 +factor(smk)+factor(htn_bp)+factor(ucd_cr3)+factor(alc4)+factor(phyact1)+factor(hei3)+chol+a1c+factor(educ)+factor(white_l)+I(bmi/5)+I(packyrs/10), bpdsn)
summary(svycox.logpat.q91.cvd2)
