## Preliminary setup.
# Change variable mypath to where you unzipped the OSF folder
# Change this to where you unzipped the data package
mypath <- "C:/Users/chjh/SURFdrive/Cyberball Meta-analysis"
setwd("~/Desktop/R Github/Multiverse")

dat<-read.csv("data_cleaned.csv",sep=";",dec=".")
require(metafor);require(xtable);require(lattice)
## Function for writing the results of each analysis into a data frame for your pleasure
name<-NULL;k<-NULL;QE<-NULL;QEp<-NULL;I2<-NULL;QM<-NULL;QMp<-NULL;tau2<-NULL;tau2low<-NULL;tau2up<-NULL;ES<-NULL;ESse<-NULL;ESp<-NULL;ESLOW<-NULL;ESUP<-NULL;tabres<-NULL
write_res_tab<-function(name,object_results){
  Name<-name
  k<-object_results$k
  QE<-object_results$QE
  QE_p<-object_results$QEp
  I_2<-object_results$I2
  QM<-object_results$QM
  QM_p<-object_results$QMp
  Tau2<-object_results$tau2
  Tau2_CI_95_low<-confint(object_results)$random[1,2]
  Tau2_CI_95_up<-confint(object_results)$random[1,3]
  ES<-object_results$b[1]
  ES_SE<-object_results$se[1]
  ES_p<-object_results$pval[1]
  CI_95_low<-object_results$ci.lb[1]
  CI_95_up<-object_results$ci.ub[1]
  tabres<-round(data.frame(k,QE,QE_p,QM,QM_p,Tau2,Tau2_CI_95_low,Tau2_CI_95_up,ES,ES_SE,ES_p,CI_95_low,CI_95_up),3)
  tabres<-cbind(name,tabres)
  tabres
}

##############################################################################
# Cyberball meta-analysis source code to
# "Temporal Effects of Ostracism: A Meta-Analysis of Cyberball Studies." 
# R code primarily written by Chris HJ Hartgerink

#   Contact info:
#   c.hartgerink@gmail.com
##############################################################################
# Note: analyses are described in order of the paper. Supplemental analyses or
#      code used to format pictures, tables, etc. can be found at the end.

##########
# Method #
##########
# Descriptives sample
length(unique(temp$Paper))
length(unique(dat$Paper_nr))
length(dat$Unique_study_nr)
round(mean(dat$N),1);median(dat$N);round(sum(dat$N),0)

###########
# Results #
###########
# Confirmatory hypotheses
# Hypothesis 1 (main effect ostracism)
conf_ostr_main_T1_H1 <- rma((dat$ostr_nomod_simple_main_T1),(dat$ostr_nomod_simple_main_svar_T1),method="REML",digits=4); 
tab_conf_ostr_main_T1_H1 <- write_res_tab("Main effect ostracism T1",conf_ostr_main_T1_H1)
regtest.rma(conf_ostr_main_T1_H1,model="rma",predictor="sei")
conf_ostr_main_T2_H1 <- rma((dat$ostr_nomod_simple_main_T2),(dat$ostr_nomod_simple_main_svar_T2),mod=~as.numeric(dat$Time_est),method="REML",digits=4)
tab_conf_ostr_main_T2_H1 <- write_res_tab("Main effect ostracism T2",conf_ostr_main_T2_H1)
regtest.rma(conf_ostr_main_T2_H1,model="rma",predictor="sei")
mean(dat$Time_est)
predict(conf_ostr_main_T2_H1,mean(dat$Time_est))
rma((dat$ostr_nomod_simple_main_T2),(dat$ostr_nomod_simple_main_svar_T2),method="REML",digits=4)

# Hypothesis 2 (interaction effect)
conf_ostrmod_int_T1_H2 <- rma((dat$int_T1),(dat$int_svar_T1),method="REML",digits=4)
tab_conf_ostrmod_int_T1_H2 <- write_res_tab("Interaction effect ostracism T1",conf_ostrmod_int_T1_H2)
regtest.rma(conf_ostrmod_int_T1_H2,model="rma",predictor="sei")
conf_ostrmod_int_T2_H2 <- rma((dat$int_T2),(dat$int_svar_T2),mod=~as.numeric(dat$Time_est),method="REML",digits=4)
tab_conf_ostrmod_int_T2_H2 <- write_res_tab("Interaction effect ostracism T2",conf_ostrmod_int_T2_H2)
regtest.rma(conf_ostrmod_int_T2_H2,model="rma",predictor="sei")
rma((dat$int_T2),(dat$int_svar_T2),method="REML",digits=4)

tab_conf <- rbind(tab_conf_ostr_main_T1_H1[1,],
                  tab_conf_ostr_main_T2_H1[1,],
                  tab_conf_ostrmod_int_T1_H2[1,],
                  tab_conf_ostrmod_int_T2_H2[1,])

# Sensitivity analyses interaction on last measure (only one mentioned in the paper)
# conf_ostrmod_int_T2_H2
# Regression test (Egger's test)
funnel(conf_ostrmod_int_T2_H2,main="conf_ostrmod_int_T2_H2",addtau2=T);
# Sensitivity tau-2 estimate
tau2_conf_ostrmod_int_T2_H2 <- rma((dat$int_T2),(dat$int_svar_T2),tau2=confint(conf_ostrmod_int_T2_H2)$random[1,3]
                                   ,mod=~as.numeric(dat$Time_est)
                                   ,digits=4)
# Influencing outliers taken out
outl_conf_ostrmod_int_T2_H2 <- rma((dat$int_T2[!influence(conf_ostrmod_int_T2_H2)$inf[,1]>=2]),(dat$int_svar_T2[!influence(conf_ostrmod_int_T2_H2)$inf[,1]>=2])
                                   ,mod=~as.numeric(dat$Time_est)[!influence(conf_ostrmod_int_T2_H2)$inf[,1]>=2]
                                   ,tau2=confint(conf_ostrmod_int_T2_H2)$random[1,3],digits=4)
qqnorm(conf_ostrmod_int_T2_H2,type="rstudent")
sens_ana_conf_ostrmod_int_T2_H2 <- rbind(tab_conf_ostrmod_int_T2_H2[1,],
                                         write_res_tab("tau2_conf_ostrmod_int_T2_H2",tau2_conf_ostrmod_int_T2_H2)[1,],
                                         write_res_tab("outl_conf_ostrmod_int_T2_H2",outl_conf_ostrmod_int_T2_H2)[1,])
sens_ana_conf_ostrmod_int_T2_H2
## See the effect does become significant after removing outlier, with full Tau2.
## No asymmetry in funnel plot

#####################################################################
### Not shown in paper because these yield highly similar results ###
# Statistical Sensitivity analyses
# conf_ostr_main_T1_H1
# No different results 
par(mfrow=c(2,1))
# conf_ostr_main_T1_H1
# Regression test (Egger's test)
funnel(conf_ostr_main_T1_H1,main="Conf_ostr_main_T1_H1",addtau2=T);
# Sensitivity tau-2 estimate
tau2_conf_ostr_main_T1_H1 <- rma((dat$ostr_nomod_simple_main_T1),(dat$ostr_nomod_simple_main_svar_T1),tau2=confint(conf_ostr_main_T1_H1)$random[1,3],digits=4)
# Influencing outliers taken out
outl_conf_ostr_main_T1_H1 <- rma((dat$ostr_nomod_simple_main_T1[!influence(conf_ostr_main_T1_H1)$inf[,1]>=2]),(dat$ostr_nomod_simple_main_svar_T1[!influence(conf_ostr_main_T1_H1)$inf[,1]>=2]),tau2=confint(conf_ostr_main_T1_H1)$random[1,3],digits=4)
qqnorm(conf_ostr_main_T1_H1,type="rstudent")
sens_ana_conf_ostr_main_T1_H1 <- rbind(tab_conf_ostr_main_T1_H1[1,],
                                       write_res_tab("tau2_conf_ostr_main_T1_H1",tau2_conf_ostr_main_T1_H1)[1,],
                                       write_res_tab("outl_conf_ostr_main_T1_H1",outl_conf_ostr_main_T1_H1)[1,])
sens_ana_conf_ostr_main_T1_H1
## No different conclusions based on sensitivity analyses
## Asymmetry in funnel plot

# conf_ostr_main_T2_H1
# Regression test (Egger's test)
funnel(conf_ostr_main_T2_H1,main="conf_ostr_main_T2_H1",addtau2=T);
# Sensitivity tau-2 estimate
tau2_conf_ostr_main_T2_H1 <- rma((dat$ostr_nomod_simple_main_T2),(dat$ostr_nomod_simple_main_svar_T2),tau2=confint(conf_ostr_main_T2_H1)$random[1,3],mod=~as.numeric(dat$Time_est),digits=4)
# Influencing outliers taken out
outl_conf_ostr_main_T2_H1 <- rma((dat$ostr_nomod_simple_main_T2[!influence(conf_ostr_main_T2_H1)$inf[,1]>=2]),(dat$ostr_nomod_simple_main_svar_T2[!influence(conf_ostr_main_T2_H1)$inf[,1]>=2])
                                 ,mod=~as.numeric(dat$Time_est)[!influence(conf_ostr_main_T2_H1)$inf[,1]>=2]
                                 ,tau2=confint(conf_ostr_main_T2_H1)$random[1,3],digits=4)
qqnorm(conf_ostr_main_T2_H1,type="rstudent")
sens_ana_conf_ostr_main_T2_H1 <- rbind(tab_conf_ostr_main_T2_H1[1,],
                                       write_res_tab("tau2_conf_ostr_main_T2_H1",tau2_conf_ostr_main_T2_H1)[1,],
                                       write_res_tab("outl_conf_ostr_main_T2_H1",outl_conf_ostr_main_T2_H1)[1,])
sens_ana_conf_ostr_main_T2_H1
## No different conclusions based on sensitivity analyses
## No asymmetry in funnel plot

# conf_ostrmod_int_T1_H2
# Regression test (Egger's test)
funnel(conf_ostrmod_int_T1_H2,main="conf_ostrmod_int_T1_H2",addtau2=T);
# Sensitivity tau-2 estimate
tau2_conf_ostrmod_int_T1_H2 <- rma((dat$int_T1),(dat$int_svar_T1),tau2=confint(conf_ostrmod_int_T1_H2)$random[1,3],digits=4)
# Influencing outliers taken out
outl_conf_ostrmod_int_T1_H2 <- rma((dat$int_T1[!influence(conf_ostrmod_int_T1_H2)$inf[,1]>=2]),(dat$int_svar_T1[!influence(conf_ostrmod_int_T1_H2)$inf[,1]>=2]),tau2=confint(conf_ostrmod_int_T1_H2)$random[1,3],digits=4)
qqnorm(conf_ostrmod_int_T1_H2,type="rstudent")
sens_ana_conf_ostrmod_int_T1_H2 <- rbind(tab_conf_ostrmod_int_T1_H2[1,],
                                         write_res_tab("tau2_conf_ostrmod_int_T1_H2",tau2_conf_ostrmod_int_T1_H2)[1,],
                                         write_res_tab("outl_conf_ostrmod_int_T1_H2",outl_conf_ostrmod_int_T1_H2)[1,])
sens_ana_conf_ostrmod_int_T1_H2
## No different conclusions based on sensitivity analyses
## No asymmetry in funnel plot
### End of not shown analyses ###
#################################

# Exploratory
# overlapping fundamental needs with model
fundMod1 <- dat$Model_immediate==1 & dat$Needs1>0
fundMod2 <- dat$Model_delayed==1 & dat$Needs2>0
rma((dat$int_T1[fundMod1]),(dat$int_svar_T1[fundMod1]),method="REML",digits=4)
rma((dat$int_T2[fundMod2]),(dat$int_svar_T2[fundMod2]),method="REML",digits=4)

# Composition
# Adjust for minor error
# See also codebook
dat$NeedScale <- ifelse(dat$NeedScale == 5, 2, dat$NeedScale)
# T1
select <- !is.na(dat$Country) & !is.na(dat$PropSex) & !is.na(dat$AgeM) & !is.na(dat$PlayersGame) & !is.na(dat$ExclusionSess) & !is.na(dat$Time_est) & !is.na(dat$ostr_main_T1) &
  !influence(conf_ostr_main_T1_H1)$inf[,1]>=2
expl_design_select_ES <- dat$ostr_nomod_simple_main_T1[select]
expl_design_select_svar <- dat$ostr_nomod_simple_main_svar_T1[select]
x <- rma((expl_design_select_ES),(expl_design_select_svar),method="REML",digits=4,
    mods=~
      as.factor(dat$Country[select])
    +dat$PropSex[select]
    +dat$AgeM[select]
    +as.factor(dat$PlayersGame[select])
    +as.factor(dat$ExclusionSess[select])
    +as.factor(dat$NeedScale[select])
    +dat$ThrowLast[select]
    )
x
confint(x)
round(min(x$pval[-1]),3)

# T2
select <- !is.na(dat$Country) & !is.na(dat$PropSex) & !is.na(dat$AgeM) & !is.na(dat$PlayersGame) & !is.na(dat$ExclusionSess) & !is.na(dat$Time_est) & !is.na(dat$ostr_main_T1) &
  !influence(conf_ostr_main_T2_H1)$inf[,1]>=2
expl_design_select_ES <- dat$ostr_nomod_simple_main_T2[select]
expl_design_select_svar <- dat$ostr_nomod_simple_main_svar_T2[select]
x <- rma((expl_design_select_ES),(expl_design_select_svar),method="REML",digits=4,
    mods=~
      as.factor(dat$Country[select])
    +dat$PropSex[select]
    +dat$AgeM[select]
    +as.factor(dat$PlayersGame[select])
    +as.factor(dat$ExclusionSess[select])
    +as.factor(dat$NeedScale[select])
    +dat$ThrowLast[select]
)
x
confint(x)

# Homogeneity?
select_homogeneous<-(dat$Needs1==1|dat$Needs1==2)&dat$PlayersGame==3&dat$ExclusionSess==1&dat$ThrowLast==30 & dat$Model_immediate==1
# Homogeneous subset of Cyberball studies
rma(dat$ostr_nomod_simple_main_T1[select_homogeneous],dat$ostr_nomod_simple_main_svar_T1[select_homogeneous],method="REML",digits=4)

##############
# Discussion #
##############
# Descriptive stats from our data
# Number of players
round(sort(table(dat$PlayersGame)/120),2)
# Gender
mean(na.omit(dat$PropSex))
# Culture/country
round(sort(table(na.omit(dat$Country))/length(na.omit(dat$Country))),2)
# Age
mean(na.omit(dat$AgeM))
range(na.omit(dat$AgeM))
# Throws
round(sort(table(dat$ThrowLast)/120),2)

# Difference index
diff <-  dat$Dvlast-dat$Dvfirst
rma(dat$ostr_nomod_simple_main_T2,dat$ostr_nomod_simple_main_svar_T2,method="REML",digits=4,mods=~diff)


# substantive differences in measure proportions
x <- matrix(ncol=2,nrow=6)
x[1,1] <- sum(dat$Dvfirst==1)/length(dat$Dvfirst)
x[2,1] <- sum(dat$Dvfirst==2)/length(dat$Dvfirst)
x[3,1] <- sum(dat$Dvfirst==3)/length(dat$Dvfirst)
x[4,1] <- sum(dat$Dvfirst==4)/length(dat$Dvfirst)
x[5,1] <- sum(dat$Dvfirst==5)/length(dat$Dvfirst)
x[6,1] <- sum(dat$Dvfirst==6)/length(dat$Dvfirst)
selDvlast <- !is.na(dat$Dvlast)
x[1,2] <- sum(dat$Dvlast[selDvlast]==1)/length(dat$Dvlast[selDvlast])
x[2,2] <- sum(dat$Dvlast[selDvlast]==2)/length(dat$Dvlast[selDvlast])
x[3,2] <- sum(dat$Dvlast[selDvlast]==3)/length(dat$Dvlast[selDvlast])
x[4,2] <- sum(dat$Dvlast[selDvlast]==4)/length(dat$Dvlast[selDvlast])
x[5,2] <- sum(dat$Dvlast[selDvlast]==5)/length(dat$Dvlast[selDvlast])
x[6,2] <- sum(dat$Dvlast[selDvlast]==6)/length(dat$Dvlast[selDvlast])
max(abs(x[,1]-x[,2]))

# Calculating mean sample size for all full factorial studies
mean(na.omit(dat$N[dat$Moderator==1]))


#############
# Footnotes #
#############
# Footnote 5
# Univariate discrepancies
dat$NeedScale <- ifelse(dat$NeedScale == 5, 2, dat$NeedScale)
# T1 - number of throws
select <- !is.na(dat$Country) & !is.na(dat$PropSex) & !is.na(dat$AgeM) & !is.na(dat$PlayersGame) & !is.na(dat$ExclusionSess) & !is.na(dat$Time_est) & !is.na(dat$ostr_main_T1) &
  !influence(conf_ostr_main_T1_H1)$inf[,1]>=2
expl_design_select_ES <- dat$ostr_nomod_simple_main_T1[select]
expl_design_select_svar <- dat$ostr_nomod_simple_main_svar_T1[select]
rma((expl_design_select_ES),(expl_design_select_svar),method="REML",digits=4,
         mods=~dat$ThrowLast[select]
)
# T2 - Nr of players
select <- !is.na(dat$Country) & !is.na(dat$PropSex) & !is.na(dat$AgeM) & !is.na(dat$PlayersGame) & !is.na(dat$ExclusionSess) & !is.na(dat$Time_est) & !is.na(dat$ostr_main_T1) &
  !influence(conf_ostr_main_T2_H1)$inf[,1]>=2
expl_design_select_ES <- dat$ostr_nomod_simple_main_T2[select]
expl_design_select_svar <- dat$ostr_nomod_simple_main_svar_T2[select]
rma((expl_design_select_ES),(expl_design_select_svar),method="REML",digits=4,
         mods=~as.factor(dat$PlayersGame[select])
)

# Footnote 6
# Composite needs only
rma((dat$int_T1[dat$SelFinal1==1]),
    (dat$int_svar_T1[dat$SelFinal1==1]),method="REML",digits=4)
rma((dat$int_T2[dat$SelFinal2==1]),
    (dat$int_svar_T2[dat$SelFinal2==1]),
    mod=~as.numeric(dat$Time_est[dat$SelFinal2==1]),method="REML",digits=4)



#################
# Table 3, 4, 5 #
#################
# Table 2 is not in here because that was a direct copy paste from the masterfile.xlsx
# Table 3 #
tab3 <- cbind(dat$int_T1,dat$int_svar_T1,dat$int_T2,dat$int_svar_T2)
rma((tab3[,1]),(tab3[,2]),method="REML",digits=4)
rma((tab3[,3]),(tab3[,4]),mod=~as.numeric(dat$Time_est),method="REML",digits=4)
# Listwise deleted, not shown in table. #
rma((tab3[,1][complete.cases(tab3)]),(tab3[,2][complete.cases(tab3)]),method="REML",digits=4)

# fundamental
rma((tab3[,1])[dat$Needs1>0 & !is.na(dat$Needs1)],(tab3[,2])[dat$Needs1>0 & !is.na(dat$Needs1)],method="REML",digits=4)
rma((tab3[,3])[dat$Needs2>0 & !is.na(dat$Needs2)],(tab3[,4])[dat$Needs2>0 & !is.na(dat$Needs2)],mod=~as.numeric(dat$Time_est)[dat$Needs2>0 & !is.na(dat$Needs2)],method="REML",digits=4)
# Listwise deleted, not shown in table. #
rma((tab3[,1])[complete.cases(tab3) & dat$Needs2>0 & !is.na(dat$Needs2)],
    (tab3[,2])[complete.cases(tab3) & dat$Needs2>0 & !is.na(dat$Needs2)],method="REML",digits=4)

# Intra
rma((tab3[,1])[(dat$Dvfirst==1|dat$Dvfirst==2|dat$Dvfirst==3) & !is.na(dat$Dvfirst)],(tab3[,2])[(dat$Dvfirst==1|dat$Dvfirst==2|dat$Dvfirst==3) & !is.na(dat$Dvfirst)],method="REML",digits=4)
rma((tab3[,3])[(dat$Dvlast==1|dat$Dvlast==2|dat$Dvlast==3) & !is.na(dat$Dvlast)],(tab3[,4])[(dat$Dvlast==1|dat$Dvlast==2|dat$Dvlast==3) & !is.na(dat$Dvlast)],mod=~as.numeric(dat$Time_est)[(dat$Dvlast==1|dat$Dvlast==2|dat$Dvlast==3) & !is.na(dat$Dvlast)],method="REML",digits=4)
# Listwise deleted, not shown in table. #
rma((tab3[,1])[complete.cases(tab3) & (dat$Dvlast==1|dat$Dvlast==2|dat$Dvlast==3) & !is.na(dat$Dvlast)],
    (tab3[,2])[complete.cases(tab3) & (dat$Dvlast==1|dat$Dvlast==2|dat$Dvlast==3) & !is.na(dat$Dvlast)],method="REML",digits=4)

# Inter
rma((tab3[,1])[(dat$Dvfirst==4|dat$Dvfirst==5) & !is.na(dat$Dvfirst)],(tab3[,2])[(dat$Dvfirst==4|dat$Dvfirst==5) & !is.na(dat$Dvfirst)],method="REML",digits=4)
rma((tab3[,3])[(dat$Dvlast==4|dat$Dvlast==5) & !is.na(dat$Dvlast)],(tab3[,4])[(dat$Dvlast==4|dat$Dvlast==5) & !is.na(dat$Dvlast)],mod=~as.numeric(dat$Time_est)[(dat$Dvlast==4|dat$Dvlast==5) & !is.na(dat$Dvlast)],method="REML",digits=4)
# Listwise deleted, SHOWN in table
rma((tab3[,1])[complete.cases(tab3) & (dat$Dvlast==4|dat$Dvlast==5) & !is.na(dat$Dvlast)],
    (tab3[,2])[complete.cases(tab3) & (dat$Dvlast==4|dat$Dvlast==5) & !is.na(dat$Dvlast)],method="REML",digits=4)

# Model 
rma((tab3[,1])[dat$Model_immediate==1],(tab3[,2])[dat$Model_immediate==1],method="REML",digits=4)
rma((tab3[,3])[dat$Model_delayed==1],(tab3[,4])[dat$Model_delayed==1],mod=~as.numeric(dat$Time_est)[dat$Model_delayed==1],method="REML",digits=4)
# Listwise deleted, not shown in table. #
rma((tab3[,1])[complete.cases(tab3) & dat$Model_delayed==1],
    (tab3[,2])[complete.cases(tab3) & dat$Model_delayed==1],method="REML",digits=4)

# Table 4
select <- !is.na(dat$Country) & !is.na(dat$PropSex) & !is.na(dat$AgeM) & !is.na(dat$PlayersGame) & !is.na(dat$ExclusionSess) & !is.na(dat$Time_est) & !is.na(dat$ostr_main_T1) &
  !influence(conf_ostr_main_T1_H1)$inf[,1]>=2
expl_design_select_ES <- dat$ostr_nomod_simple_main_T1[select]
expl_design_select_svar <- dat$ostr_nomod_simple_main_svar_T1[select]
rma((expl_design_select_ES),(expl_design_select_svar),method="REML",digits=4,
         mods=~
           as.factor(dat$Country[select])
         +dat$PropSex[select]
         +dat$AgeM[select]
         +as.factor(dat$PlayersGame[select])
         +as.factor(dat$ExclusionSess[select])
         +as.factor(dat$NeedScale[select])
         +dat$ThrowLast[select])


# Table 5
select <- !is.na(dat$Country) & !is.na(dat$PropSex) & !is.na(dat$AgeM) & !is.na(dat$PlayersGame) & !is.na(dat$ExclusionSess) & !is.na(dat$Time_est) & !is.na(dat$ostr_main_T1) &
  !influence(conf_ostr_main_T2_H1)$inf[,1]>=2
expl_design_select_ES <- dat$ostr_nomod_simple_main_T2[select]
expl_design_select_svar <- dat$ostr_nomod_simple_main_svar_T2[select]
rma((expl_design_select_ES),(expl_design_select_svar),method="REML",digits=4,
         mods=~
           as.factor(dat$Country[select])
         +dat$PropSex[select]
         +dat$AgeM[select]
         +as.factor(dat$PlayersGame[select])
         +as.factor(dat$ExclusionSess[select])
         +as.factor(dat$NeedScale[select])
         +dat$ThrowLast[select])

###########
# Figures #
###########
# Figure 2
# Dotplot of all simple effects
# Running all models to get the information for the plots
# Exploring the simple effects of the interaction in an exploratory fashion
expl_ostr_mod_T1 <- rma((dat$ostr_mod_simple_main_T1),(dat$ostr_mod_simple_main_svar_T1),method="REML",digits=4);
# tab_expl_ostr_mod_T1 <- write_res_tab("Ostracism T1 (moderator level)",expl_ostr_mod_T1)
expl_ostr_nomod_T1 <- rma((dat$ostr_nomod_simple_main_T1[!is.na(dat$ostr_mod_simple_main_T1)]),(dat$ostr_nomod_simple_main_svar_T1[!is.na(dat$ostr_mod_simple_main_T1)]),method="REML",digits=4);
# tab_expl_ostr_nomod_T1 <- write_res_tab("Ostracism T1 (no-moderator level)",expl_ostr_nomod_T1)
expl_mod_ostr_T1 <- rma((dat$mod_ostr_simple_main_T1),(dat$mod_ostr_simple_main_svar_T1),method="REML",digits=4);
# tab_expl_mod_ostr_T1 <- write_res_tab("Moderater T1 (ostracism level)",expl_mod_ostr_T1)
expl_mod_incl_T1 <- rma((dat$mod_incl_simple_main_T1),(dat$mod_incl_simple_main_svar_T1),method="REML",digits=4);
# tab_expl_mod_incl_T1 <- write_res_tab("Moderator T1 (inclusion level)",expl_mod_incl_T1)

expl_ostr_mod_T2 <- rma((dat$ostr_mod_simple_main_T2),(dat$ostr_mod_simple_main_svar_T2),method="REML",digits=4);
# tab_expl_ostr_mod_T2 <- write_res_tab("Ostracism T2 (moderator level)",expl_ostr_mod_T2)
expl_ostr_nomod_T2 <- rma((dat$ostr_nomod_simple_main_T2[!is.na(dat$ostr_mod_simple_main_T2)]),(dat$ostr_nomod_simple_main_svar_T2[!is.na(dat$ostr_mod_simple_main_T2)]),method="REML",digits=4);
# tab_expl_ostr_nomod_T2 <- write_res_tab("Ostracism T2 (no-moderator level)",expl_ostr_nomod_T2)
expl_mod_ostr_T2 <- rma((dat$mod_ostr_simple_main_T2),(dat$mod_ostr_simple_main_svar_T2),method="REML",digits=4);
# tab_expl_mod_ostr_T2 <- write_res_tab("Moderater T2 (ostracism level)",expl_mod_ostr_T2)
expl_mod_incl_T2 <- rma((dat$mod_incl_simple_main_T2),(dat$mod_incl_simple_main_svar_T2),method="REML",digits=4);
# tab_expl_mod_incl_T2 <- write_res_tab("Moderator T2 (inclusion level)",expl_mod_incl_T2)
# tab_expl_simple_effects <- rbind(tab_expl_ostr_mod_T1,
#                                  tab_expl_ostr_nomod_T1,
#                                  tab_expl_mod_ostr_T1,
#                                  tab_expl_mod_incl_T1,
#                                  tab_expl_ostr_mod_T2,
#                                  tab_expl_ostr_nomod_T2,
#                                  tab_expl_mod_ostr_T2,
#                                  tab_expl_mod_incl_T2)                     
# xtable(tab_expl_simple_effects)

# Simple effects for fundamental needs
expl_ostr_mod_T1_fund <- rma((dat$ostr_mod_simple_main_T1)[dat$Needs1>0 & !is.na(dat$Needs1)],(dat$ostr_mod_simple_main_svar_T1)[dat$Needs1>0 & !is.na(dat$Needs1)],method="REML",digits=4);
# tab_expl_ostr_mod_T1_fund <- write_res_tab("Ostracism T1 (moderator level)",expl_ostr_mod_T1_fund)
expl_ostr_nomod_T1_fund <- rma((dat$ostr_nomod_simple_main_T1)[dat$Needs1>0 & !is.na(dat$Needs1) & !is.na(dat$ostr_mod_simple_main_T1)],(dat$ostr_nomod_simple_main_svar_T1)[dat$Needs1>0 & !is.na(dat$Needs1) & !is.na(dat$ostr_mod_simple_main_T1)],method="REML",digits=4);
# tab_expl_ostr_nomod_T1_fund <- write_res_tab("Ostracism T1 (no-moderator level)",expl_ostr_nomod_T1_fund)
expl_mod_ostr_T1_fund <- rma((dat$mod_ostr_simple_main_T1)[dat$Needs1>0 & !is.na(dat$Needs1)],(dat$mod_ostr_simple_main_svar_T1)[dat$Needs1>0 & !is.na(dat$Needs1)],method="REML",digits=4);
# tab_expl_mod_ostr_T1_fund <- write_res_tab("Moderator T1 (ostracism level)",expl_mod_ostr_T1_fund)
expl_mod_incl_T1_fund <- rma((dat$mod_incl_simple_main_T1)[dat$Needs1>0 & !is.na(dat$Needs1)],(dat$mod_incl_simple_main_svar_T1)[dat$Needs1>0 & !is.na(dat$Needs1)],method="REML",digits=4);
# tab_expl_mod_incl_T1_fund <- write_res_tab("Moderator T1 (inclusion level)",expl_mod_incl_T1_fund)
expl_ostr_mod_T2_fund <- rma((dat$ostr_mod_simple_main_T2)[dat$Needs2>0 & !is.na(dat$Needs2)],(dat$ostr_mod_simple_main_svar_T2)[dat$Needs2>0 & !is.na(dat$Needs2)],method="REML",digits=4);
# tab_expl_ostr_mod_T2_fund <- write_res_tab("Ostracism T2 (moderator level)",expl_ostr_mod_T2_fund)
expl_ostr_nomod_T2_fund <- rma((dat$ostr_nomod_simple_main_T2)[dat$Needs2>0 & !is.na(dat$Needs2) & !is.na(dat$ostr_mod_simple_main_T2)],(dat$ostr_nomod_simple_main_svar_T2)[dat$Needs2>0 & !is.na(dat$Needs2) & !is.na(dat$ostr_mod_simple_main_T2)],method="REML",digits=4);
# tab_expl_ostr_nomod_T2_fund <- write_res_tab("Ostracism T2 (no-moderator level)",expl_ostr_nomod_T2_fund)
expl_mod_ostr_T2_fund <- rma((dat$mod_ostr_simple_main_T2)[dat$Needs2>0 & !is.na(dat$Needs2)],(dat$mod_ostr_simple_main_svar_T2)[dat$Needs2>0 & !is.na(dat$Needs2)],method="REML",digits=4);
# tab_expl_mod_ostr_T2_fund <- write_res_tab("Moderater T2 (ostracism level)",expl_mod_ostr_T2_fund)
expl_mod_incl_T2_fund <- rma((dat$mod_incl_simple_main_T2)[dat$Needs2>0 & !is.na(dat$Needs2)],(dat$mod_incl_simple_main_svar_T2)[dat$Needs2>0 & !is.na(dat$Needs2)],method="REML",digits=4,control=list(tau2.init=1,threshold=10,maxiter=1000));
# tab_expl_mod_incl_T2_fund <- write_res_tab("Moderator T2 (inclusion level)",expl_mod_incl_T2_fund)
# rbind(tab_expl_ostr_mod_T1_fund
#       ,tab_expl_ostr_nomod_T1_fund
#       ,tab_expl_mod_ostr_T1_fund
#       ,tab_expl_mod_incl_T1_fund
#       ,tab_expl_ostr_mod_T2_fund
#       ,tab_expl_ostr_nomod_T2_fund
#       ,tab_expl_mod_ostr_T2_fund
#       ,tab_expl_mod_incl_T2_fund)

# Simple effects for intrapersonal measures
expl_ostr_mod_T1_intra <- rma((dat$ostr_mod_simple_main_T1)[(dat$Dvfirst==1|dat$Dvfirst==2|dat$Dvfirst==3) & !is.na(dat$Dvfirst)],(dat$ostr_mod_simple_main_svar_T1)[(dat$Dvfirst==1|dat$Dvfirst==2|dat$Dvfirst==3) & !is.na(dat$Dvfirst)],method="REML",digits=4);
# tab_expl_ostr_mod_T1_intra <- write_res_tab("Ostracism T1 (moderator level)",expl_ostr_mod_T1_intra)
expl_ostr_nomod_T1_intra <- rma((dat$ostr_nomod_simple_main_T1)[(dat$Dvfirst==1|dat$Dvfirst==2|dat$Dvfirst==3) & !is.na(dat$Dvfirst) & !is.na(dat$ostr_mod_simple_main_T1)],(dat$ostr_nomod_simple_main_svar_T1)[(dat$Dvfirst==1|dat$Dvfirst==2|dat$Dvfirst==3) & !is.na(dat$Dvfirst) & !is.na(dat$ostr_mod_simple_main_T1)],method="REML",digits=4);
# tab_expl_ostr_nomod_T1_intra <- write_res_tab("Ostracism T1 (no-moderator level)",expl_ostr_nomod_T1_intra)
expl_mod_ostr_T1_intra <- rma((dat$mod_ostr_simple_main_T1)[(dat$Dvfirst==1|dat$Dvfirst==2|dat$Dvfirst==3) & !is.na(dat$Dvfirst)],(dat$mod_ostr_simple_main_svar_T1)[(dat$Dvfirst==1|dat$Dvfirst==2|dat$Dvfirst==3) & !is.na(dat$Dvfirst)],method="REML",digits=4);
# tab_expl_mod_ostr_T1_intra <- write_res_tab("Moderater T1 (ostracism level)",expl_mod_ostr_T1_intra)
expl_mod_incl_T1_intra <- rma((dat$mod_incl_simple_main_T1)[(dat$Dvfirst==1|dat$Dvfirst==2|dat$Dvfirst==3) & !is.na(dat$Dvfirst)],(dat$mod_incl_simple_main_svar_T1)[(dat$Dvfirst==1|dat$Dvfirst==2|dat$Dvfirst==3) & !is.na(dat$Dvfirst)],method="REML",digits=4);
# tab_expl_mod_incl_T1_intra <- write_res_tab("Moderator T1 (inclusion level)",expl_mod_incl_T1_intra)
expl_ostr_mod_T2_intra <- rma((dat$ostr_mod_simple_main_T2)[(dat$Dvlast==1|dat$Dvlast==2|dat$Dvlast==3) & !is.na(dat$Dvlast)],(dat$ostr_mod_simple_main_svar_T2)[(dat$Dvlast==1|dat$Dvlast==2|dat$Dvlast==3) & !is.na(dat$Dvlast)],method="REML",digits=4);
# tab_expl_ostr_mod_T2_intra <- write_res_tab("Ostracism T2 (moderator level)",expl_ostr_mod_T2_intra)
expl_ostr_nomod_T2_intra <- rma((dat$ostr_nomod_simple_main_T2)[(dat$Dvlast==1|dat$Dvlast==2|dat$Dvlast==3) & !is.na(dat$Dvlast) & !is.na(dat$ostr_mod_simple_main_T2)],(dat$ostr_nomod_simple_main_svar_T2)[(dat$Dvlast==1|dat$Dvlast==2|dat$Dvlast==3) & !is.na(dat$Dvlast) & !is.na(dat$ostr_mod_simple_main_T2)],method="REML",digits=4);
# tab_expl_ostr_nomod_T2_intra <- write_res_tab("Ostracism T2 (no-moderator level)",expl_ostr_nomod_T2_intra)
expl_mod_ostr_T2_intra <- rma((dat$mod_ostr_simple_main_T2)[(dat$Dvlast==1|dat$Dvlast==2|dat$Dvlast==3) & !is.na(dat$Dvlast)],(dat$mod_ostr_simple_main_svar_T2)[(dat$Dvlast==1|dat$Dvlast==2|dat$Dvlast==3) & !is.na(dat$Dvlast)],method="REML",digits=4);
# tab_expl_mod_ostr_T2_intra <- write_res_tab("Moderater T2 (ostracism level)",expl_mod_ostr_T2_intra)
expl_mod_incl_T2_intra <- rma((dat$mod_incl_simple_main_T2)[(dat$Dvlast==1|dat$Dvlast==2|dat$Dvlast==3) & !is.na(dat$Dvlast)],(dat$mod_incl_simple_main_svar_T2)[(dat$Dvlast==1|dat$Dvlast==2|dat$Dvlast==3) & !is.na(dat$Dvlast)],method="REML",digits=4,control=list(tau2.init=1,threshold=10,maxiter=1000));
# tab_expl_mod_incl_T2_intra <- write_res_tab("Moderator T2 (inclusion level)",expl_mod_incl_T2_intra)
# rbind(tab_expl_ostr_mod_T1_intra
#       ,tab_expl_ostr_nomod_T1_intra
#       ,tab_expl_mod_ostr_T1_intra
#       ,tab_expl_mod_incl_T1_intra
#       ,tab_expl_ostr_mod_T2_intra
#       ,tab_expl_ostr_nomod_T2_intra
#       ,tab_expl_mod_ostr_T2_intra
#       ,tab_expl_mod_incl_T2_intra)

# Simple effects for interpersonal measures
expl_ostr_mod_T1_inter <- rma((dat$ostr_mod_simple_main_T1)[(dat$Dvfirst==4|dat$Dvfirst==5) & !is.na(dat$Dvfirst)],(dat$ostr_mod_simple_main_svar_T1)[(dat$Dvfirst==4|dat$Dvfirst==5) & !is.na(dat$Dvfirst)],method="REML",digits=4);
# tab_expl_ostr_mod_T1_inter <- write_res_tab("Ostracism T1 (moderator level)",expl_ostr_mod_T1_intra)
expl_ostr_nomod_T1_inter <- rma((dat$ostr_nomod_simple_main_T1)[(dat$Dvfirst==4|dat$Dvfirst==5) & !is.na(dat$Dvfirst)& !is.na(dat$ostr_mod_simple_main_T1)],(dat$ostr_nomod_simple_main_svar_T1)[(dat$Dvfirst==4|dat$Dvfirst==5) & !is.na(dat$Dvfirst)& !is.na(dat$ostr_mod_simple_main_T1)],method="REML",digits=4);
# tab_expl_ostr_nomod_T1_inter <- write_res_tab("Ostracism T1 (no-moderator level)",expl_ostr_nomod_T1_inter)
expl_mod_ostr_T1_inter <- rma((dat$mod_ostr_simple_main_T1)[(dat$Dvfirst==4|dat$Dvfirst==5) & !is.na(dat$Dvfirst)],(dat$mod_ostr_simple_main_svar_T1)[(dat$Dvfirst==4|dat$Dvfirst==5) & !is.na(dat$Dvfirst)],method="REML",digits=4);
# tab_expl_mod_ostr_T1_inter <- write_res_tab("Moderater T1 (ostracism level)",expl_mod_ostr_T1_inter)
expl_mod_incl_T1_inter <- rma((dat$mod_incl_simple_main_T1)[(dat$Dvfirst==4|dat$Dvfirst==5) & !is.na(dat$Dvfirst)],(dat$mod_incl_simple_main_svar_T1)[(dat$Dvfirst==4|dat$Dvfirst==5) & !is.na(dat$Dvfirst)],method="REML",digits=4);
# tab_expl_mod_incl_T1_inter <- write_res_tab("Moderator T1 (inclusion level)",expl_mod_incl_T1_inter)
expl_ostr_mod_T2_inter <- rma((dat$ostr_mod_simple_main_T2)[(dat$Dvlast==4|dat$Dvlast==5) & !is.na(dat$Dvlast)],(dat$ostr_mod_simple_main_svar_T2)[(dat$Dvlast==4|dat$Dvlast==5) & !is.na(dat$Dvlast)],method="REML",digits=4);
# tab_expl_ostr_mod_T2_inter <- write_res_tab("Ostracism T2 (moderator level)",expl_ostr_mod_T2_inter)
expl_ostr_nomod_T2_inter <- rma((dat$ostr_nomod_simple_main_T2)[(dat$Dvlast==4|dat$Dvlast==5) & !is.na(dat$Dvlast)& !is.na(dat$ostr_mod_simple_main_T2)],(dat$ostr_nomod_simple_main_svar_T2)[(dat$Dvlast==4|dat$Dvlast==5) & !is.na(dat$Dvlast)& !is.na(dat$ostr_mod_simple_main_T2)],method="REML",digits=4);
# tab_expl_ostr_nomod_T2_inter <- write_res_tab("Ostracism T2 (no-moderator level)",expl_ostr_nomod_T2_inter)
expl_mod_ostr_T2_inter <- rma((dat$mod_ostr_simple_main_T2)[(dat$Dvlast==4|dat$Dvlast==5) & !is.na(dat$Dvlast)],(dat$mod_ostr_simple_main_svar_T2)[(dat$Dvlast==4|dat$Dvlast==5) & !is.na(dat$Dvlast)],method="REML",digits=4);
# tab_expl_mod_ostr_T2_inter <- write_res_tab("Moderater T2 (ostracism level)",expl_mod_ostr_T2_inter)
expl_mod_incl_T2_inter <- rma((dat$mod_incl_simple_main_T2)[(dat$Dvlast==4|dat$Dvlast==5) & !is.na(dat$Dvlast)],(dat$mod_incl_simple_main_svar_T2)[(dat$Dvlast==4|dat$Dvlast==5) & !is.na(dat$Dvlast)],method="REML",digits=4,control=list(tau2.init=1,threshold=10,maxiter=1000));
# tab_expl_mod_incl_T2_inter <- write_res_tab("Moderator T2 (inclusion level)",expl_mod_incl_T2_inter)
# rbind(tab_expl_ostr_mod_T1_inter
#       ,tab_expl_ostr_nomod_T1_inter
#       ,tab_expl_mod_ostr_T1_inter
#       ,tab_expl_mod_incl_T1_inter
#       ,tab_expl_ostr_mod_T2_inter
#       ,tab_expl_ostr_nomod_T2_inter
#       ,tab_expl_mod_ostr_T2_inter
#       ,tab_expl_mod_incl_T2_inter)

# Simple effects for subset coded as being immediate as defined in the model
expl_ostr_mod_T1_model <- rma((dat$ostr_mod_simple_main_T1)[dat$Model_immediate==1],(dat$ostr_mod_simple_main_svar_T1)[dat$Model_immediate==1],method="REML",digits=4);
# tab_expl_ostr_mod_T1_model <- write_res_tab("Ostracism T1 (moderator level)",expl_ostr_mod_T1)
expl_ostr_nomod_T1_model <- rma((dat$ostr_nomod_simple_main_T1)[dat$Model_immediate==1& !is.na(dat$ostr_mod_simple_main_T1)],(dat$ostr_nomod_simple_main_svar_T1)[dat$Model_immediate==1& !is.na(dat$ostr_mod_simple_main_T1)],method="REML",digits=4);
# tab_expl_ostr_nomod_T1_model <- write_res_tab("Ostracism T1 (no-moderator level)",expl_ostr_nomod_T1)
expl_mod_ostr_T1_model <- rma((dat$mod_ostr_simple_main_T1)[dat$Model_immediate==1],(dat$mod_ostr_simple_main_svar_T1)[dat$Model_immediate==1],method="REML",digits=4);
# tab_expl_mod_ostr_T1_model <- write_res_tab("Moderater T1 (ostracism level)",expl_mod_ostr_T1)
expl_mod_incl_T1_model <- rma((dat$mod_incl_simple_main_T1)[dat$Model_immediate==1],(dat$mod_incl_simple_main_svar_T1)[dat$Model_immediate==1],method="REML",digits=4);
# tab_expl_mod_incl_T1_model <- write_res_tab("Moderator T1 (inclusion level)",expl_mod_incl_T1)
expl_ostr_mod_T2_model <- rma((dat$ostr_mod_simple_main_T2)[dat$Model_delayed==1],(dat$ostr_mod_simple_main_svar_T2)[dat$Model_delayed==1],method="REML",digits=4);
# tab_expl_ostr_mod_T2_model <- write_res_tab("Ostracism T2 (moderator level)",expl_ostr_mod_T2)
expl_ostr_nomod_T2_model <- rma((dat$ostr_nomod_simple_main_T2)[dat$Model_delayed==1& !is.na(dat$ostr_mod_simple_main_T2)],(dat$ostr_nomod_simple_main_svar_T2)[dat$Model_delayed==1& !is.na(dat$ostr_mod_simple_main_T2)],method="REML",digits=4);
# tab_expl_ostr_nomod_T2_model <- write_res_tab("Ostracism T2 (no-moderator level)",expl_ostr_nomod_T2)
expl_mod_ostr_T2_model <- rma((dat$mod_ostr_simple_main_T2)[dat$Model_delayed==1],(dat$mod_ostr_simple_main_svar_T2)[dat$Model_delayed==1],method="REML",digits=4);
# tab_expl_mod_ostr_T2_model <- write_res_tab("Moderater T2 (ostracism level)",expl_mod_ostr_T2)
expl_mod_incl_T2_model <- rma((dat$mod_incl_simple_main_T2)[dat$Model_delayed==1],(dat$mod_incl_simple_main_svar_T2)[dat$Model_delayed==1],method="REML",digits=4);
# tab_expl_mod_incl_T2_model <- write_res_tab("Moderator T2 (inclusion level)",expl_mod_incl_T2) 
# rbind(tab_expl_ostr_mod_T1_model
#       ,tab_expl_ostr_nomod_T1_model
#       ,tab_expl_mod_ostr_T1_model
#       ,tab_expl_mod_incl_T1_model
#       ,tab_expl_ostr_mod_T2_model
#       ,tab_expl_ostr_nomod_T2_model
#       ,tab_expl_mod_ostr_T2_model
#       ,tab_expl_mod_incl_T2_model)

# Preparing the plots
# Traditional ostracism effect T1
dotplot_simple_effects <- c(expl_ostr_nomod_T1$b[1,1]
                            ,expl_ostr_nomod_T1_fund$b[1,1]
                            ,expl_ostr_nomod_T1_intra$b[1,1]
                            ,expl_ostr_nomod_T1_inter$b[1,1]
                            ,expl_ostr_nomod_T1_model$b[1,1])
dotplot_CI_low <- c(expl_ostr_nomod_T1$ci.lb
                    ,expl_ostr_nomod_T1_fund$ci.lb
                    ,expl_ostr_nomod_T1_intra$ci.lb
                    ,expl_ostr_nomod_T1_inter$ci.lb
                    ,expl_ostr_nomod_T1_model$ci.lb
)
dotplot_CI_up <- c(expl_ostr_nomod_T1$ci.ub
                   ,expl_ostr_nomod_T1_fund$ci.ub
                   ,expl_ostr_nomod_T1_intra$ci.ub
                   ,expl_ostr_nomod_T1_inter$ci.ub
                   ,expl_ostr_nomod_T1_model$ci.ub
)
dotplot_names <- c("All"
                   ,"Fundamental"
                   ,"Intrapersonal"
                   ,"Interpersonal"
                   ,"Model")
order <- c(1
           ,2
           ,3
           ,4
           ,5)
linetyp <- c(1
             ,1
             ,1
             ,1
             ,1)
numStud1 <- c(paste0("k=",expl_ostr_nomod_T1$k)
              ,paste0("k=",expl_ostr_nomod_T1_fund$k)
              ,paste0("k=",expl_ostr_nomod_T1_intra$k)
              ,paste0("k=",expl_ostr_nomod_T1_inter$k)
              ,paste0("k=",expl_ostr_nomod_T1_model$k))
dotplot_df1 <- data.frame(estimate=dotplot_simple_effects,ord=order,labels=dotplot_names,CIlow=dotplot_CI_low,CIup=dotplot_CI_up)
ostr_nomod_T1 <- xyplot(reorder(labels,sort(ord,decreasing=T))~estimate, data=dotplot_df1, xlim=c(-3,1.5),xlab="Effect size",ylab="Subset",panel=function(x,y){
  panel.xyplot(x, y, pch=18 ,cex=1.5, col = "black")
  panel.abline(v=0, col="black", lty=2)
  panel.segments(dotplot_df1$CIlow, as.numeric(y), dotplot_df1$CIup, as.numeric(y), lwd=2,lty=linetyp, col="black")
  panel.text(-2.5,5,labels=numStud1[1],cex=.75)
  panel.text(-2.5,4,labels=numStud1[2],cex=.75)
  panel.text(-2.5,3,labels=numStud1[3],cex=.75)
  panel.text(-2.5,2,labels=numStud1[4],cex=.75)
  panel.text(-2.5,1,labels=numStud1[5],cex=.75)
},main="(1) Traditional ostracism effect T1")

# Moderated ostracism effect T1 
dotplot_simple_effects <- c(expl_ostr_mod_T1$b[1,1]
                            ,expl_ostr_mod_T1_fund$b[1,1]
                            ,expl_ostr_mod_T1_intra$b[1,1]
                            ,expl_ostr_mod_T1_inter$b[1,1]
                            ,expl_ostr_mod_T1_model$b[1,1])
dotplot_CI_low <- c(expl_ostr_mod_T1$ci.lb
                    ,expl_ostr_mod_T1_fund$ci.lb
                    ,expl_ostr_mod_T1_intra$ci.lb
                    ,expl_ostr_mod_T1_inter$ci.lb
                    ,expl_ostr_mod_T1_model$ci.lb
)
dotplot_CI_up <- c(expl_ostr_mod_T1$ci.ub
                   ,expl_ostr_mod_T1_fund$ci.ub
                   ,expl_ostr_mod_T1_intra$ci.ub
                   ,expl_ostr_mod_T1_inter$ci.ub
                   ,expl_ostr_mod_T1_model$ci.ub
)
dotplot_names <- c("All"
                   ,"Fundamental"
                   ,"Intrapersonal"
                   ,"Interpersonal"
                   ,"Model")
order <- c(1
           ,2
           ,3
           ,4
           ,5)
linetyp <- c(1
             ,1
             ,1
             ,1
             ,1)
numStud2 <- c(paste0("k=",expl_ostr_mod_T1$k)
              ,paste0("k=",expl_ostr_mod_T1_fund$k)
              ,paste0("k=",expl_ostr_mod_T1_intra$k)
              ,paste0("k=",expl_ostr_mod_T1_inter$k)
              ,paste0("k=",expl_ostr_mod_T1_model$k))
dotplot_df2 <- data.frame(estimate=dotplot_simple_effects,ord=order,labels=dotplot_names,CIlow=dotplot_CI_low,CIup=dotplot_CI_up)
ostr_mod_T1 <- xyplot(reorder(labels,sort(ord,decreasing=T))~estimate, data=dotplot_df2, xlim=c(-3,1.5),xlab="Effect size",ylab="Subset",panel=function(x,y){
  panel.xyplot(x, y, pch=18,cex=1.5, col = "black")
  panel.abline(v=0, col="black", lty=2)
  panel.segments(dotplot_df2$CIlow, as.numeric(y), dotplot_df2$CIup, as.numeric(y), lwd=2,lty=linetyp, col="black")
  panel.text(1,5,labels=numStud2[1],cex=.75)
  panel.text(1,4,labels=numStud2[2],cex=.75)
  panel.text(1,3,labels=numStud2[3],cex=.75)
  panel.text(1,2,labels=numStud2[4],cex=.75)
  panel.text(1,1,labels=numStud2[5],cex=.75)
},main="(2) Moderated ostracism effect T1")

# Moderator effect ostracism level T1
dotplot_simple_effects <- c(expl_mod_ostr_T1$b[1,1]
                            ,expl_mod_ostr_T1_fund$b[1,1]
                            ,expl_mod_ostr_T1_intra$b[1,1]
                            ,expl_mod_ostr_T1_inter$b[1,1]
                            ,expl_mod_ostr_T1_model$b[1,1])
dotplot_CI_low <- c(expl_mod_ostr_T1$ci.lb
                    ,expl_mod_ostr_T1_fund$ci.lb
                    ,expl_mod_ostr_T1_intra$ci.lb
                    ,expl_mod_ostr_T1_inter$ci.lb
                    ,expl_mod_ostr_T1_model$ci.lb
)
dotplot_CI_up <- c(expl_mod_ostr_T1$ci.ub
                   ,expl_mod_ostr_T1_fund$ci.ub
                   ,expl_mod_ostr_T1_intra$ci.ub
                   ,expl_mod_ostr_T1_inter$ci.ub
                   ,expl_mod_ostr_T1_model$ci.ub
)
dotplot_names <- c("All"
                   ,"Fundamental"
                   ,"Intrapersonal"
                   ,"Interpersonal"
                   ,"Model")
order <- c(1
           ,2
           ,3
           ,4
           ,5)
linetyp <- c(1
             ,1
             ,1
             ,1
             ,1)
numStud3 <- c(paste0("k=",expl_mod_ostr_T1$k)
              ,paste0("k=",expl_mod_ostr_T1_fund$k)
              ,paste0("k=",expl_mod_ostr_T1_intra$k)
              ,paste0("k=",expl_mod_ostr_T1_inter$k)
              ,paste0("k=",expl_mod_ostr_T1_model$k))
dotplot_df3 <- data.frame(estimate=dotplot_simple_effects,ord=order,labels=dotplot_names,CIlow=dotplot_CI_low,CIup=dotplot_CI_up)
mod_ostr_T1 <- xyplot(reorder(labels,sort(ord,decreasing=T))~estimate, data=dotplot_df3, xlim=c(-3,1.5),xlab="Effect size",ylab="Subset",panel=function(x,y){
  panel.xyplot(x, y, pch=18 ,cex=1.5, col = "black")
  panel.abline(v=0, col="black", lty=2)
  panel.segments(dotplot_df3$CIlow, as.numeric(y), dotplot_df3$CIup, as.numeric(y), lwd=2,lty=linetyp, col="black")
  panel.text(-2.5,5,labels=numStud3[1],cex=.75)
  panel.text(-2.5,4,labels=numStud3[2],cex=.75)
  panel.text(-2.5,3,labels=numStud3[3],cex=.75)
  panel.text(-2.5,2,labels=numStud3[4],cex=.75)
  panel.text(-2.5,1,labels=numStud3[5],cex=.75)
},main="(5) Moderator effect within ostracism level T1")

# Moderator effect inclusion level T1
dotplot_simple_effects <- c(expl_mod_incl_T1$b[1,1]
                            ,expl_mod_incl_T1_fund$b[1,1]
                            ,expl_mod_incl_T1_intra$b[1,1]
                            ,expl_mod_incl_T1_inter$b[1,1]
                            ,expl_mod_incl_T1_model$b[1,1])
dotplot_CI_low <- c(expl_mod_incl_T1$ci.lb
                    ,expl_mod_incl_T1_fund$ci.lb
                    ,expl_mod_incl_T1_intra$ci.lb
                    ,expl_mod_incl_T1_inter$ci.lb
                    ,expl_mod_incl_T1_model$ci.lb
)
dotplot_CI_up <- c(expl_mod_incl_T1$ci.ub
                   ,expl_mod_incl_T1_fund$ci.ub
                   ,expl_mod_incl_T1_intra$ci.ub
                   ,expl_mod_incl_T1_inter$ci.ub
                   ,expl_mod_incl_T1_model$ci.ub
)
dotplot_names <- c("All"
                   ,"Fundamental"
                   ,"Intrapersonal"
                   ,"Interpersonal"
                   ,"Model")
order <- c(1
           ,2
           ,3
           ,4
           ,5)
linetyp <- c(1
             ,1
             ,1
             ,1
             ,1)
numStud4 <- c(paste0("k=",expl_mod_incl_T1$k)
              ,paste0("k=",expl_mod_incl_T1_fund$k)
              ,paste0("k=",expl_mod_incl_T1_intra$k)
              ,paste0("k=",expl_mod_incl_T1_inter$k)
              ,paste0("k=",expl_mod_incl_T1_model$k))
dotplot_df4 <- data.frame(estimate=dotplot_simple_effects,ord=order,labels=dotplot_names,CIlow=dotplot_CI_low,CIup=dotplot_CI_up)
mod_incl_T1 <- xyplot(reorder(labels,sort(ord,decreasing=T))~estimate, data=dotplot_df4, xlim=c(-3,1.5),xlab="Effect size",ylab="Subset",panel=function(x,y){
  panel.xyplot(x, y, pch=18 ,cex=1.5, col = "black")
  panel.abline(v=0, col="black", lty=2)
  panel.segments(dotplot_df4$CIlow, as.numeric(y), dotplot_df4$CIup, as.numeric(y), lwd=2,lty=linetyp, col="black")
  panel.text(1,5,labels=numStud4[1],cex=.75)
  panel.text(1,4,labels=numStud4[2],cex=.75)
  panel.text(1,3,labels=numStud4[3],cex=.75)
  panel.text(1,2,labels=numStud4[4],cex=.75)
  panel.text(1,1,labels=numStud4[5],cex=.75)
},main="(6) Moderator effect within inclusion level T1")

# Traditional ostracism effect T2
dotplot_simple_effects <- c(expl_ostr_nomod_T2$b[1,1]
                            ,expl_ostr_nomod_T2_fund$b[1,1]
                            ,expl_ostr_nomod_T2_intra$b[1,1]
                            ,expl_ostr_nomod_T2_inter$b[1,1]
                            ,expl_ostr_nomod_T2_model$b[1,1])
dotplot_CI_low <- c(expl_ostr_nomod_T2$ci.lb
                    ,expl_ostr_nomod_T2_fund$ci.lb
                    ,expl_ostr_nomod_T2_intra$ci.lb
                    ,expl_ostr_nomod_T2_inter$ci.lb
                    ,expl_ostr_nomod_T2_model$ci.lb
)
dotplot_CI_up <- c(expl_ostr_nomod_T2$ci.ub
                   ,expl_ostr_nomod_T2_fund$ci.ub
                   ,expl_ostr_nomod_T2_intra$ci.ub
                   ,expl_ostr_nomod_T2_inter$ci.ub
                   ,expl_ostr_nomod_T2_model$ci.ub
)
dotplot_names <- c("All"
                   ,"Fundamental"
                   ,"Intrapersonal"
                   ,"Interpersonal"
                   ,"Model")
order <- c(1
           ,2
           ,3
           ,4
           ,5)
linetyp <- c(1
             ,1
             ,1
             ,1
             ,1)
numStud5 <- c(paste0("k=",expl_ostr_nomod_T2$k)
              ,paste0("k=",expl_ostr_nomod_T2_fund$k)
              ,paste0("k=",expl_ostr_nomod_T2_intra$k)
              ,paste0("k=",expl_ostr_nomod_T2_inter$k)
              ,paste0("k=",expl_ostr_nomod_T2_model$k))
dotplot_df5 <- data.frame(estimate=dotplot_simple_effects,ord=order,labels=dotplot_names,CIlow=dotplot_CI_low,CIup=dotplot_CI_up)
ostr_nomod_T2 <- xyplot(reorder(labels,sort(ord,decreasing=T))~estimate, data=dotplot_df5, xlim=c(-3,1.5),xlab="Effect size",ylab="Subset",panel=function(x,y){
  panel.xyplot(x, y, pch=18,cex=1.5, col = "black")
  panel.abline(v=0, col="black", lty=2)
  panel.segments(dotplot_df5$CIlow, as.numeric(y), dotplot_df5$CIup, as.numeric(y), lwd=2,lty=linetyp, col="black")
  panel.text(-2.5,5,labels=numStud5[1],cex=.75)
  panel.text(-2.5,4,labels=numStud5[2],cex=.75)
  panel.text(-2.5,3,labels=numStud5[3],cex=.75)
  panel.text(-2.5,2,labels=numStud5[4],cex=.75)
  panel.text(-2.5,1,labels=numStud5[5],cex=.75)
},main="(3) Traditional ostracism effect T2")

# Moderated ostracism effect T2 
dotplot_simple_effects <- c(expl_ostr_mod_T2$b[1,1]
                            ,expl_ostr_mod_T2_fund$b[1,1]
                            ,expl_ostr_mod_T2_intra$b[1,1]
                            ,expl_ostr_mod_T2_inter$b[1,1]
                            ,expl_ostr_mod_T2_model$b[1,1])
dotplot_CI_low <- c(expl_ostr_mod_T2$ci.lb
                    ,expl_ostr_mod_T2_fund$ci.lb
                    ,expl_ostr_mod_T2_intra$ci.lb
                    ,expl_ostr_mod_T2_inter$ci.lb
                    ,expl_ostr_mod_T2_model$ci.lb
)
dotplot_CI_up <- c(expl_ostr_mod_T2$ci.ub
                   ,expl_ostr_mod_T2_fund$ci.ub
                   ,expl_ostr_mod_T2_intra$ci.ub
                   ,expl_ostr_mod_T2_inter$ci.ub
                   ,expl_ostr_mod_T2_model$ci.ub
)
dotplot_names <- c("All"
                   ,"Fundamental"
                   ,"Intrapersonal"
                   ,"Interpersonal"
                   ,"Model")
order <- c(1
           ,2
           ,3
           ,4
           ,5)
linetyp <- c(1
             ,1
             ,1
             ,1
             ,1)
numStud6 <- c(paste0("k=",expl_ostr_mod_T2$k)
              ,paste0("k=",expl_ostr_mod_T2_fund$k)
              ,paste0("k=",expl_ostr_mod_T2_intra$k)
              ,paste0("k=",expl_ostr_mod_T2_inter$k)
              ,paste0("k=",expl_ostr_mod_T2_model$k))
dotplot_df6 <- data.frame(estimate=dotplot_simple_effects,ord=order,labels=dotplot_names,CIlow=dotplot_CI_low,CIup=dotplot_CI_up)
ostr_mod_T2 <- xyplot(reorder(labels,sort(ord,decreasing=T))~estimate, data=dotplot_df6, xlim=c(-3,1.5),xlab="Effect size",ylab="Subset",panel=function(x,y){
  panel.xyplot(x, y, pch=18 ,cex=1.5, col = "black")
  panel.abline(v=0, col="black", lty=2)
  panel.segments(dotplot_df6$CIlow, as.numeric(y), dotplot_df6$CIup, as.numeric(y), lwd=2,lty=linetyp, col="black")
  panel.text(1,5,labels=numStud6[1],cex=.75)
  panel.text(1,4,labels=numStud6[2],cex=.75)
  panel.text(1,3,labels=numStud6[3],cex=.75)
  panel.text(1,2,labels=numStud6[4],cex=.75)
  panel.text(1,1,labels=numStud6[5],cex=.75)
},main="(4) Moderated ostracism effect T2")

# Moderator effect ostracism level T2
dotplot_simple_effects <- c(expl_mod_ostr_T2$b[1,1]
                            ,expl_mod_ostr_T2_fund$b[1,1]
                            ,expl_mod_ostr_T2_intra$b[1,1]
                            ,expl_mod_ostr_T2_inter$b[1,1]
                            ,expl_mod_ostr_T2_model$b[1,1])
dotplot_CI_low <- c(expl_mod_ostr_T2$ci.lb
                    ,expl_mod_ostr_T2_fund$ci.lb
                    ,expl_mod_ostr_T2_intra$ci.lb
                    ,expl_mod_ostr_T2_inter$ci.lb
                    ,expl_mod_ostr_T2_model$ci.lb
)
dotplot_CI_up <- c(expl_mod_ostr_T2$ci.ub
                   ,expl_mod_ostr_T2_fund$ci.ub
                   ,expl_mod_ostr_T2_intra$ci.ub
                   ,expl_mod_ostr_T2_inter$ci.ub
                   ,expl_mod_ostr_T2_model$ci.ub
)
dotplot_names <- c("All"
                   ,"Fundamental"
                   ,"Intrapersonal"
                   ,"Interpersonal"
                   ,"Model")
order <- c(1
           ,2
           ,3
           ,4
           ,5)
linetyp <- c(1
             ,1
             ,1
             ,1
             ,1)
numStud7 <- c(paste0("k=",expl_mod_ostr_T2$k)
              ,paste0("k=",expl_mod_ostr_T2_fund$k)
              ,paste0("k=",expl_mod_ostr_T2_intra$k)
              ,paste0("k=",expl_mod_ostr_T2_inter$k)
              ,paste0("k=",expl_mod_ostr_T2_model$k))
dotplot_df7 <- data.frame(estimate=dotplot_simple_effects,ord=order,labels=dotplot_names,CIlow=dotplot_CI_low,CIup=dotplot_CI_up)
mod_ostr_T2 <- xyplot(reorder(labels,sort(ord,decreasing=T))~estimate, data=dotplot_df7, xlim=c(-3,1.5),xlab="Effect size",ylab="Subset",panel=function(x,y){
  panel.xyplot(x, y, pch=18 ,cex=1.5, col = "black")
  panel.abline(v=0, col="black", lty=2)
  panel.segments(dotplot_df7$CIlow, as.numeric(y), dotplot_df7$CIup, as.numeric(y), lwd=2,lty=linetyp, col="black")
  panel.text(-2.5,5,labels=numStud7[1],cex=.75)
  panel.text(-2.5,4,labels=numStud7[2],cex=.75)
  panel.text(-2.5,3,labels=numStud7[3],cex=.75)
  panel.text(-2.5,2,labels=numStud7[4],cex=.75)
  panel.text(-2.5,1,labels=numStud7[5],cex=.75)
},main="(7) Moderator effect within ostracism level T2")

# Moderator effect inclusion level T2
dotplot_simple_effects <- c(expl_mod_incl_T2$b[1,1]
                            ,expl_mod_incl_T2_fund$b[1,1]
                            ,expl_mod_incl_T2_intra$b[1,1]
                            ,expl_mod_incl_T2_inter$b[1,1]
                            ,expl_mod_incl_T2_model$b[1,1])
dotplot_CI_low <- c(expl_mod_incl_T2$ci.lb
                    ,expl_mod_incl_T2_fund$ci.lb
                    ,expl_mod_incl_T2_intra$ci.lb
                    ,expl_mod_incl_T2_inter$ci.lb
                    ,expl_mod_incl_T2_model$ci.lb
)
dotplot_CI_up <- c(expl_mod_incl_T2$ci.ub
                   ,expl_mod_incl_T2_fund$ci.ub
                   ,expl_mod_incl_T2_intra$ci.ub
                   ,expl_mod_incl_T2_inter$ci.ub
                   ,expl_mod_incl_T2_model$ci.ub
)
dotplot_names <- c("All"
                   ,"Fundamental"
                   ,"Intrapersonal"
                   ,"Interpersonal"
                   ,"Model")
order <- c(1
           ,2
           ,3
           ,4
           ,5)
linetyp <- c(1
             ,1
             ,1
             ,1
             ,1)
numStud8 <- c(paste0("k=",expl_mod_incl_T2$k)
              ,paste0("k=",expl_mod_incl_T2_fund$k)
              ,paste0("k=",expl_mod_incl_T2_intra$k)
              ,paste0("k=",expl_mod_incl_T2_inter$k)
              ,paste0("k=",expl_mod_incl_T2_model$k))
dotplot_df8 <- data.frame(estimate=dotplot_simple_effects,ord=order,labels=dotplot_names,CIlow=dotplot_CI_low,CIup=dotplot_CI_up)
mod_incl_T2 <- xyplot(reorder(labels,sort(ord,decreasing=T))~estimate, data=dotplot_df8, xlim=c(-3,1.5),xlab="Effect size",ylab="Subset",panel=function(x,y){
  panel.xyplot(x, y, pch=18 ,cex=1.5, col = "black")
  panel.abline(v=0, col="black", lty=2)
  panel.segments(dotplot_df8$CIlow, as.numeric(y), dotplot_df8$CIup, as.numeric(y), lwd=2,lty=linetyp, col="black")
  panel.text(1,5,labels=numStud8[1],cex=.75)
  panel.text(1,4,labels=numStud8[2],cex=.75)
  panel.text(1,3,labels=numStud8[3],cex=.75)
  panel.text(1,2,labels=numStud8[4],cex=.75)
  panel.text(1,1,labels=numStud8[5],cex=.75)
},main="(8) Moderator effect within inclusion level T2")

# Actually plotting 
tiff(filename = 'Fig2.tiff', width=1000, height=1123, res = 300)
print(ostr_nomod_T1, position = c(0, .75, .5, 1), more = T)
print(ostr_nomod_T2, position = c(0, .5, .5, .75), more = T)
print(ostr_mod_T1, position = c(0.5, .75, 1, 1), more = T)
print(ostr_mod_T2, position = c(.5, .5, 1, .75), more = T)
print(mod_ostr_T1, position = c(0, .25, .5, .5), more = T)
print(mod_ostr_T2, position = c(0, 0, .5, .25), more = T)
print(mod_incl_T1, position = c(.5, .25, 1, .5), more = T)
print(mod_incl_T2, position = c(.5, 0, 1, .25), more = F)
dev.off()
