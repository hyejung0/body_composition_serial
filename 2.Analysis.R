#2024 Dec 11 (Wed)

#Well I lost my quarto file
#that I worked on for 9 hrs yesterday
#so I will make a backup copy here.


rm(list=ls())
library(data.table)
# data --------------------------------------------------------------------

dat<-readRDS("clean_dat_new_BMI.rds")


dat[,PatientID_fct:=factor(PatientID)] #Must change to factor to use in my regression as random effect. 


#save column names
#change the column names
outcomes<-c("new_OS_time","new_OS_ind","new_OS_time_mth","new_OS_time_yr")
covar_cols<-c("ageAtMetDiag"="Age at Metastatic Diagnosis",
              "gender"="Gender",
              "smoking"="Smoking",
              "raceEthnicity"="Race/Ethnicity",
              "trt_type"="Treatment Type",
              "bmi"="BMI (kg/m2)",
              "bmi_cateogrical_new"="BMI (categorical)",
              "tumor_hist"="Tumor Histology",
              "Gene"="Mutation Status" )

body_cols<-colnames(dat)[!colnames(dat) %in% c("mrn","PatientID","PatientID_fct","is_baseline","new_study_date","StudyDate","times","times_mth","times_yr","earliest_studydate","metastaticDate",names(covar_cols),outcomes)]
body_cols_names<-
  sapply(body_cols,function(x){
    
    pieces_x<-strsplit(x,";")[[1]]
    
    if("cross_sectional_area_cm2" %in% pieces_x){
      paste(pieces_x[2],"cm2",sep="_")
      
    }else if("volume_cm3" %in% pieces_x){
      paste(pieces_x[2],"cm3",sep="_")
      
    }else if("HU_mean" %in% pieces_x){
      paste(pieces_x[2],"mean_HU",sep="_")
      
    }else if(stringr::str_detect(string = x,pattern = "whole")){
      pieces_x2<-strsplit(x,"_")[[1]]
      
      paste(pieces_x2[1],"whole_body_cm3",sep="_")
    }
    
  })

setnames(dat,
         old=names(body_cols_names),
         new=body_cols_names)


# Separate Cox modeling ---------------------------------------------------

#generate time-dependent dataset


dat[,is_baseline:=0]
dat[StudyDate==earliest_studydate,is_baseline:=1]
baseline_dat<-dat[is_baseline==1,]


baseline_dat2<-copy(baseline_dat)
grab_thes<-colnames(baseline_dat2)[!colnames(baseline_dat2) %in% names(body_cols_names)]
baseline_dat2<-baseline_dat2[,.SD,.SDcols = grab_thes]
tmp2<-tmerge(baseline_dat2,
             baseline_dat2,
             id=PatientID,endpt = event(new_OS_time_yr, new_OS_ind))



#Yes it works.
#So just tmerge first and replace with rpoper values.
tmp3<-tmerge(
  tmp2, 
  dat, 
  id=PatientID, 
  x1 = tdc(times_yr, get(body_cols_names[1], dat)),
  x2 = tdc(times_yr, get(body_cols_names[2], dat)),
  x3 = tdc(times_yr, get(body_cols_names[3], dat)),
  x4 = tdc(times_yr, get(body_cols_names[4], dat)),
  x5 = tdc(times_yr, get(body_cols_names[5], dat)),
  x6 = tdc(times_yr, get(body_cols_names[6], dat)),
  x7 = tdc(times_yr, get(body_cols_names[7], dat)),
  x8 = tdc(times_yr, get(body_cols_names[8], dat)),
  x9 = tdc(times_yr, get(body_cols_names[9], dat)),
  x10 = tdc(times_yr, get(body_cols_names[10], dat)),
  x11 = tdc(times_yr, get(body_cols_names[11], dat)),
  x12 = tdc(times_yr, get(body_cols_names[12], dat)),
  x13 = tdc(times_yr, get(body_cols_names[13], dat)),
  x14 = tdc(times_yr, get(body_cols_names[14], dat)),
  x15 = tdc(times_yr, get(body_cols_names[15], dat)),
  x16 = tdc(times_yr, get(body_cols_names[16], dat))
)

tmp3<-data.table(tmp3)

setnames(tmp3,
         old=paste("x",1:16,sep = ""),
         new=body_cols_names)
tmp3[,.(PatientID,tstart,tstop,IMAT_cm2)]
#Replace with correct values
(tmp3$PatientID==dat$PatientID) |> table()
(tmp3$tstart==dat$times_yr) |> table()
#The two data tables are in the same order. 
#Just replace the values

for( this.body in body_cols_names){
  tmp3[,c(this.body):=dat[,get(this.body)]]
}


tmp3[,.(PatientID,tstart,tstop,SAT_whole_body_cm3)]
dat[,.(PatientID,times_yr,SAT_whole_body_cm3)]
#Fit time-dependent cox model
cox_naive<-
  lapply(body_cols_names,function(x){
    my_formula<-paste("Surv(tstart,tstop,endpt)~","`",x,"`","+ageAtMetDiag + gender + smoking + raceEthnicity + tumor_hist + Gene + bmi_cateogrical_new", sep="")
    my_formula<-as.formula(my_formula)
    coxph(my_formula, data=tmp3, x=T)
  })
names(cox_naive)<-body_cols_names
cox_naive$IMAT_cm2

#return the estimate
source("coxph.ci.coef.R")
cox_naive_HR<-lapply(cox_naive,coxph.ci.coef, digits=3)
cox_naive_HR
cox_naive_HR<-do.call(rbind,cox_naive_HR)
cox_naive_HR<-data.table(cox_naive_HR)
setnames(cox_naive_HR,
         old=colnames(cox_naive_HR),
         new=c("Body Composition",
               "Age",
               "Gender (Female)",
               "Smoking (No history of smoking)",
               "Race/Ethnicity (Other/Missing)",
               "Tumor histology (Adenocarcinoma)",
               "Tumor histology (Squamous cell carcinoma, NOS)",
               "Genetic Mutation (mutated)",
               "BMI (Overweight (>= 25kg/m2 & < 30kg/m^2))",
               "BMI (Obesity class I/II (>= 30kg/m2 & < 40kg/m2))"))
cox_naive_HR<-data.table(Outcome=names(cox_naive),
                         cox_naive_HR)
cox_naive_HR



# separate longitudinal model ---------------------------------------------



library(nlme)

rm(list=ls())
dat<-readRDS("clean_dat_new_BMI.rds")


dat[,PatientID_fct:=factor(PatientID)] #Must change to factor to use in my regression as random effect. 


#save column names
#change the column names
outcomes<-c("new_OS_time","new_OS_ind","new_OS_time_mth","new_OS_time_yr")
covar_cols<-c("ageAtMetDiag"="Age at Metastatic Diagnosis",
              "gender"="Gender",
              "smoking"="Smoking",
              "raceEthnicity"="Race/Ethnicity",
              "trt_type"="Treatment Type",
              "bmi"="BMI (kg/m2)",
              "bmi_cateogrical_new"="BMI (categorical)",
              "tumor_hist"="Tumor Histology",
              "Gene"="Mutation Status" )

body_cols<-colnames(dat)[!colnames(dat) %in% c("mrn","PatientID","PatientID_fct","is_baseline","new_study_date","StudyDate","times","times_mth","times_yr","earliest_studydate","metastaticDate",names(covar_cols),outcomes)]

#We get reid of all special characters
body_cols_names<-
  sapply(body_cols,function(x){
    
    pieces_x<-strsplit(x,";")[[1]]
    
    if("cross_sectional_area_cm2" %in% pieces_x){
      paste(pieces_x[2],"cm2",sep="_")
      
    }else if("volume_cm3" %in% pieces_x){
      paste(pieces_x[2],"cm3",sep="_")
      
    }else if("HU_mean" %in% pieces_x){
      paste(pieces_x[2],"mean_HU",sep="_")
      
    }else if(stringr::str_detect(string = x,pattern = "whole")){
      pieces_x2<-strsplit(x,"_")[[1]]
      
      paste(pieces_x2[1],"whole_body_cm3",sep="_")
    }
    
  })

setnames(dat,
         old=names(body_cols_names),
         new=body_cols_names)


#Below allows for debugging in Quarto.
myCatch <- function(# The expression to execute.
  expr,
  # Further arguments to tryCatch().
  ...,
  # User-defined function to extract diagnostic info from
  # warning object, based on output that resulted from expr.
  custom_fun = function(result, w){return(w)}) {
  ######################
  ## Default Settings ##
  ######################
  
  # Defaults to NULL results and empty list of diagnostics.
  DEFAULT_RESULTS <- NULL
  DEFAULT_DIAGNOSTICS <- NULL
  
  # Defaults to standard R error message, rather than a ponderous traceback
  # through the error handling stacks themselves; also returns the error object
  # itself as the results.
  DEFAULT_ERROR <- function(e){
    message("Error in ", deparse(e$call), " : ", e$message)
    return(e)
  }
  
  
  ################
  ## Initialize ##
  ################
  
  # Initialize output to default settings.
  res <- DEFAULT_RESULTS
  diag <- DEFAULT_DIAGNOSTICS
  err <- DEFAULT_ERROR
  
  # Adjust error handling if specified by user.
  if("error" %in% names(list(...))) {
    err <- list(...)$error
  }
  
  
  #######################
  ## Handle Expression ##
  #######################
  
  res <- tryCatch(
    expr = {
      withCallingHandlers(
        expr = expr,
        # If expression throws a warning, record diagnostics without halting,
        # so as to store the result of the expression.
        warning = function(w){
          parent <- parent.env(environment())
          parent$diag <- w
        }
      )
    },
    error = err,
    ...
  )
  
  
  ############
  ## Output ##
  ############
  
  # Package the results as desired.
  return(list(result = res,
              diagnostics = custom_fun(res, diag)))
}



#Let's fit linear term on time and compare that with flexible model.
lme_out_0<-lapply(body_cols_names,function(this_body){ #this fits times_yr as is
  fixed_lme_formula<-paste("`",this_body,"` ~ ageAtMetDiag + gender + smoking + raceEthnicity + tumor_hist + Gene + bmi_cateogrical_new+times_yr", sep="")
  fixed_lme_formula<-as.formula(fixed_lme_formula)
  myCatch(lme(fixed_lme_formula, 
              random = ~ times_yr | PatientID_fct, 
              method = "ML", 
              data = dat,
              control =list(msMaxIter = 1000, msMaxEval = 1000)))
}) #error

lme_out_0<-lapply(lme_out_0,"[[","result")
names(lme_out_0)<-body_cols_names[names(lme_out_0)]

error_lme_out_0<-
sapply(lme_out_0,function(x){ #Return TRUE if the fit had an error
  class_of_fit<-class(x) 
  if(any(class_of_fit=="error")){
    TRUE
  }else{
    FALSE
  }
})
which(error_lme_out_0==TRUE)
body_cols_names[error_lme_out_0] #error fit on IMAT_whole_body_cm3

#Try fitting simplified random effect.
#Below fixes the same covariances off diagonal, but different variance in the diagonal. 
#I purposefully chose pdDiag b/c jointModel() can hadnle it fine. page 60 of textbook by Pizopoulos.
#Regarding the covariance matrix of the random eﬀects, 
#by default, jointModel() assumes it to be unstructured and estimates all of its free parameters. 
#However, it also allows for a diagonal covariance matrix, which can be speciﬁed using function pdDiag() in the random argument of lme().

test<-lme(IMAT_whole_body_cm3~ageAtMetDiag + gender + smoking + raceEthnicity + tumor_hist + Gene + bmi_cateogrical_new+times_yr, 
           random=list(PatientID_fct = pdDiag(~ 1 + times_yr)), #this is a diagonal correlation structure, removing relationship between random intercept and slope.
           method = "ML", 
           data = dat,
           control =list(msMaxIter = 1000, msMaxEval = 1000))
summary(test)
#converged just fine. replace it.
lme_out_0$IMAT_whole_body_cm3<-test

saveRDS(lme_out_0,"lme_out_0.rds")


#Try plotting one fit
#For that one we will just fit some splines model. 
ggplot(data=dat)+
  geom_point(aes(x=times_yr, y=IMAT_whole_body_cm3,group = PatientID)) +
  geom_line(aes(x=times_yr, y=IMAT_whole_body_cm3,group = PatientID), alpha=0.2)+
  stat_smooth(
    data=dat,
    mapping=aes(x=times_yr, y=IMAT_whole_body_cm3),
    geom="smooth",
    method="lm",
    formula=y~I(x^3), 
    col="red")+
  xlab("Time (year)")




#Fit non-linear model 
lme_out_poly3<-lapply(body_cols_names,function(this_body){
  fixed_lme_formula<-paste("`",this_body,"` ~ ageAtMetDiag + gender + smoking + raceEthnicity + tumor_hist + Gene + bmi_cateogrical_new+times_yr+times_yr+I(times_yr^2)+I(times_yr^3)", sep="")
  fixed_lme_formula<-as.formula(fixed_lme_formula)
  myCatch(lme(fixed_lme_formula, 
              random = ~ times_yr+I(times_yr^2)+I(times_yr^3) | PatientID_fct, 
              method = "ML", 
              data = dat,
              control =list(msMaxIter = 1000, msMaxEval = 1000)))
})

lme_out_poly3<-lapply(lme_out_poly3,"[[","result")
names(lme_out_poly3)<-body_cols_names[names(lme_out_poly3)]

error_lme_out_poly3<-
  sapply(lme_out_poly3,function(x){ #Return TRUE if the fit had an error
    class_of_fit<-class(x) 
    if(any(class_of_fit=="error")){
      TRUE
    }else{
      FALSE
    }
  })
body_cols_names[error_lme_out_poly3] #error fit IMAT_whole_cm3

#fit the simplified model again
test<-lme(IMAT_whole_body_cm3~ageAtMetDiag + gender + smoking + raceEthnicity + tumor_hist + Gene + bmi_cateogrical_new+times_yr+I(times_yr^2)+I(times_yr^3), 
          random=list(PatientID_fct = pdDiag(~ 1 + times_yr+I(times_yr^2)+I(times_yr^3))), #this is a diagonal correlation structure, removing relationship between random intercept and slope.
          method = "ML", 
          data = dat,
          control =list(msMaxIter = 1000, msMaxEval = 1000))
summary(test) #fitted fine
#replace the fit with this one
lme_out_poly3$IMAT_whole_body_cm3<-test
saveRDS(lme_out_poly3,"lme_out_poly3.rds")



tmpdat<-data.table(
  dat[,.(times_yr,PatientID,SKM_whole_body_cm3)],
  data.table(lme_out_poly3$SKM_whole_body_cm3$fitted)
)


ggplot(data=tmpdat)+
  geom_point(aes(x=times_yr, y=SKM_whole_body_cm3,group = PatientID)) +
  geom_line(aes(x=times_yr, y=SKM_whole_body_cm3,group = PatientID), alpha=0.2)+
  # geom_line(aes(x=times_yr, y=fixed,group = PatientID), alpha=0.2, col="red")+
  # geom_line(aes(x=times_yr, y=PatientID_fct,group = PatientID), alpha=0.2, col="red")+
  stat_smooth(
    data=tmpdat,
    mapping=aes(x=times_yr, y=SKM_whole_body_cm3),
    geom="smooth",
    method="lm",
    formula=y~ns(x, df=3), 
    col="red")+  
  stat_smooth(
    data=tmpdat,
    mapping=aes(x=times_yr, y=PatientID_fct),
    geom="smooth",
    method="lm",
    formula=y~ns(x, df=3), 
    col="blue")+  
  stat_smooth(
    data=tmpdat,
    mapping=aes(x=times_yr, y=fixed),
    geom="smooth",
    method="lm",
    formula=y~ns(x, df=3), 
    col="orange")+ 
  xlab("Time (year)")


anova_out<-list()
for(i in 1:length(lme_out_0)){
  anova_out[[i]]<-anova(lme_out_0[[i]], lme_out_poly3[[i]])
}
names(anova_out)<-body_cols_names
anova.p<-sapply(anova_out,function(x)x[2,"p-value"]) 
anova.p.bh<-p.adjust(anova.p,method = "BH")
sum(anova.p.bh<=0.05) #15. all but one has significantly different fit. 
#which one is insignificant?
names(anova.p.bh)[anova.p.bh>0.05] #IMAT_cm3
anova_out$IMAT_cm3
anova.p.bh['IMAT_cm3']


tmpdat<-data.table(
  dat[,.(times_yr,PatientID,IMAT_cm3)],
  fixed_nonlinear=lme_out_poly3$IMAT_cm3$fitted[,"fixed"],
  PatientID_fct_nonlinear=lme_out_poly3$IMAT_cm3$fitted[,"PatientID_fct"],
  data.table(lme_out_0$IMAT_cm3$fitted)
)


ggplot(data=tmpdat)+
  geom_point(aes(x=times_yr, y=IMAT_cm3,group = PatientID)) +
  geom_line(aes(x=times_yr, y=IMAT_cm3,group = PatientID), alpha=0.2)+
  # geom_line(aes(x=times_yr, y=PatientID_fct,group = PatientID), alpha=0.4, col="red")+
  # geom_line(aes(x=times_yr, y=PatientID_fct_nonlinear,group = PatientID), alpha=0.4, col="red")+
  geom_line(aes(x=times_yr, y=fixed_nonlinear,group = PatientID), alpha=0.4, col="red")+
  xlab("Time (year)")

summary(lme_out_poly3[[1]])

#To understand which term explains the variation in the outcome
#try fitting the orthogonal polynomial. 
#Currently, I fitted polynomial on the raw times_yr value.

lme_out_poly3_ortho<-lapply(body_cols_names,function(this_body){
  fixed_lme_formula<-paste("`",this_body,"` ~ ageAtMetDiag + gender + smoking + raceEthnicity + tumor_hist + Gene + bmi_cateogrical_new+stats::poly(times_yr,3)", sep="")
  fixed_lme_formula<-as.formula(fixed_lme_formula)
  myCatch(lme(fixed_lme_formula, 
              random = ~ stats::poly(times_yr,3)| PatientID_fct, 
              method = "ML", 
              data = dat,
              control =list(msMaxIter = 1000, msMaxEval = 1000)))
})











#What about perfect correlation?
lapply(lme_out,function(x){
  tmp<-VarCorr(x)
  tmp[,"Corr"]
})
#Perfect correlation for IMAT_whole_cm3, SAT_whole_cm3
#Obtaining a random effect correlation estimate of +1 or -1 means that the optimization algorithm hit "a boundary": correlations cannot be higher than +1 or lower than -1
#this usually means that there are not enough data to estimate all the parameters reliably. Matuschek et al. 2017 say that in this situation the power can be compromised.



#Change the covariance structure to involve less number of paramters.
getVarCov(lme_out$IMAT_cm2, individuals = 1,type="marginal") #this is unstructured. That's what we have for now
#Compound symmetric (only random intercept)
fixed_lme_formula<-paste("`","IMAT_cm2","` ~ ageAtMetDiag + gender + smoking + raceEthnicity + tumor_hist + Gene + bmi_cateogrical_new+times_yr", sep="")
fixed_lme_formula<-as.formula(fixed_lme_formula)
test<-lme(fixed_lme_formula, 
            random =~1|PatientID_fct, 
            method = "ML", 
            data = dat,
            control =list(msMaxIter = 1000, msMaxEval = 1000))
summary(test)
getVarCov(test, individuals = 1,type="marginal") #this is unstructured. 



#Below fixes the same covariances off diagonal, but different variance in the diagonal. I purposefully chose pdDiag b/c jointModel() can hadnle it fine. page 60 of textbook by Pizopoulos.
#Regarding the covariance matrix of the random eﬀects, by default jointModel() assumes it to be unstructured and estimates all of its free parameters. 
#However, it also allows for a diagonal covariance matrix, which can be speciﬁed using function pdDiag() in the random argument of lme().

test2<-lme(fixed_lme_formula, 
          random=list(PatientID_fct = pdDiag(~ 1 + times_yr)), #this is a diagonal correlation structure, removing relationship between random intercept and slope.
          method = "ML", 
          data = dat,
          control =list(msMaxIter = 1000, msMaxEval = 1000))
summary(test2)
getVarCov(test2, individuals = 1,type="marginal") #this is unstructured. 
#I like this more than compound symetric.. I think.. 
#Let's go with this. 


#IMAT_whole_cm3, 
fixed_lme_formula<-paste("IMAT_whole_body_cm3 ~ ageAtMetDiag + gender + smoking + raceEthnicity + tumor_hist + Gene + bmi_cateogrical_new+splines::bs(times_yr,knots=0.25,degree=1)", sep="")
fixed_lme_formula<-as.formula(fixed_lme_formula)
out<-lme(fixed_lme_formula, 
         random = list(PatientID_fct = pdDiag(~ 1 + splines::bs(times_yr,knots=0.25,degree=1) )), 
         method = "ML", 
        data = dat,
          control =list(msMaxIter = 1000, msMaxEval = 1000))
summary(out)
#No more perfect correlation but the redicual deviance is high.
#What if I change to compound symmetric?
out2<-lme(fixed_lme_formula, 
         random =~1|PatientID_fct, 
         method = "ML", 
         data = dat,
         control =list(msMaxIter = 1000, msMaxEval = 1000))
summary(out2) #STill high residual standard deviation. 

lme_out$IMAT_whole_body_cm3<-out




#SAT_whole_cm3
fixed_lme_formula<-paste("SAT_whole_body_cm3 ~ ageAtMetDiag + gender + smoking + raceEthnicity + tumor_hist + Gene + bmi_cateogrical_new+times_yr", sep="")
fixed_lme_formula<-as.formula(fixed_lme_formula)
out<-lme(fixed_lme_formula, 
         random = list(PatientID_fct = pdDiag(~ 1 + times_yr)), 
         method = "ML", 
         data = dat,
         control =list(msMaxIter = 1000, msMaxEval = 1000))
summary(out)
lme_out$SAT_whole_body_cm3<-out


saveRDS(lme_out,"lme_out_clean.rds") #all models have good convergence of paramters




# joint modeling ----------------------------------------------------------

## JointModel function in the JM package fits joint models. It requires the fitted objects of 
## separate longitudinal and survival (Cox proportional hazards model with baseline covariate) models. 

# RRT

library(JM)
library(nlme)
library(survival)


#JointModel has error 
#with time issue for piecewise-PH-aGH. 
#Scale to monthly or yearly
#https://github.com/drizopoulos/JM/issues/10
#I had an issue with that so that's why I creaed times_yr and new_OS_times_yr.


baseline_dat<-dat[is_baseline==1,]

#Fit baseline cox model. Need it to fit JM.
cox.fit <- coxph(Surv(new_OS_time,new_OS_ind) ~ ageAtMetDiag + gender + smoking + raceEthnicity + tumor_hist + Gene + bmi_cateogrical_new, data=baseline_dat, x=TRUE)


source("debug_JM.R") #I tried to fit method = "Cox-PH-aGH" as it is most unstructured baseline risk.
#But I realized that there's an error with the current function.
#jointModel() calls initial.surv(). initial.surv() uses basehaz().
#basehaz(data,FALSE). is what's written by the package developer.
#but basehaz(fit, newdata, centered=TRUE) <- this is form. So the functio thinks I'm inputting newdata=FALSE.
#So I fixed it. 
#That's what I source.

joint_out1<-list()
for(i in 1:length(body_cols_names)){
  set.seed(55)
  joint_out1[[i]]<-jointModel_debug_Hyejung(
    lme_out[[i]], 
    cox.fit, 
    timeVar = "times_yr", 
    method = "Cox-PH-aGH",   #nspeciﬁedbaseline risk function. This assumption turns out to b e equivalent toassuming that h0(·) is discrete with point-masses at the unique eventtime
    verbose = T,
    control=c(iter.EM=3000)
             )
} #stops at 1st list due to missing value where TRUE/FALSE needed.


#change to Cox-PH-GH
joint_out1<-list()
for(i in 1:length(body_cols_names)){
  set.seed(55)
  joint_out1[[i]]<-jointModel_debug_Hyejung(
    lme_out[[i]], 
    cox.fit, 
    timeVar = "times_yr", 
    method = "Cox-PH-GH",   #nspeciﬁedbaseline risk function. This assumption turns out to b e equivalent toassuming that h0(·) is discrete with point-masses at the unique eventtime
    verbose = T,
    control=c(iter.EM=3000)
  )
} #stops at 1st list due to missing value where TRUE/FALSE needed.


#Go to piecewise-PH-aGH
#We will fix number of intenal knot location to be == 12. 
#That allows for different baseline hazard 
joint_out1<-list()
for(i in 1:length(body_cols_names)){
  set.seed(55)
  joint_out1[[i]]<-jointModel(
    lme_out[[i]], 
    cox.fit, 
    timeVar = "times_yr", 
    method = "piecewise-PH-aGH",   #nspeciﬁedbaseline risk function. This assumption turns out to b e equivalent toassuming that h0(·) is discrete with point-masses at the unique eventtime
    verbose = T,
    control=c(iter.EM=3000, lng.in.kn=12)
  )
} 
#error at i=8.
#Error in if (t1 || t2) { : missing value where TRUE/FALSE needed
#Reduce lng.in.kn=12 to smaller
length(joint_out1)
set.seed(55)
tmp<-jointModel(
  lme_out[[8]], 
  cox.fit, 
  timeVar = "times_yr", 
  method = "piecewise-PH-aGH",   #nspeciﬁedbaseline risk function. This assumption turns out to b e equivalent toassuming that h0(·) is discrete with point-masses at the unique eventtime
  verbose = T,
  control=c(iter.EM=3000, lng.in.kn=5)
)
#Doesn't work until lng.in.kn=5
#still getting
#Error in if (t1 || t2) { : missing value where TRUE/FALSE needed
#change the lme model to simplify random effect structure. SEe if that helps.
formula(lme_out[[8]])
out<-lme(formula(lme_out[[8]]), 
         random = list(PatientID_fct = pdDiag(~ 1 + times_yr)), 
         method = "ML", 
         data = dat,
         control =list(msMaxIter = 1000, msMaxEval = 1000))
summary(out)
set.seed(55)
tmp<-jointModel(
  out, 
  cox.fit, 
  timeVar = "times_yr", 
  method = "piecewise-PH-aGH",   #nspeciﬁedbaseline risk function. This assumption turns out to b e equivalent toassuming that h0(·) is discrete with point-masses at the unique eventtime
  verbose = T,
  control=c(iter.EM=3000, lng.in.kn=16)
)#NaN for D. in the verbose...
#Check for multicolinearity in my data
library(car)
vif(out) #Low multicolinearity. 
lme_out[[8]]$data

set.seed(55)
tmp<-jointModel(
  out, 
  cox.fit, 
  timeVar = "times_yr", 
  method = "piecewise-PH-aGH",   #nspeciﬁedbaseline risk function. This assumption turns out to b e equivalent toassuming that h0(·) is discrete with point-masses at the unique eventtime
  verbose = T,
  control=c(iter.EM=3000, lng.in.kn=10,
            optim = list(
              method = "BFGS",
              maxit = 1000,
              parscale = rep(0.01, length(coef(out)) + length(coef(cox.fit)))
            )))

#change to compound symmetric.
out2<-lme(formula(lme_out[[8]]), 
         random = ~1|PatientID_fct, 
         method = "ML", 
         data = dat,
         control =list(msMaxIter = 1000, msMaxEval = 1000))

tmp<-jointModel(
  out2, 
  cox.fit, 
  timeVar = "times_yr", 
  method = "piecewise-PH-aGH",   #nspeciﬁedbaseline risk function. This assumption turns out to b e equivalent toassuming that h0(·) is discrete with point-masses at the unique eventtime
  verbose = T,
  control=c(iter.EM=3000, lng.in.kn=12))#No... it doens' twork.


#Change to piecewise-PH-GH
tmp<-jointModel(
  lme_out[[8]], 
  cox.fit, 
  timeVar = "times_yr", 
  method = "piecewise-PH-GH",   #nspeciﬁedbaseline risk function. This assumption turns out to b e equivalent toassuming that h0(·) is discrete with point-masses at the unique eventtime
  verbose = T,
  control=c(iter.EM=3000, lng.in.kn=6))#No... it doens't work.
#lng.in.kn = 12,10,6 didn't work.

tmp<-jointModel(
  out,
  cox.fit, 
  timeVar = "times_yr", 
  method = "piecewise-PH-GH",   #nspeciﬁedbaseline risk function. This assumption turns out to b e equivalent toassuming that h0(·) is discrete with point-masses at the unique eventtime
  verbose = T,
  control=c(iter.EM=3000, lng.in.kn=6))#No... it doens' twork.
#lng.in.kn = 12,10,6 didn't work.



#In the textbook, they speak about when log-likelihood=0
#Look at coefficient of estimate in LME output.
#For them, their effect was too large for age
#So they scaled AGE
lme_out[[8]] |> summary()
#we don't have such issue.
#Take out fixed effect.
summary(cox.fit)
#I think we can take out GeneMutation
#That might be hard to get
#Expensive and invasive
cox.fit2<-coxph(formula = Surv(new_OS_time, new_OS_ind) ~ ageAtMetDiag + 
                  gender + smoking + raceEthnicity + tumor_hist + bmi_cateogrical_new, 
                data = baseline_dat, x = TRUE)
anova(cox.fit2,cox.fit) #Not different at all

#What about tumor histology?
cox.fit3<-coxph(formula = Surv(new_OS_time, new_OS_ind) ~ ageAtMetDiag + 
                  gender + smoking + raceEthnicity + bmi_cateogrical_new, 
                data = baseline_dat, x = TRUE)
anova(cox.fit3,cox.fit) #Not different at all....
#I could take both out.
#what about in LME model?

lme_out_8_1<-lme(SKM_cm3 ~ ageAtMetDiag + gender + smoking + raceEthnicity  +tumor_hist+ bmi_cateogrical_new + times_yr, 
                 random = ~times_yr|PatientID_fct, 
                 method = "ML", 
                 data = dat,
                 control =list(msMaxIter = 1000, msMaxEval = 1000))
lme_out_8_2<-lme(SKM_cm3 ~ ageAtMetDiag + gender + smoking + raceEthnicity+ bmi_cateogrical_new + times_yr, 
                 random = ~times_yr|PatientID_fct, 
                 method = "ML", 
                 data = dat,
                 control =list(msMaxIter = 1000, msMaxEval = 1000))
anova(lme_out_8_1,lme_out[[8]]) #Not different at all....
anova(lme_out_8_2,lme_out[[8]]) #Not different at all....

#Well, we could take out both Gene and TumorHist if wanted to . Let's give them a try. 

out<-jointModel(
  lme_out_8_1, 
  cox.fit2, 
  timeVar = "times_yr", 
  method = "piecewise-PH-aGH",   #nspeciﬁedbaseline risk function. This assumption turns out to b e equivalent toassuming that h0(·) is discrete with point-masses at the unique eventtime
  verbose = T,
  control=c(iter.EM=3000, lng.in.kn=6)
)#no...


out<-jointModel(
  lme_out_8_2, 
  cox.fit3, 
  timeVar = "times_yr", 
  method = "piecewise-PH-aGH",   #nspeciﬁedbaseline risk function. This assumption turns out to b e equivalent toassuming that h0(·) is discrete with point-masses at the unique eventtime
  verbose = T,
  control=c(iter.EM=3000, lng.in.kn=7)
)#No...





lme_out_8_3<-lme(SKM_cm3 ~ ageAtMetDiag + gender + smoking + raceEthnicity  +tumor_hist+ bmi_cateogrical_new + times_yr, 
                 random = list(PatientID_fct = pdDiag(~ 1 + times_yr)), 
                 method = "ML", 
                 data = dat,
                 control =list(msMaxIter = 1000, msMaxEval = 1000))
lme_out_8_4<-lme(SKM_cm3 ~ ageAtMetDiag + gender + smoking + raceEthnicity+ bmi_cateogrical_new + times_yr, 
                 random = list(PatientID_fct = pdDiag(~ 1 + times_yr)), 
                 method = "ML", 
                 data = dat,
                 control =list(msMaxIter = 1000, msMaxEval = 1000))

out<-jointModel(
  lme_out_8_3, 
  cox.fit2, 
  timeVar = "times_yr", 
  method = "piecewise-PH-aGH",   #nspeciﬁedbaseline risk function. This assumption turns out to b e equivalent toassuming that h0(·) is discrete with point-masses at the unique eventtime
  verbose = T,
  control=c(iter.EM=3000, lng.in.kn=6)
)#no...



out<-jointModel(
  lme_out_8_4, 
  cox.fit3, 
  timeVar = "times_yr", 
  method = "piecewise-PH-aGH",   #nspeciﬁedbaseline risk function. This assumption turns out to b e equivalent toassuming that h0(·) is discrete with point-masses at the unique eventtime
  verbose = T,
  control=c(iter.EM=3000, lng.in.kn=6)
)#no...


#WHY..
#What's going on

ggplot(data=dat)+
  geom_point(aes(x=times_yr, y=SKM_cm3,group = PatientID)) +
  geom_line(aes(x=times_yr, y=SKM_cm3,group = PatientID), alpha=0.2)+
  stat_smooth(
    data=dat,
    mapping=aes(x=times_yr, y=SKM_cm3),
    geom="smooth",
    method="lm",
    formula=y~ splines::ns(x,df=3), 
    col="red")+
  xlab("Time (year)")



out<-jointModel(
  lme_out_8_1, 
  cox.fit2, 
  timeVar = "times_yr", 
  method = "piecewise-PH-aGH",   #nspeciﬁedbaseline risk function. This assumption turns out to b e equivalent toassuming that h0(·) is discrete with point-masses at the unique eventtime
  verbose = T,
  control=c(iter.EM=3000, lng.in.kn=10)
)#no...




for(i in 1:length(body_cols_names)){
  set.seed(55)
  joint_out1[[i]]<-jointModel(
    lme_out[[i]], 
    cox.fit, 
    timeVar = "times_yr", 
    method = "piecewise-PH-aGH",   #nspeciﬁedbaseline risk function. This assumption turns out to b e equivalent toassuming that h0(·) is discrete with point-masses at the unique eventtime
    verbose = T,
    control=c(iter.EM=3000, lng.in.kn=12)
  )
} 





# analysis ----------------------------------------------------------------



set.seed(1284)
library(lcmm)
fit1<-Jointlcmm(IMAT_cm2~times*ageAtMetDiag,
                random=~times,
                subject="PatientID_fct",
                #mixture=~baseline_to_scan,
                survival=Surv(new_OS_time, new_OS_ind) ~ smoking+ageAtMetDiag,
                hazard="3-quant-splines",
                #hazardtype="PH",
                #ng=2,
                data=dat) 


