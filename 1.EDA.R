#2024 Nov 04 (Mon)
#Explore the data

rm(list=ls())
library(data.table)
library(table1)
# data --------------------------------------------------------------------
dat<-readRDS("clean_dat.rds")

dat[, bmi_cateogrical_new:=as.character(bmi_cateogrical)]
dat[bmi_cateogrical_new %in% c("Underweight (< 18kg/m2)", "Normal weight (>= 18kg/m2 & < 25kg/m^2)"), bmi_cateogrical_new:="Underweight/Normal weight (<25kg/m2)"]
dat[bmi_cateogrical_new %in% c("Obesity class I (>= 30kg/m2 & < 35kg/m^2)", "Obesity class II (>= 35kg/m2 & < 40kg/m^2)"), bmi_cateogrical_new:="Obesity class I/II (>= 30kg/m2 & < 40kg/m2)"]
dat[,
    bmi_cateogrical_new:=
      factor(bmi_cateogrical_new,
             levels = c("Underweight/Normal weight (<25kg/m2)","Overweight (>= 25kg/m2 & < 30kg/m^2)","Obesity class I/II (>= 30kg/m2 & < 40kg/m2)"),
             labels = c("Underweight/Normal weight (<25kg/m2)","Overweight (>= 25kg/m2 & < 30kg/m^2)","Obesity class I/II (>= 30kg/m2 & < 40kg/m2)"))]

dat[,bmi_cateogrical:=NULL]
saveRDS(dat,"clean_dat_new_BMI.rds")

# Table1 ------------------------------------------------------------------

# function that generates Median [Min, Max] in Table 1 for continuous variables 
my.render.cont <- function(x) {
  with(stats.default(x), c("", "Median [Min, Max]" =
                             sprintf("%0.1f [%0.1f, %0.1f]", MEDIAN, MIN, MAX)))
}

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

body_cols<-colnames(dat)[!colnames(dat) %in% c("mrn","PatientID","PatientID_fct","is_baseline","new_study_date","StudyDate","times","times_mth","times_yr","earliest_studydate","metastaticDate",covar_cols,outcomes)]

body_cols_names<-
  sapply(body_cols,function(x){
    
    pieces_x<-strsplit(x,";")[[1]]
    
    if("cross_sectional_area_cm2" %in% pieces_x){
      paste(pieces_x[2],"(cm2)",sep=" ")
      
    }else if("volume_cm3" %in% pieces_x){
      paste(pieces_x[2],"(cm3)",sep=" ")
      
    }else if("HU_mean" %in% pieces_x){
      paste(pieces_x[2],"(mean HU)",sep=" ")
      
    }else if(stringr::str_detect(string = x,pattern = "whole")){
      pieces_x2<-strsplit(x,"_")[[1]]
      
      paste(pieces_x2[1],"(whole body cm3)",sep=" ")
    }
    
  })
body_cols_names<-do.call(c,body_cols_names)


#Just use the baseline data
dat[,is_baseline:=0]
dat[StudyDate==earliest_studydate,is_baseline:=1]
tmp<-dat[is_baseline==1,]
tmp<-tmp[,.SD,.SDcols = c("PatientID",names(covar_cols),names(body_cols_names))] |> unique()
setnames(tmp,
         old=c(names(covar_cols),names(body_cols_names)),
         new=c(covar_cols,body_cols_names))
all_cols_for_table1<-c(covar_cols,body_cols_names)
my_formula<-paste("`",all_cols_for_table1,"`", sep="")
my_formula<-paste(my_formula,collapse = "+")
my_formula<-paste("~",my_formula)
my_formula<-as.formula(my_formula)
t1_out<-table1(my_formula, render.continuous=my.render.cont, data = tmp)
t1_out

# check 1L therapy date ---------------------------------------------------

#2024 Dec 11 (Wed)
#Adriana and I chatted.
#We want to use 1L thearpy base baseline
#If everyone started within -60/+30 days of diagnosis date
#So check the date.

new_trt3<-readRDS("../../Body composition/code/trt_type_data.rds")
dat<-readRDS("clean_dat.rds")


#Restrict to 49 patients of interest
new_trt3<-new_trt3[mrn %in% dat$mrn,]

#metastaticDate the same?
tmp<-merge(new_trt3[,.(mrn,metastaticDate)],dat[,.(mrn,metastaticDate)],by="mrn")
tmp[,metastaticDate.x==metastaticDate.y] |> table() #TRUE. Yes.

#Now check the therapy date.
library(lubridate)
new_trt3[,dig_to_1L:=as.numeric(interval(metastaticDate,initial_trt_date),"days")] 
hist(new_trt3$dig_to_1L,breaks=20, main="# days since diag to 1L")
new_trt3[dig_to_1L>30,.N] #28 people started 30 days after diagnosis
#so we can't use 1L therapy type as a baseline covariate.

#check if there are scans post 1L initiation
tmp<-merge(dat[,.(mrn,times)],new_trt3[,.(mrn,dig_to_1L)])
tmp[,max_time:=max(times),by=mrn]
tmp<-tmp[times==max_time,]
tmp[dig_to_1L>=max_time,] #3 patients who received 1L have scans only before or on the day of 1L therapy.

#Anyways, we can't use 1L therapy as baseline covariate.

#42 people have scans post 1L initiation.
#Do we have scans 1 year post 1L therapy?
tmp[,yr_after_1L:=dig_to_1L+365]

#Replace yr_after_1L with death date. B/c if they died, then we definitely won't have scans 
tmp2<-merge(tmp[,.(mrn,max_time,yr_after_1L)],unique(dat[,.(mrn,new_OS_time)]),by="mrn", all.x=T, all.y = F)
tmp2[yr_after_1L>new_OS_time,yr_after_1L:=new_OS_time] 
tmp2[,.(max_time,yr_after_1L)] |> View()

#Okay the scans were definitely collected 1yr post baselien.
#Adriana told methat via chat as I was working on this.

# plot scan over time -----------------------------------------------------

library(ggplot2)
library(gridExtra)

all_plots<-
  lapply(1:length(body_cols_names),function(i){
    
    name_x<-names(body_cols_names[i])
    pretty_x<-body_cols_names[i]
    ggplot(data=dat)+
      geom_point(aes(x=times_yr, y=get(name_x),group = PatientID)) +
      geom_line(aes(x=times_yr, y=get(name_x),group = PatientID), alpha=0.2)+
      stat_smooth(
        data=dat,
        mapping=aes(x=times_yr, y=get(name_x)),
        geom="smooth",
        method="lm",
        formula=y~ splines::ns(x,df=3), 
        col="red")+
      ylab(pretty_x)+
      xlab("Time (year)")
    
  })


#Plot for all body parts
tmp<-arrangeGrob(grobs = all_plots, nrow=4, ncol=4)
pdf("plot1.pdf", width=22,height = 20)
grid.arrange(tmp,ncol=1,nrow=1)
dev.off()

#It's a bit hard to look at the slope b/c the intercepts differ so much!
#So when I fit repeated model, I should put random intercept model. 
#all plots look linear except for SKM(cm3),IMAT (whole body cm3),SAT (whole body cm3),SKM (whole body cm3)

all_plots[[8]]
all_plots[[13]]
all_plots[[14]]
all_plots[[15]]

#Why is slope increasing??


# debug increasing slope --------------------------------------------------


raw_dat<-readRDS("dat_new_death_times.rds") #this is non-imputed data
body_cols_names<-readRDS("body_cols_names.rds")
#We will only use the individual body part
#The whole scan that adds other parts of body are not plotted
#b/c it's a composite outcome. Any patient can miss any of those compoenets, 
#Which adds too much noise. 
names_body_cols_names<-stringr::str_subset(names(body_cols_names),";")
body_cols_names<-body_cols_names[names_body_cols_names]

#Just check the people who are missing any one of them
missing_index<-
  lapply(names(body_cols_names),function(x){
    raw_dat[,which(is.na(get(x)))]
  })
# missing_index
#This looks like the same scan is missing the value
#for all L3 section.


all_plots<-
  lapply(1:length(body_cols_names),function(i){
    
    pretty_name<-body_cols_names[i]
    x<-names(body_cols_names)[i]
    
    #To fit lmer, we have to take non-missing values
    raw_dat_tmp<-na.omit(raw_dat, cols=x)
    
    #Save the number of patients <- add to plot title 
    n<-length(unique(raw_dat_tmp$PatientID))
    impute_n<-length(unique(dat$PatientID))
    
    
    # no_impute_plot<-
    # ggplot(data=raw_dat_tmp)+
    #   geom_point(aes(x=times, y=get(x),group = PatientID)) +
    #   geom_line(aes(x=times, y=get(x),group = PatientID), alpha=0.2)+
    #   stat_smooth(
    #     data=dat,
    #     mapping=aes(x=times, y=get(x)),
    #     geom="smooth",
    #     method="lm",
    #     formula=y~ splines::ns(x,df=3), 
    #     col="red")+
    #   ylab(pretty_name)+
    #   xlab("Time (days)")+
    #   labs(title="raw data",
    #           subtitle = paste("N=", n,", sample size=",nrow(raw_dat_tmp), sep=""))
    # 
    # 
    # impute_plot<-
    # ggplot(data=dat)+
    #   geom_point(aes(x=times, y=get(x),group = PatientID)) +
    #   geom_line(aes(x=times, y=get(x),group = PatientID), alpha=0.2)+
    #   stat_smooth(
    #     data=dat,
    #     mapping=aes(x=times, y=get(x)),
    #     geom="smooth",
    #     method="lm",
    #     formula=y~ splines::ns(x,df=3), 
    #     col="red")+
    #   ylab(NULL)+
    #   xlab("Time (days)")+
    #   labs(title="imputed data",
    #           subtitle=paste("N =", impute_n,", sample size=",nrow(dat), sep=""))
    
    
    # 
    # this.plt<-ggarrange(no_impute_plot,impute_plot,
    #                     ncol=2, nrow=1)
    
    
    
    total_dat<-data.table(
      dat_type=rep(c("Raw data","Imputed data"),times=c(nrow(raw_dat_tmp), nrow(dat))),
      value=c(raw_dat_tmp[,get(x)], dat[,get(pretty_name)]),
      PatientID=c(raw_dat_tmp$PatientID, dat$PatientID),
      `Time (year)`=c(raw_dat_tmp$times_yr, dat$times_yr)
    )
    
    total_dat[dat_type=="Raw data",dat_type:=paste("Raw data","\n","N=", n,", sample size=",nrow(raw_dat_tmp), sep="")]
    total_dat[dat_type=="Imputed data",dat_type:=paste("Imputed data","\n","N=", impute_n,", sample size=",nrow(dat), sep="")]
    
    p<-ggplot(total_dat)+
      geom_point(aes(x=`Time (days)`, y=value,group = PatientID)) +
      geom_line(aes(x=`Time (days)`, y=value,group = PatientID) ,alpha=0.2)+
      stat_smooth(
        data=total_dat,
        mapping=aes(x=`Time (days)`, y=value),
        geom="smooth",
        method="lm",
        formula=y~ splines::ns(x,df=3),
        col="red")+
      facet_wrap(~dat_type)+
      theme(legend.position = "none")+
      ggtitle(pretty_name)+
      # ggtitle(paste(clean_name,"N =", n))+
      ylab(sub("\\).*", "", sub(".*\\(", "", pretty_name)) ) #Extract within the parenthesis
    return(p)
    
  })



# #Change to clean names
# setnames(raw_dat,
#          old=names(body_cols_names),
#          new=body_cols_names)
## Don't this above! the nlme::lme function cannot take variable names with space even if I wrap it with backquote. 
## nlme::lme cannot seem to take ; as well, so let's get rid of all special characters
body_cols_names_lme<-stringr::str_replace_all(names(body_cols_names),";","__")
names(body_cols_names_lme)<-names(body_cols_names)
setnames(raw_dat,
         old=names(body_cols_names_lme),
         new=body_cols_names_lme)

all_plots<-
lapply(body_cols_names_lme,function(x){
  
  #Get the name we want to replace with later for plotting
  clean_name<-body_cols_names[body_cols_names_lme==x]

  #To fit lmer, we have to take non-missing values
  raw_dat_tmp<-na.omit(raw_dat, cols=x)

  #Save the number of patients <- add to plot title 
  n<-length(unique(raw_dat_tmp$PatientID))
            
            
  fixed_formula<-paste(x,"~ageAtMetDiag+gender+smoking+splines::bs(times_yr, degree=1,df = 3)", sep="")
  fixed_formula<-as.formula(fixed_formula)
  fit<-lme(fixed_formula,random= ~ times_yr | PatientID_fct,data=raw_dat_tmp)
  
  
  
  
  # #How do we get predicted value?
  # fit$fitted #<-this is already available.
  # 
  # #To handcode the fitted value, we do the following:
  # newdat<-
  #   data.table(raw_dat_tmp[,.(PatientID_fct,ageAtMetDiag,gender,smoking,times_yr)],
  #              splines::bs(raw_dat_tmp$times_yr, degree=1,df = 3, intercept = F))
  # setnames(newdat,
  #          old=c("1","2","3"),
  #          new=paste("splines::bs(times_yr, degree = 1, df = 3)",1:3,sep=""))
  # 
  # 
  # fixed_coef<-t(t(fit$coefficients$fixed))
  # 
  # newdat2<-data.table(`(Intercept)`=1,newdat[,.(ageAtMetDiag,gender,smoking,times_yr,`splines::bs(times_yr, degree = 1, df = 3)1`,`splines::bs(times_yr, degree = 1, df = 3)2`,`splines::bs(times_yr, degree = 1, df = 3)3`)])
  # newdat2[,genderM:=ifelse(gender=="M",1,0)]
  # newdat2[,`smokingNo history of smoking`:=ifelse(smoking=="No history of smokin",1,0)]
  # newdat2[,c("gender","smoking"):=.(NULL)]
  # newdat2<-newdat2[,.(`(Intercept)`,ageAtMetDiag,genderM,`smokingNo history of smoking`,`splines::bs(times_yr, degree = 1, df = 3)1`,`splines::bs(times_yr, degree = 1, df = 3)2`,`splines::bs(times_yr, degree = 1, df = 3)3`)]
  # newdat2<-as.matrix(newdat2)
  # newdat2  %*%  fixed_coef |> head()  #this is the fixed effect 
  # fit$fitted |> head() #same as the above in the 'fixed' column
  # 
  # #What about if we include random effect?
  # fit_random<-fit$coefficients$random$PatientID_fct
  # order_id<-rownames(fit_random)
  # fit_random<-data.table(fit_random)
  # fit_random[,PatientID_fct:=order_id]
  # fit_random<-merge(raw_dat_tmp[,.(PatientID_fct)],fit_random,by="PatientID_fct",all.x = T)
  # rand_effect<-
  #   sapply(1:nrow(fit_random),function(i){
  #     fit_random[i,`(Intercept)`] + fit_random[i,times_yr] * raw_dat_tmp$times_yr[i]
  #   })
  # rand_effect<-t(t(rand_effect))
  # (newdat2  %*%  fixed_coef  + rand_effect) |> head() #this is the estimated effect accounting for the random estimate 
  # fit$fitted |> head() #same as the above in the 'PatientID_fct' column
  # 
  
  #Okay we will use the individual fitted value 
  # fit_dat<-data.table(fit$data[,.SD,.SDcols=c("PatientID_fct","times_yr")],
  #                     observed=raw_dat_tmp[,get(x)],
  #                     # fitted=fit$fitted[,"fixed"] ,
  #                     fitted_individual=fit$fitted[,"PatientID_fct"])
  # fit_dat<-melt(fit_dat,
  #      id.vars=c("PatientID_fct","times_yr"), 
  #      measure.vars = c("observed","fitted_individual"))
  # # setnames(fit_dat,old="value",new=clean_name)
  # fit_dat[variable=="observed",variable:="Observed"]
  # fit_dat[variable=="fitted_individual",variable:="Regression Fit"]
  # 
  # p<-ggplot(fit_dat,aes(x=times_yr, y=value,group = PatientID_fct, col=PatientID_fct))+
  #   geom_point() +
  #   geom_line(alpha=0.2)+
  #   facet_wrap(~variable)+
  #   theme(legend.position = "none")+
  #   ggtitle(clean_name)+
  #   # ggtitle(paste(clean_name,"N =", n))+
  #   ylab(clean_name)+
  #   xlab("Times (year)")
  # return(p)  
  fit_dat<-data.table(fit$data[,.SD,.SDcols=c("PatientID_fct","times_yr",x)],
                      # fitted=fit$fitted[,"fixed"] ,
                      fitted_individual=fit$fitted[,"PatientID_fct"])


  #Plot before & after fitting linear mixed effects model
  library(gridExtra)
  p1<-  ggplot(raw_dat_tmp)+
    geom_point(aes(x=times_yr, y=.data[[x]],group = PatientID, col=PatientID)) +
    geom_line(aes(x=times_yr, y=.data[[x]],group = PatientID, col=PatientID), alpha=0.2)+
    ylab(clean_name)+
    theme(legend.position = "none")+
    ggtitle("Connecting observations")+
    xlab("Times (year)")
  
  p2<-  ggplot(fit_dat)+
    geom_point(aes(x=times_yr, y=.data[[x]],group = PatientID_fct, col=PatientID_fct)) +
    geom_line(aes(x=times_yr, y=fitted_individual,group = PatientID_fct, col=PatientID_fct), alpha=0.5)+
    theme(legend.position = "none")+
    ylab(clean_name)+
    ggtitle("Regression lines")+
    xlab("Times (year)")
  
  

  this.plt<-ggarrange(p1,p2,
                      ncol=2, nrow=1)

  #add the overall plot label
  annotate_figure(
    this.plt,
    top = text_grob(label=clean_name,
                    color = "black", face = "bold", size = 20))
  

  
})


all_plots[[2]]
library(ggpubr)
pdf("LMER_plots_tmp.pdf", onefile = T, width = 10, height = 4)
all_plots

# for( i in 1:length(all_plots)){
# 
#   this.plt<-ggarrange(all_plots[[i]]$this.plot[[1]],all_plots[[i]]$this.plot[[2]],
#                       ncol=2, nrow=1)
#   # grid.arrange(
#   #   arrangeGrob(grobs=all_plots[[i]]$this.plot, ncol=2,widths = c(2,1)),
#   #   nrow=1
#   # )
# 
# 
#   annotate_figure(
#     this.plt, 
#     top = text_grob(label=paste("# of unique patients =", all_plots[[i]]$n), 
#                     color = "black", face = "bold", size = 14))
#   
#   grid::grid.newpage()
# }
dev.off()

#Look at these plots. They look like they are 
# KM plot -----------------------------------------------------------------


rm(list=ls())
dat<-readRDS("clean_dat_new_BMI.rds")

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

body_cols<-colnames(dat)[!colnames(dat) %in% c("mrn","PatientID","PatientID_fct","is_baseline","new_study_date","StudyDate","times","times_mth","times_yr","earliest_studydate","metastaticDate",covar_cols,outcomes)]

body_cols_names<-
  sapply(body_cols,function(x){
    
    pieces_x<-strsplit(x,";")[[1]]
    
    if("cross_sectional_area_cm2" %in% pieces_x){
      paste(pieces_x[2],"(cm2)",sep=" ")
      
    }else if("volume_cm3" %in% pieces_x){
      paste(pieces_x[2],"(cm3)",sep=" ")
      
    }else if("HU_mean" %in% pieces_x){
      paste(pieces_x[2],"(mean HU)",sep=" ")
      
    }else if(stringr::str_detect(string = x,pattern = "whole")){
      pieces_x2<-strsplit(x,"_")[[1]]
      
      paste(pieces_x2[1],"(whole body cm3)",sep=" ")
    }
    
  })
body_cols_names<-do.call(c,body_cols_names)

setnames(dat,
         old=names(body_cols_names),
         new=body_cols_names)



#Look at median survival


dat[,is_baseline:=0]
dat[StudyDate==earliest_studydate,is_baseline:=1]
baseline_dat<-dat[is_baseline==1,]
fit1 <- survfit(Surv(new_OS_time_yr, new_OS_ind) ~ 1, data = baseline_dat)
summary_fit1<-summary(fit1)$table
source("survfit_median_CI.R")
survfit_median_CI(fit1)



#Kaplan meier curve
source("return_km_plot.R")

km_plot<-return_km_plot(data = baseline_dat, time="new_OS_time_yr",xlab_name = "Times (year)",status ="new_OS_ind" )

#Remove legend
# Assuming your plot is saved as 'x'
km_plot$plot <- km_plot$plot + theme(legend.position = "none")
pdf("km_plot.pdf", width=5,height=4, onefile = FALSE)
km_plot
dev.off()

