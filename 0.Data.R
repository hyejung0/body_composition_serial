#2024 Nov 04 (Mon)
#Generate data we need to use

rm(list=ls())
library(data.table)

# data --------------------------------------------------------------------


#Read in the patient baseline data
dat1_revised<-readRDS("../../Body composition/code/dat1_revised.rds") #demographic data from body composition analysis (not serialized)
colnames(dat1_revised)

trt_type_data<-readRDS("../../Body composition/code/trt_type_data.rds")
BMI<-readRDS("../../Body composition/code/BMI.rds")
tumor_histology<-readRDS("../../Body composition/code/tumor_histology.rds")

old_dat<-
  merge(dat1_revised[,.(mrn, `Subject-ID`,metastaticDate,ageAtMetDiag,gender,smoking,raceEthnicity,Gene,new_OS_time,new_OS_ind)],
        trt_type_data[,.(mrn,trt_type)],
        by="mrn",
        all.x=TRUE)
old_dat[is.na(trt_type),trt_type:="Not treated"]
old_dat<-merge(old_dat,BMI[,.(mrn,bmi,bmi_cateogrical)])
old_dat<-merge(old_dat,tumor_histology[,.(mrn,tumor_hist_breakdwon)],by="mrn")
length(unique(old_dat$`Subject-ID`))
colnames(old_dat)
#remove individual data set
rm(list=c("dat1_revised","trt_type_data","BMI","tumor_histology"))


#Below here are the serialized scan data
#comprehensive file
# comp_scan<-fread("../../../../LCC Pilot Grant Project 2024/Data/Comprehensive Analysis/Publish_by_scan_2024-02-23_16.57.54.389.csv",nThread=4) 
comp_scan<-fread("../../../../LCC Pilot Grant Project 2024/Data/Comprehensive Analysis/Updated Results 4-22-2024/Publish_by_scan_2024-11-08_06.49.32.807_JY.csv",nThread=4) 
comp_id<-comp_scan$PatientID|> unique()

#The data is too big. We will subset the data so we don't have such as big data. 
#Adriana provided columns that we need to use. 
dat_dict<-fread("../Documents/variables_dictionary.csv",na.strings = c("", NA))
library(zoo)
dat_dict<-na.locf(dat_dict)
stringr::str_subset(colnames(comp_scan),";", negate = F) |> length() #16368 body composition names
select_comp_scan<-c(stringr::str_subset(colnames(comp_scan),";", negate = T), #Column names that are not body parts scans.
                    unique(dat_dict$`Variable Name`))
comp_scan<-comp_scan[,.SD,.SDcols=select_comp_scan]
#convert to numeric
scan_cols<- unique(dat_dict$`Variable Name`)
comp_scan[,(scan_cols):=lapply(.SD,as.numeric),.SDcols=scan_cols]

# Modify data -------------------------------------------------------------



length(comp_id)
#Exclude one person who was changed to non-metastatic by Wally on 2024 Sept 25. PatientID=="108345_00056"
comp_scan[PatientID=="108345_00056",]
#Why isn't it here?

#Read in the file used to identify valid patients with paired chest and abdomen scan.
Hyejung_new<-fread("../../../../LCC Pilot Grant Project 2024/Data/paired_scan_dates Final.csv",nThread=4) 
Hyejung_new[,.(mrn,`Subject-ID`)]
Hyejung_new[`Subject-ID`=="108345_00056",.(mrn,`Subject-ID`)]
#I guess it was excluded before...?
#check that all serial scan patients are in my old dataset.
comp_id %in% old_dat$`Subject-ID` |> table() #One isn't in my previous dataset
comp_id[!(comp_id %in% old_dat$`Subject-ID`)] #108345_00076
#What MRN is that?
Hyejung_new[`Subject-ID`=="108345_00076", mrn] #20895589
#This is the patient who Wally changed metastatic date to 2015-05-20
#Note that this person got excluded in revised original analysis (non-serialized analysis)
#because number of days from metastatic diagnosis date to baseline scan date fell outside -60 to +30 days.
#so exclude the person. 
#Okay so remove this person from the serialized dataset.
comp_scan<-comp_scan[PatientID!="108345_00076",]




#Merge the new serial scan data to old covariate data
#first of all, check that all PatientID is in the old data
table(unique(comp_scan$PatientID) %in% old_dat$`Subject-ID`)
#TRUE 
#  50 
#Yes. Merge.
rm(Hyejung_new) #I don't need this data anymore.

comp_scan<-
  merge(comp_scan,
        old_dat,
        by.x="PatientID",
        by.y="Subject-ID",
        all.x=T,
        allow.cartesian = T)
rm(old_dat)
setnames(
  comp_scan,
  old="tumor_hist_breakdwon",
  new="tumor_hist"
)


# quality check -----------------------------------------------------------

#Baseline scan dates different?
comp_scan[,.(PatientID,orig_path,StudyDate,StudyID,BodyPartExamined)]
unique(comp_scan$BodyPartExamined)
comp_scan[PatientID=="108345_00009" & StudyDate==20160608,]


#change the StudyDate to Date
StudyDate_IDate<-comp_scan$StudyDate |> as.character() |> as.IDate("%Y%m%d")
comp_scan[,StudyDate:=StudyDate_IDate]

#2024 Nov 08 (Fri)
#Talked with Jeff to figure out how to get the data.

#I should have max 2 observation per patient per study date. Check. 
comp_scan[,.N,by=.(PatientID,StudyDate)][,range(N)] # 1 2. YES. I'm right.
#Make sure that 2 copies don't include CHESTABDOMEN.
comp_scan[,count:=.N,by=.(PatientID,StudyDate)]
comp_scan[count==2,unique(BodyPartExamined)] #doesn't include CHESTABDOMEN. Good!
#IF there are two observations, ever, then it's only ABDOMEN and CHEST combination
#What about study that involves just one observation? Is it just CHESTABDOMEN?
comp_scan[count==1,unique(BodyPartExamined)]  #"CHESTABDOMEN" "ABDOMEN"      "CHEST". There are just single abdomen and chest. 

comp_scan[count==1 & BodyPartExamined=="ABDOMEN",.N]  #8 people have ABDOMEN but not the paired CHEST
comp_scan[count==1 & BodyPartExamined=="CHEST",.N]  #6 people have ABDOMEN but not the paired CHEST


#get the column names of interest, the location of scan.
spinal_cord_location<-substr(scan_cols, start = 1, stop = 2) #this returns the spinal cord location. Except.a feew have Sa. What are those?
spinal_cord_location |> unique()
scan_cols[spinal_cord_location=="Sa"] #Sacrummid
#I googled Sacrummid. It was Sacrum mid. It's below L. Should be part of ABDOMEN scan. 
L_and_SA_cols<-scan_cols[stringr::str_detect(spinal_cord_location,"T",negate = T)]
T_cols<-scan_cols[stringr::str_detect(spinal_cord_location,"T",negate = F)]



get_scan_value<-function(dat){
  #Per patient per studydate information
  
  #We should have either ABDOMEN and CHEST or CHESTABDOMEN.
  if(length(dat$BodyPartExamined)>1){#This is the case where I have two copies of scan on the particular studydate. This is only the case if we have ABDOMEN and CHEST data
    #If we have both ABDOMEN and CHEST copies.

    #If any L_and_SA_cols are missing in the ABDOMEN, fill that with CHEST
    temp<-rbind(
      dat[BodyPartExamined=="CHEST",..L_and_SA_cols],
      dat[BodyPartExamined=="ABDOMEN",..L_and_SA_cols]
    )
    complete_L_and_SA_cols<-na.locf(temp,na.rm = FALSE) #na.rm=FALSE will keep the columns even if both rows are NA. 
    complete_L_and_SA_cols<-complete_L_and_SA_cols[2,] #2nd row is the BodyPartExamined=="ABDOMEN" row. We should keep that
    
    #Likewise, do that for CHEST variables.
    temp<-rbind(
      dat[BodyPartExamined=="ABDOMEN",..T_cols],
      dat[BodyPartExamined=="CHEST",..T_cols]
    )
    complete_T_cols<-na.locf(temp,na.rm = FALSE)
    complete_T_cols<-complete_T_cols[2,] #2nd row is the BodyPartExamined=="CHEST" row. We should keep that
    

    #Then we should just return one row. 
    #the complete_L_and_SA_cols columns should be returned from row BodyPartExamined=="ABDOMEN",
    #and  columns for T_cols columns should should be returned from row BodyPartExamined=="CHEST",
    #And then also append most superior_vert_in_scan and inferior_vert_in_scan
    #the most superior should come from Chest scan b/c it's located above ABdomen
    #the most inferior should come from abdomen scan b/c it's located below chest
    
    output<-cbind(complete_L_and_SA_cols,complete_T_cols)
    output<-data.table(output)
    output[,most_superior_vert_in_scan:=dat[BodyPartExamined=="CHEST",most_superior_vert_in_scan]]
    output[,most_inferior_vert_in_scan:=dat[BodyPartExamined=="ABDOMEN",most_superior_vert_in_scan]]

    #Also indicate bodyparts examined
    output[,BodyPartExamined:="ABDOMEN/CHEST"]
  }else{
    #If there's only one line available, just return the whole data
    dat[,.SD,.SDcols = c(scan_cols,"most_superior_vert_in_scan","most_inferior_vert_in_scan","BodyPartExamined")]
  }
}

comp_scan2<-comp_scan[,get_scan_value(.SD), by=.(PatientID,StudyDate)]
dim(comp_scan2)
#check for missing entries
comp_scan2[,sapply(.SD,function(x)sum(is.na(x)))]#Lots of missing...
#now merge all other columns
cols_not_in_comp_scan2<-colnames(comp_scan)[!colnames(comp_scan) %in%colnames(comp_scan2)]
#remove the columns not needed
cols_not_in_comp_scan2
cols_not_in_comp_scan2<-cols_not_in_comp_scan2[!cols_not_in_comp_scan2 %in% c("scan_folder","orig_path","verts_in_scan","verts_with_mid_annot_in_scan","segmented_tissues","voxel_size_WxLxH","Modality","SeriesDescription","SeriesDate","SeriesNumber","count")]
#We need PatientID and StudyDate. So include them
cols_not_in_comp_scan2<-c(cols_not_in_comp_scan2,"PatientID","StudyDate")
tmp<-comp_scan[,.SD,.SDcols = cols_not_in_comp_scan2] |> unique()
dim(tmp)
dim(comp_scan2)
#same number of rows!
comp_scan3<-merge(tmp,comp_scan2,by=c("PatientID","StudyDate"),all.y=T)
comp_scan4<-merge(tmp,comp_scan2,by=c("PatientID","StudyDate"),all.y=T,all.x = T,allow.cartesian = T)
#They should have same patientID and studyID. Check
in_4<-comp_scan4[,paste(PatientID,StudyDate,sep ="__")]
in_3<-comp_scan3[,paste(PatientID,StudyDate,sep ="__")]
identical(in_4,in_3) #YES! there's no one who got missed out when creating comp_scan2
rm(comp_scan4)
#Use comp_scan3


#check if the earliest scan date is the same results as the baseline scan used for non-serialized analysis.
comp_scan3[,earliest_studydate:=min(StudyDate),by=PatientID]
length(unique(comp_scan3$PatientID)) #50 people
comp_earliest<-comp_scan3[StudyDate==earliest_studydate,]
nrow(comp_earliest) #50



#Grab the baseline scan date by extracting earliest study date per patient
comp_baseline<-comp_scan3[,min(StudyDate),by=PatientID]
setnames(comp_baseline,"V1","New_baseline_scan_date")
#Compare it to my baseline scan date
dat1_revised<-readRDS("../../Body composition/code/dat1_revised.rds") #demographic data from body composition analysis (not serialized)
dat1_revised<-dat1_revised[`Subject-ID` %in% comp_baseline$PatientID,]
comp_baseline<-merge(dat1_revised[,.(`Subject-ID`,metastaticDate,`Scan Date`)],
                     comp_baseline,
                     by.x="Subject-ID",
                     by.y="PatientID")
str(comp_baseline)
comp_baseline[,`Scan Date`:=as.IDate(`Scan Date`)]
comp_baseline[,metastaticDate:=as.IDate(metastaticDate)]
comp_baseline[,agree:=(`Scan Date`==New_baseline_scan_date)]
comp_baseline[agree==FALSE,]
#Key: <Subject-ID>
#      Subject-ID metastaticDate  Scan Date New_baseline_scan_date  agree
#          <char>         <IDat>     <IDat>                 <IDat> <lgcl>
# 1: 108345_00011     2017-03-09 2017-03-17             2017-03-08  FALSE
# 2: 108345_00012     2017-03-31 2017-03-20             2017-03-30  FALSE
# 3: 108345_00036     2016-03-24 2016-03-20             2016-03-30  FALSE
# 4: 108345_00038     2016-09-23 2016-10-05             2016-09-22  FALSE
# 5: 108345_00042     2016-03-24 2016-03-28             2016-03-23  FALSE
# 6: 108345_00062     2016-04-29 2016-05-23             2016-05-26  FALSE
# 7: 108345_00068     2016-10-12 2016-09-26             2017-03-09  FALSE
# 8: 108345_00078     2016-04-13 2016-04-01             2016-04-11  FALSE
# 9: 108345_00081     2016-12-15 2017-01-10             2016-12-16  FALSE
#10: 108345_00084     2016-11-04 2016-11-30             2016-11-03  FALSE

#The scan date from my revised version from non-serialized data anlysis (Scan Date)
#is different from the new serialized scan date for 10/50 people. 
#But everyone's baseline date from this serialized data (New_baseline_scan_date) 
#is closer to the metastatic date compared to the baseline date from the non-serialized data (Scan Date)
#except for 108345_00068
#Exclude the 108345_00068, b/c the earliest scan available is sitting outside the -60/+30 days window.
comp_scan3<-comp_scan3[PatientID!="108345_00068"]
comp_baseline[agree==FALSE & `Subject-ID` !="108345_00068",`New_baseline_scan_date`-metastaticDate] 
#[1] -1 -1  6 -1 -1 27 -2  1 -1
#all within -60/+30 days window!
#I thought about grabbing the baseline data from non-serialized data
#But Jeff said that that can bring different source of error so I am just going to use the new dataset
comp_scan3[,metastaticDate:=as.IDate(metastaticDate)]


# change first study date to metastatic date ------------------------------

#We assumed that thee first scan reflect that of metastatic date.
#we should change the first study date to the metastatic date 
#because if we keep the actual study date, some patients seem to have scan post death. For example:

comp_scan3[,times:=StudyDate-earliest_studydate]

#check that the patient's death date does not come before the scan date
comp_scan3[,.(PatientID,new_OS_ind,new_OS_time,times)]
tmptmp<-comp_scan3[,new_OS_time-times]
(tmptmp<0)|> table() #TRUE==7. Meaning, anywhere from one to seven people had death date prior to scan date. Which are these?
comp_scan3[tmptmp<0,.(PatientID,new_OS_ind,new_OS_time,times, metastaticDate,StudyDate)]
comp_scan3[tmptmp<0,.(mrn,new_OS_ind,new_OS_time,times, metastaticDate,StudyDate)]
#        mrn new_OS_ind new_OS_time times metastaticDate  StudyDate
#      <int>      <int>       <int> <int>         <IDat>     <IDat>
#1: 20826271          1         275   277     2016-06-21 2017-03-12
#2:  7044951          1         353   368     2016-03-25 2017-03-07
#3: 16666703          0         141   181     2017-02-01 2017-07-31
#4: 16666703          0         141   265     2017-02-01 2017-10-23
#5: 16666703          0         141   331     2017-02-01 2017-12-28
#6: 14695027          0         125   203     2017-02-24 2017-09-27
#7:  7868896          0          30    55     2016-06-21 2016-07-20
#here, you can see that for MRN==20826271, the person died at 275 post metastatic date
#But then had scans two days after. 
comp_scan3[mrn==20826271, .(StudyDate,metastaticDate)]
#but as you can see here, the study date that's closest to metastatic date is 
#2016-06-08, which is prior to metastaticDate. 
#The scan that showed 277 days was 2017-03-12, 
as.IDate("2017-03-12") - as.IDate("2016-06-08") #277
#If we changed the earliest study date to metastatic date (2016-06-21),
as.IDate("2017-03-12") - as.IDate("2016-06-21")  #264
#It's now shorter length than time to death. 

#We assumed that the first scan is reflective of the body component at baseline, so just change the studydate.
comp_scan3[,new_study_date:=StudyDate]
comp_scan3[StudyDate==earliest_studydate,new_study_date:=metastaticDate]
# dat[,.(PatientID,metastaticDate,StudyDate,new_study_date)] |> View() #good

#new_study_date column has modified first scan date, where the first scan is the metastaticDate
#now check again, if any scan is performed after death
library(lubridate)
comp_scan3[,times:=as.numeric(interval(metastaticDate,new_study_date),"days")]
comp_scan3[,times_mth:=as.numeric(interval(metastaticDate,new_study_date),"months")]
comp_scan3[,times_yr:=as.numeric(interval(metastaticDate,new_study_date),"years")]

tmptmp<-comp_scan3[,new_OS_time-times]
(tmptmp<0)|> table() #TRUE==5
comp_scan3[tmptmp<0,.(mrn,new_OS_ind,new_OS_time,times, metastaticDate,new_study_date)]
#        mrn new_OS_ind new_OS_time times metastaticDate new_study_date
#      <int>      <int>       <int> <int>         <IDat>         <IDat>
#1: 16666703          0         141   180     2017-02-01     2017-07-31
#2: 16666703          0         141   264     2017-02-01     2017-10-23
#3: 16666703          0         141   330     2017-02-01     2017-12-28
#4: 14695027          0         125   136     2017-02-24     2017-07-10
#5: 14695027          0         125   215     2017-02-24     2017-09-27

#two patients have censoring date that are shorter than the scan date.
#I will let Adriana know that I will change the censoring date to the last scan date for these two people.
comp_scan3[mrn==16666703,new_OS_time:=max(times)]
comp_scan3[mrn==14695027,new_OS_time:=max(times)]
tmptmp<-comp_scan3[,new_OS_time-times]
(tmptmp<0)|> table() #all FALSE. Meaning, no one's scan dates come after death dates.

saveRDS(comp_scan3,"dat.rds")


# Death dates -------------------------------------------------------------


#2024 Dec 4 (Wed)
#I requested Adriana to get the most updated death dates for people. 

dat<-readRDS("dat.rds")
death<-readxl::read_excel("../Data/RISR-4361.xlsx", sheet = "DeathDates")
death<-data.table(death)
str(death)
length(unique(death$mrn))
death[,.N,by=mrn][N>1,]
#        mrn     N
#      <num> <int>
#1: 20585729     2
#2: 16665630     2
#3: 21146661     2
#These three people have 2 copies of record. Why?
death[mrn %in% c(20585729,16665630,21146661)]
#        mrn  deathDate lastEdwActivity
#      <num>     <POSc>          <POSc>
#1: 20585729 2018-10-31            <NA>
#2: 20585729 2018-10-19            <NA>
#3: 16665630 2018-03-26            <NA>
#4: 16665630 2018-03-16            <NA>
#5: 21146661 2017-08-05            <NA>
#6: 21146661 2017-08-03            <NA>
#2024 Dec 05 (Thur)
#email back from Jeffrey 
#use the latest recorded date as date of death.
#Multiple dates comes from pulling different sources 
death[,deathDate:=as.IDate(deathDate)]
death[,lastEdwActivity:=as.IDate(lastEdwActivity)]

death<-death[!(mrn==20585729& deathDate==as.IDate("2018-10-19")),]
death<-death[!(mrn==16665630& deathDate==as.IDate("2018-03-16")),]
death<-death[!(mrn==21146661& deathDate==as.IDate("2017-08-03")),]


#Now that clean_dat has mrn, merge death and clean_dat datasets together.
death<-merge(death, unique(dat[,.(mrn,metastaticDate)]),by="mrn",all.x=T)

#generate death indicator
death[,new_OS_ind:=1]
death[is.na(deathDate),new_OS_ind:=0]

#collect the censoring column to death column to generate 
#one variable that has all times
death[is.na(deathDate),deathDate:=lastEdwActivity]
#Remove censoring time column
death[,lastEdwActivity:=NULL]
#Calculate time to event. 
library(lubridate)
death[,new_OS_time:=as.numeric(interval(metastaticDate,deathDate),"days")]
death[,new_OS_time_mth:=as.numeric(interval(metastaticDate,deathDate),"months")]
death[,new_OS_time_yr:=as.numeric(interval(metastaticDate,deathDate),"years")]

#Now save time and indicator as separate file
saveRDS(death,"updated_death.rds")


#merge the file to the clean_dat
dat2<-merge(dat,death,by="mrn",all.x = T)
stringr::str_subset(colnames(dat2),".x")
dat2[,c("new_OS_time.x","new_OS_time.y","new_OS_ind.x","new_OS_ind.y")] |> View()
#throw out old new_OS_ind and new_OS_time
#It looks like most of the new death/censor times are way longer...
#hmm....
#Well, I will take that.

dat2[,c("new_OS_time.x","new_OS_ind.x"):=NULL]
setnames(dat2,
        old=c("new_OS_time.y","new_OS_ind.y"),
        new=c("new_OS_time","new_OS_ind"))
#What about metastatic Date? It should be the exact same.
dat2[,metastaticDate.x==metastaticDate.y] |> table() #Yes. So remove one. 
dat2[,metastaticDate.x:=NULL]
setnames(dat2,"metastaticDate.y","metastaticDate")

saveRDS(dat2,"dat_new_death_times.rds")

# imputation --------------------------------------------------------------


rm(list=ls())
library(data.table)
dat<-readRDS("dat_new_death_times.rds")

dat_dict<-fread("../Documents/variables_dictionary.csv",na.strings = c("", NA))
library(zoo)
dat_dict<-na.locf(dat_dict)

body_cols<-dat_dict$`Variable Name` |> unique() #these are the column names I need that's the scan measures.
covar_cols<-c("ageAtMetDiag","gender","smoking","raceEthnicity","trt_type","bmi","bmi_cateogrical","tumor_hist","Gene" )
outcome_cols<-c("new_OS_time","new_OS_ind","new_OS_time_mth","new_OS_time_yr" )


#Check how many missing in each column
dat[,lapply(.SD,function(x)sum(is.na(x)))]

#Check how many missing for outcomes
dat[,lapply(.SD,function(x)sum(is.na(x))), .SDcols = outcome_cols] #none


#Check how many missing for covariances
dat[,lapply(.SD,function(x)sum(is.na(x))), .SDcols = covar_cols] #Lots missing raceEthnicity
tmp<-dat[,.(PatientID,raceEthnicity)] |> unique()
tmp[,table(raceEthnicity, useNA ='ifany')]  
#Hispanic or Latino     Non-Hisp Asian     Non-Hisp White              Other               <NA> 
#                 1                  1                 31                  1                 15 
#15 people are missing race/ethnicity. Majority is Non-Hisp White.
#Due to small sample size in all categories besides Non-Hisp White, we will combine them to "Other"
#People missing will be combined to Other b/c we only have 3 Other race. 
#We will make them "Other/Missing"
dat[raceEthnicity!="Non-Hisp White" | is.na(raceEthnicity),raceEthnicity:="Other/Missing"]
dat$raceEthnicity |> table(useNA = "ifany")




#Exclude columns that are not needed. That doesn't have any statistical values for imputation
str(dat)
#Take out StudyID and earliest_studydate



#change characters to factor before imputation
char_col_names<-colnames(dat)[dat[,sapply(.SD,is.character)]]
char_col_names<-char_col_names[char_col_names!="StudyID"]
#We decided to keep PatientID for imputation because we want to let the algorithm know that the set of observations come from the same patient.

str(dat[,.SD,.SDcols = char_col_names])
dat[,lapply(.SD,function(x)sum(is.na(x))), .SDcols = char_col_names]
#Now really chance to factors 
lapply(char_col_names,function(x){ #change to factor
  unique_x<-dat[,unique(get(x))]
  #remove NA from the unique entry
  unique_x<-unique_x[!is.na(unique_x)]
  dat[,(x):=factor(get(x),levels = unique_x,labels = unique_x)]
})
str(dat)


#For imputation, we will exclude mrn b/c it is exactly same meaning as PatientID


#We will perform imputation. 
library(missRanger)
test<-
missRanger(
  data=dat,
  formula=.~.-mrn-StudyID-earliest_studydate-times_yr-times_mth-new_OS_times_mth-new_OS_times_yr,
  pmm.k = 5L,
  maxiter = 10L, 
  seed = 111,
  verbose = 1,
  returnOOB = FALSE,
  keep_forests = TRUE,
  num.trees = 500
)
saveRDS(test,"impute_missRanger_dat.rds")



#Subset to columns we need
dat_dict$`Variable Category` |> unique()
#[1] "IMAT at L3"                                                   "SAT at L3"                                                   
#[3] "SKM at L3"                                                    "VAT at L3"                                                   
#[5] "IMAT- Whole abdomen + chest volume  (sum of these variables)" "SAT- Whole abdomen + chest volume  (sum of these variables)" 
#[7] "SKM- Whole abdomen + chest volume  (sum of these variables)"  "VAT- Whole abdomen + chest volume  (sum of these variables)" 

#Here, We first subset our data to include variables in "IMAT at L3", "SAT at L3","SKM at L3","VAT at L3"
outcome_cols<-c("new_OS_time","new_OS_ind", "new_OS_time_mth","new_OS_time_yr")

clean_dat<-test$data[
  ,.SD
  ,.SDcols = c("mrn",
               "PatientID", 
               "StudyDate",
               "new_study_date",
               "times",
               "times_mth",
               "times_yr",
               "earliest_studydate",
               "metastaticDate",
               outcome_cols,
               covar_cols,
               dat_dict[`Variable Category` %in% c("IMAT at L3","SAT at L3","SKM at L3","VAT at L3"),`Variable Name`])]
#now we create the composite variables 

#"IMAT- Whole abdomen + chest volume  (sum of these variables)"
IMAT_ab_chest_cols<-dat_dict[`Variable Category`=="IMAT- Whole abdomen + chest volume  (sum of these variables)", `Variable Name`]
IMAT_whole_cm3<-test$data[,.SD,.SDcols = IMAT_ab_chest_cols] |> rowSums()
clean_dat[,IMAT_whole_cm3:=IMAT_whole_cm3]

#"SAT- Whole abdomen + chest volume  (sum of these variables)" 
SAT_ab_chest_cols<-dat_dict[`Variable Category`=="SAT- Whole abdomen + chest volume  (sum of these variables)", `Variable Name`]
SAT_whole_cm3<-test$data[,.SD,.SDcols = SAT_ab_chest_cols] |> rowSums()
clean_dat[,SAT_whole_cm3:=SAT_whole_cm3]

#"SKM- Whole abdomen + chest volume  (sum of these variables)"
SKM_ab_chest_cols<-dat_dict[`Variable Category`=="SKM- Whole abdomen + chest volume  (sum of these variables)", `Variable Name`]
SKM_whole_cm3<-test$data[,.SD,.SDcols = SKM_ab_chest_cols] |> rowSums()
clean_dat[,SKM_whole_cm3:=SKM_whole_cm3]

#"VAT- Whole abdomen + chest volume  (sum of these variables)"
VAT_ab_chest_cols<-dat_dict[`Variable Category`=="VAT- Whole abdomen + chest volume  (sum of these variables)", `Variable Name`]
VAT_whole_cm3<-test$data[,.SD,.SDcols = VAT_ab_chest_cols] |> rowSums()
clean_dat[,VAT_whole_cm3:=VAT_whole_cm3]


saveRDS(clean_dat,"clean_dat.rds")



