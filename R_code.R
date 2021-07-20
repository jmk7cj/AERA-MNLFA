#--------------------------------------------------------------------------------------------#
# Author: Joseph Kush (jkush1@jhu.edu) 
# 
# Title: Utilizing Moderated Nonlinear Factor Analysis Models for Integrative Data Analysis
#
# Date: 7/19/2021
#
# Purpose: Master .R file to conduct integrative data analyses using moderated 
#          non-linear factor analysis
#          Step 0: Load packages, set working directory, import data, etc.
#          Step 1: Create Mplus input files for EFAs, estimate models, examine output
#          Step 2: Create Mplus input files for MNLFA model building, estimate models, 
#                  examine output
#          Step 3: Conduct LRT between each item-model and the baseline model
#          Step 4: Remove remaining non-significant parameters, estimate next-to-last and 
#                  final MNLFA model
#          Step 5: Merge estimated factor scores to be used in subsequent analyses
#--------------------------------------------------------------------------------------------#





#--------------------------------------------------------------------------------------------#
# Step 0: Load packages, set working directory, import data, etc.
#--------------------------------------------------------------------------------------------#
# Remove working environment, close any open connections
rm(list = ls()); closeAllConnections()

# Load necessary packages
library("parallel")
library("MplusAutomation")
library("MASS")

# Set working directory to folder
myfolder <- setwd("/Users/joekush/Desktop/myfolder")

# Import data, label variables
data <- read.csv("mnlfa.csv", header=T)
data <- subset(data, study_id != 0 & study_id != 2 & study_id != 3)
data$study_id <- factor(data$study_id, labels=c("LIFT_L", "PIRC1", "PIRC2", "SAFE"))
data$sex <- factor(data$sex, labels=c("female", "male"))
data$race <- factor(data$race, labels=c("white", "black"))

head(data)
#--------------------------------------------------------------------------------------------#





#--------------------------------------------------------------------------------------------#
# Step 1: Create Mplus input files for EFAs, estimate models, examine output
#--------------------------------------------------------------------------------------------#
# First, prepare datafile for Mplus
data_efa <- data
data_efa[is.na(data_efa)] <- -999 
data_efa$study_id <- as.numeric(data_efa$study_id)-1
data_efa$sex <- as.numeric(data_efa$sex)-1
data_efa$race <- as.numeric(data_efa$race)-1
head(data_efa)


# Create new 'efa' sub-folder within original working directory
dir.create(paste(myfolder,"/efa", sep=""))
setwd(paste(myfolder,"/efa", sep=""))
write.table(data_efa, "data_efa.dat", row.names=F, col.names=F, quote=F) 

# Determine number of processors to run in parallel
my_processors <- detectCores() - 1

# Create Mplus input files, in which EFAs are fit separately for each study
for(i in min(data_efa$study_id):max(data_efa$study_id)) {
input <- paste(
"title: EFA for study",i,"

data:
file = data_efa.dat;

variable: 
names = id study_id sex race x1-x9 learn;
usevariables = x1-x9;
categorical = x1-x9;
useobservations = study_id == ",i,";
missing = all (-999);

analysis:
type = efa 1 4;
estimator = wlsmv; !binary items
processors =",my_processors,";
", sep="")
write.table(input, paste("efa_study",i,".inp", sep=""), quote=F, row.names=F, col.names=F)
}

# Estimate models
runModels(replaceOutfile="never")



# Mplus input file for a final EFA for all observations (pooled across study)
input <- paste(
"title: EFA for all observations (pooled across study)

data:
file = data_efa.dat;

variable: 
names = id study_id sex race x1-x9 learn;
usevariables = x1-x9;
categorical = x1-x9;
missing = all (-999);

analysis:
type = efa 1 4;
estimator = wlsmv; !binary items
processors =",my_processors,";
", sep="")
write.table(input, "efa_allobs.inp", quote=F, row.names=F, col.names=F)

# Estimate final pooled model
runModels(target = paste(myfolder,"/efa/efa_allobs.inp", sep=""), replaceOutfile="never")


# Examine output from EFAs to determine dimensionality
out_efa_study0 <- readModels("efa_study0.out")
out_efa_study1 <- readModels("efa_study1.out")
out_efa_study2 <- readModels("efa_study2.out")
out_efa_study3 <- readModels("efa_study3.out")
out_efa_allobs <- readModels("efa_allobs.out")

# eigenvalues for each study
out_efa_study0$output[181:191] 
out_efa_study1$output[175:185] 
out_efa_study2$output[174:184]
out_efa_study3$output[175:185]

# eigenvalues for final pooled model with all observations
out_efa_allobs$output[171:181] #supports a 1-factor solution
#--------------------------------------------------------------------------------------------#





#--------------------------------------------------------------------------------------------#
# Step 2: Create Mplus input files for MNLFA model building, estimate models, examine output
#--------------------------------------------------------------------------------------------#
# First, prepare datafile for Mplus
data_mnlfa <- data_efa
head(data_mnlfa)

# Create new 'mnlfa' sub-folder within original working directory
dir.create(paste(myfolder,"/mnlfa", sep=""))
setwd(paste(myfolder,"/mnlfa", sep=""))
write.table(data_mnlfa, "data_mnlfa.dat", row.names=F, col.names=F, quote=F) 

# Baseline model allows covariates to moderate the factor mean & factor variance
baseline <- paste(
"title: Moderation of factor mean and variance

data:
file = data_mnlfa.dat;

variable: 
names = id study_id sex race x1-x9 learn; 
usevariables = study_id sex race x1-x9;
categorical = x1-x9; 
missing = all (-999);
constraint = study_id sex race;

analysis:
estimator = mlr;
link = logit;
processors = ", my_processors, ";
!estimator = wlsmv; !cannot be used with certain model constraints

model: 
Factor BY x1-x9;!measurement model

!allow covariates to moderate factor mean (linear function) 
Factor ON study_id sex race;
[Factor@0]; !constrain factor mean to zero to identify model

!factor variance implicitly set to one to identify model
Factor(factor_variance); !label for factor variance

model constraint:
new (f_study f_sex f_race); 

!allow covariates to moderate factor variance 
!(log-linear function to avoid negative values)
factor_variance = EXP(f_study*study_id + f_sex*sex + f_race*race);

output:
sampstat; 
svalues;
tech1;
", sep="")
write.table(baseline, "baseline.inp", quote=F, row.names=F, col.names=F)



# Item-models allows covariates to moderate the factor mean & factor variance, as
# well as the item intercept and item loading (done sequentially for each item)
for(i in 1:9) {
input <- paste(
"title: Moderation of factor mean and variance, as well as 
item intercept and factor loading for item x",i,"

data:
file = data_mnlfa.dat;

variable: 
names = id study_id sex race x1-x9 learn; 

usevariables = study_id sex race x1-x9;

categorical = x1-x9; 
missing = all (-999);
constraint = study_id sex race;

analysis:
estimator = mlr;
link = logit;
processors = ", my_processors, ";
!estimator = wlsmv; !cannot be used with certain model constraints

model: 
Factor BY x1-x9;!measurement model

!allow covariates to moderate factor mean (linear function) 
Factor ON study_id sex race;
[Factor@0]; !constrain factor mean to zero to identify model

!factor variance implicitly set to one to identify model
Factor(factor_variance); !label for factor variance

!allow covariates to moderate item i intercept
x",i," ON study_id sex race;

Factor BY x",i," (x",i,"_loading); !label for item i loading

model constraint:
new (f_study f_sex f_race); 
new (x_int x_study x_sex x_race);

!allow covariates to moderate factor variance 
!(log-linear function to avoid negative values)
factor_variance = EXP(f_study*study_id + f_sex*sex + f_race*race);

!allow covariates to moderate item i loading
x",i,"_loading = x_int + x_study*study_id + x_sex*sex + x_race*race;

output:
sampstat; 
svalues;
tech1;
", sep="")
write.table(input, paste("x",i,"_model.inp", sep=""), quote=F, row.names=F, col.names=F)
}


# Estimate all of the models
runModels(replaceOutfile="never") 
#--------------------------------------------------------------------------------------------#





#--------------------------------------------------------------------------------------------#
# Step 3: Conduct LRT between each item-model and the baseline model
#--------------------------------------------------------------------------------------------#
modelResults <- readModels(what=c("summaries", "parameters"))

modelLRT <- do.call("rbind", sapply(modelResults,"[", "summaries"))
modelLRT <- modelLRT[, c("Title", "Parameters", "LL", "LLCorrectionFactor")]
modelLRT$pval <- NA

for(i in 2:10) {
  L0 <- modelLRT[1,3]
  c0 <- modelLRT[1,4]
  p0 <- modelLRT[1,2]
  
  L1 <- modelLRT[i,3]
  c1 <- modelLRT[i,4]
  p1 <- modelLRT[i,2]
  
  cd <- ((p0*c0) - (p1*c1)) / (p0-p1)
  lrt <- -2 * (L0 - L1) / cd
  modelLRT[i,"pval"] <- dchisq(x=lrt, df=(p1-p0))
}
modelLRT[, 2:5]
# Results indicate each item-model (except for item-8) fits the data 
# significantly better than the baseline model (sig. p-values)

# As a result, item-models 1-7 & 9 will keep covariate moderation of an 
# intercept or loading, only if the covariate effect is significant
modelParms <- do.call("rbind", sapply(modelResults,"[", "parameters"))
baselineParms <- data.frame(do.call(cbind,modelParms[[1]]))
for(i in 2:10) {
  assign(paste0("x",i-1,"Parms", sep=""), data.frame(do.call(cbind,modelParms[[i]])))
}

# For the baseline model (factor moderation only):
# Rows 10-12 give moderation of factor mean 
# Rows 24-26 give moderation of factor variance
baselineParms
baselineParms[c(10:12, 24:26), ]

# For each item-model:
# Rows 10-12 give moderation of factor mean 
# Rows 27-29 give moderation of factor variance
# Rows 13-15 give moderation of item intercept (for item i)
# Rows 31-33 give moderation of factor loading (for item i)
x1Parms
x1Parms[c(13:15, 31:33), ] 
#for item 1:
#intercept moderators: sex
#slope moderators: race

x2Parms[c(13:15, 31:33), ]
#for item 2:
#intercept moderators: sex & race
#slope moderators: 

x3Parms[c(13:15, 31:33), ]
#for item 3:
#intercept moderators: 
#slope moderators: study & sex

x4Parms[c(13:15, 31:33), ]
#for item 4:
#intercept moderators: study, sex, & race
#slope moderators: study & race

x5Parms[c(13:15, 31:33), ]
#for item 5:
#intercept moderators: study & race
#slope moderators: study 

x6Parms[c(13:15, 31:33), ]
#for item 6:
#intercept moderators: study, sex, & race
#slope moderators: study

x7Parms[c(13:15, 31:33), ]
#for item 7:
#intercept moderators: study, sex, & race
#slope moderators: study & sex

# Item 8 model had non-significant LRT 
# No covariate moderation is explored as a result
# x8Parms[c(13:15, 31:33), ]

x9Parms[c(13:15, 31:33), ]
#for item 9:
#intercept moderators: race
#slope moderators: race
#--------------------------------------------------------------------------------------------#





#--------------------------------------------------------------------------------------------#
# Step 4: Remove remaining non-significant parameters, estimate next-to-last 
#         and final MNLFA model
#--------------------------------------------------------------------------------------------#
# Based on information above, construct the next-to-last model,
# keeping only significant moderators of item intercepts and 
# loadings (but always keeping moderation of factor mean and 
# variance, regardless of significance)
next_to_last_model <- paste(
"title: next-to-last MNLFA model

data:
file = data_mnlfa.dat;

variable: 
names = id study_id sex race
x1-x9 learn;

usevariables = study_id sex race x1-x9;

categorical = x1-x9;
missing = all (-999);
constraint = study_id sex race;

analysis:
estimator = mlr;
link = logit;
processors = ", my_processors, ";
!estimator = wlsmv; !cannot be used with certain model constraints

model:
factor BY x1*1 (x1_loading);
factor BY x2*1 (x2_loading);
factor BY x3*1 (x3_loading);
factor BY x4*1 (x4_loading);
factor BY x5*1 (x5_loading);
factor BY x6*1 (x6_loading);
factor BY x7*1 (x7_loading);
factor BY x8*1; !no label for item 8, no moderation
factor BY x9*1 (x9_loading);

!allow covariates to moderate factor mean (linear function) 
Factor ON study_id sex race;
[Factor@0]; !constrain factor mean to zero to identify model

!factor variance implicitly set to one to identify model
Factor(factor_variance); !label for factor variance

! Moderation of item intercepts (previously determined)
x1 ON sex;
x2 ON sex race;
!no moderation of x3 interceptr
x4 ON study_id sex race;
x5 ON study_id race;
x6 ON study_id sex race;
x7 ON study_id sex race;
!no moderation of x8;
x9 ON race;

model constraint:
NEW(f_study f_sex f_race);
NEW(int1 int3 int4 int5 int6 int7 int9); !intercepts for slope moderation equation

NEW(x1_race);
! no slope moderation of x2
NEW(x3_study x3_sex);
NEW(x4_study x4_race);
NEW(x5_study);
NEW(x6_study);
NEW(x7_study x7_sex);
! no moderation of x8
NEW(x9_race);

!allow covariates to moderate factor variance
!(log-linear function to avoid negative values)
factor_variance = EXP(f_study*study_id + f_sex*sex + f_race*race);

!allow covariates to moderate item loadings
x1_loading = int1 + x1_race*race;
! no slope moderation of x2
x3_loading = int3 + x3_study*study_id + x3_sex*sex;
x4_loading = int4 + x4_study*study_id + x4_race*race;
x5_loading = int5 + x5_study*study_id;
x6_loading = int6 + x6_study*study_id;
x7_loading = int7 + x7_study*study_id + x7_sex*sex;
! no moderation of x8
x9_loading = int9 + x9_race*race;

output:
sampstat;
svalues;
tech1;

", sep="")
write.table(next_to_last_model, "next_to_last_model.inp", quote=F, row.names=F, col.names=F)

# Estimate next-to-last model 
runModels("next_to_last_model.inp", replaceOutfile="never") 

# Read in the output of the next-to-last model 
out_next_to_last_model <- readModels("next_to_last_model.out") 

# Remove any non-significant moderation of item parameters
# (leave moderation of factor parameters regardless of significance)
out_next_to_last_model$parameters$unstandardized
# NOTE: This last pruning effort results in the final MNLFA model



# Build final MNLFA model 
final_MNLFA_model <- paste(
"title: final MNLFA model

data:
file = data_mnlfa.dat;

variable: 
names = id study_id sex race
x1-x9 learn;

idvariable = id;
usevariables = study_id sex race x1-x9;

categorical = x1-x9;
missing = all (-999);
constraint = study_id sex race;

analysis:
estimator = mlr;
link = logit;
processors = ", my_processors, ";
!estimator = wlsmv; !cannot be used with certain model constraints

model:
factor BY x1*1 (x1_loading);
factor BY x2*1 (x2_loading);
factor BY x3*1 (x3_loading);
factor BY x4*1 (x4_loading);
factor BY x5*1 (x5_loading);
factor BY x6*1 (x6_loading);
factor BY x7*1 (x7_loading);
factor BY x8*1; !no label for item 8, no moderation
factor BY x9*1 (x9_loading);

!allow covariates to moderate factor mean (linear function) 
Factor ON study_id sex race;
[Factor@0]; !constrain factor mean to zero to identify model

!factor variance implicitly set to one to identify model
Factor(factor_variance); !label for factor variance

! Moderation of item intercepts (previously determined)
x1 ON sex;
x2 ON sex;
x4 ON study_id;
x5 ON study_id race;
x6 ON study_id sex race;
x7 ON study_id race;

model constraint:
NEW(f_study f_sex f_race);
NEW(int1 int3 int4 int6 int7); !intercepts for slope moderation equation

NEW(x1_race);
NEW(x3_study);
NEW(x4_study);
NEW(x6_study);
NEW(x7_study);

!allow covariates to moderate factor variance
!(log-linear function to avoid negative values)
factor_variance = EXP(f_study*study_id + f_sex*sex + f_race*race);

!allow covariates to moderate item loadings
x1_loading = int1 + x1_race*race;
x3_loading = int3 + x3_study*study_id;
x4_loading = int4 + x4_study*study_id;
x6_loading = int6 + x6_study*study_id;
x7_loading = int7 + x7_study*study_id;

output:
sampstat;
svalues;
tech1;

savedata:
save = fscores; !save estimated factor scores
file = est_factor_scores.csv;

", sep="")
write.table(final_MNLFA_model, "final_MNLFA_model.inp", quote=F, row.names=F, col.names=F)

# Estimate final MNLFA model 
runModels("final_MNLFA_model.inp", replaceOutfile="never") 

# Read in the output of the final MNLFA model
out_final_MNLFA_model <- readModels("final_MNLFA_model.out") 
out_final_MNLFA_model$parameters$unstandardized
#--------------------------------------------------------------------------------------------#






#--------------------------------------------------------------------------------------------#
# Step 5: Merge estimated factor scores to be used in subsequent analyses
#--------------------------------------------------------------------------------------------#
est_factor_scores <- read.table("est_factor_scores.csv")
colnames(est_factor_scores) <- tolower(out_final_MNLFA_model$savedata_info$fileVarNames)
head(est_factor_scores)

# Keep just id and estimated factor score
est_factor_scores <- est_factor_scores[,c("id", "factor")]


# Merge estimated factor scores in with original data
data <- merge(data, est_factor_scores, by=c("id"))
head(data)


# Estimate the effect of aggressive behavior factor scores on 
# 'learns up to ability' (5-point Likert scale outcome) using 
# an ordinal logistic regression
data$learn_factor <- as.factor(data$learn)
ologit_model <- polr(learn_factor ~ factor, data, Hess=T)
summary(ologit_model)
exp(coef(ologit_model))
#--------------------------------------------------------------------------------------------#
