# by Linke, JO, Phd 06/01/21
# This code selects the relevant nuissance variables from the fMRIprep output that will be regressed from the time-series
#########################################################
# Get the necessary packages
#########################################################
if (!require("tidyverse")) {
  install.packages("tidyverse", dependencies = TRUE)
  library(tidyverse)
}

#########################################################
# Set directories
#########################################################
FPREPDIR <- ""  #Path to your fMRIprep output
LISTDIR <- ""   #Path to your participant list
REGDIR <- ""    #Path to where the output should be written

#########################################################
# Provide the IDs (numbers only) of the subjects you want to process
#########################################################
setwd(LISTDIR)
sublist <- readLines("FinalSample.txt")
sublist <- as.list(sublist)

for(subj in sublist){
  #Test if directory already exists, if not create directory
  ifelse(!dir.exists(paste(REGDIR,'/sub-s',subj,sep="")), dir.create(paste(REGDIR,'/sub-s',subj,sep="")), FALSE)
  datadir <- (paste(FPREPDIR,'/sub-s',subj,'/ses-1/func',sep=""))

  #Read all the different data files in
  setwd(datadir)
  rest <- read_delim(paste('sub-s',subj,'_ses-1_task-rest_desc-confounds_timeseries.tsv',sep=""),delim="\t", col_names = T)
  Melrest <- read_delim(paste('sub-s',subj,'_ses-1_task-rest_desc-MELODIC_mixing.tsv',sep=""),delim="\t", col_names = F)
  AROMArest <- read_delim(paste('sub-s',subj,'_ses-1_task-rest_AROMAnoiseICs.csv',sep=""),delim=",", col_names = F)
  AROMArest <- as.vector(as.numeric(AROMArest))

  #select variables from fmriprep
  rest1<-rest[c("csf", "white_matter","a_comp_cor_00","a_comp_cor_01", "a_comp_cor_02", "a_comp_cor_03", "cosine00", "cosine01", "cosine02",
                   "trans_x", "trans_y", "trans_z", "rot_x", "rot_y","rot_z")]
  #existence of cosine03 and non-steady-state differs between participants, so need to test, what's there
  if("cosine03" %in% colnames(rest1)) {rest1$cosine03<-rest1$cosine03;}
  if("non_steady_state_outlier00" %in% colnames(rest)) {rest1$non_steady_state_outlier00<-rest$non_steady_state_outlier00;}
  if("non_steady_state_outlier01" %in% colnames(rest)) {rest1$non_steady_state_outlier01<-rest$non_steady_state_outlier01;}
  if("non_steady_state_outlier02" %in% colnames(rest)) {rest1$non_steady_state_outlier02<-rest$non_steady_state_outlier02;}
  if("non_steady_state_outlier03" %in% colnames(rest)) {rest1$non_steady_state_outlier03<-rest$non_steady_state_outlier03;}
  if("non_steady_state_outlier04" %in% colnames(rest)) {rest1$non_steady_state_outlier04<-rest$non_steady_state_outlier04;}
  if("non_steady_state_outlier05" %in% colnames(rest)) {rest1$non_steady_state_outlier05<-rest$non_steady_state_outlier05;}
  if("non_steady_state_outlier06" %in% colnames(rest)) {rest1$non_steady_state_outlier06<-rest$non_steady_state_outlier06;}
  if("non_steady_state_outlier07" %in% colnames(rest)) {rest1$non_steady_state_outlier07<-rest$non_steady_state_outlier07;}
  
  
  #select AROMA components
  Melrest <- Melrest[AROMArest]
  #merge
  rest2 <- cbind(rest1,Melrest)
 
  #########################################################
  #Write to file
  #########################################################
  setwd(paste(REGDIR,"/sub-s",subj,"/",sep=""))

  OUT1 <- paste("sub-s",subj,"_ses-1_task-rest_AllConfounds.tsv",sep="")
  write.table(rest2, OUT1, quote=FALSE, sep="\t", col.names = F, row.names = F)
  }
