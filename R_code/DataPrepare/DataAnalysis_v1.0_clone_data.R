#R code to analyze the heavy chain sequencing data
#updated--- 10/18/2020
#copied from LBSeq8 and start doing new pipeline analysis. this is 
# done to read the clone data. !!!only partially updated on 10/18/202
# Feng
# -----updated 3/12/2020----
#Start doing LBSeq8
# copied from WL03R2. In this unit, we want to read the data for doing
# the clonal analysis
#--------------------------------------------------------------------------------------
# ---- update 2/26/2020
#clones. and do pre-processing as well.
# this one is built on top of the DataAnalysis_v1.0_clone.R
#Feng @BU -----5/30/2019.
#Plan 1) read the data file, Clones.txt, cloneAssignment.txt
#				summary.txt (?) not sure.
#	  2) we will get clone size, mutation freq and intraclonal divisity
#	  3) what else???
#--------------
## updated 6/5/2019
# 1) now use the newly processed data set to do the analysis. the old one was
#	not removed the non-Ig data
# 2)
#--------------
### updated 6/7/2019
# 1) add code to read the cloneAssignment.txt. so that we can do analysis using 
#		the individual sequences. (see the file in "DataAnalysis_v1.2_cloneAssign_Analysis.R")
#----------------
### update 6/10/2019
#  1)add code process the data into data frame. 
#			previously this was done in the code processing the data
#	2)also add the code the read mutation file from the VRG analysis.
#==========================
library(dplyr)
#library(ggplot2)
library(frLib)
library(here)
#library(Rtsne)
#first do gene usage.
#Read in the data
#need to figure out the file name

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#reading the Clone.txt first. 
	#	This is the clone information (no sequence information)
	#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#set working directory first
#setwd("E:\\feng\\LAB\\MSI\\MouseLungProject\\LBSeq8\\IgTotal");
#setwd("/home/feng/Feng/LAB/Wetzler_PorinB/wl03R2/newPipeline/")

data.in.dir<-"Data"
tissue<-c("BM", "SP")
name_tissue<-c("MB", "SP")
n<-1
files<-list.files(path=here(data.in.dir,tissue[n]), 
	pattern=paste0(name_tissue[n],"(.)+R2_IgG.Clones.txt"),full.names=T, recursive=F);

##go into each one and run stats
clones<-vector("list", length=length(files));
index<-1
sampleNames<-c()
s<-files[1]
for (s in files)
{
    cat( "Reading sample ", s, ".....\n");
    sname<-sub(paste0("./",tissue[n],"/"),"", s, fixed=T)
	sname<-basename(sub("_.*", "", sname, fixed=F))
    
	sampleNames<-c(sampleNames,sname)
	clone.each<-read.table(file=s,
				header=T, sep="\t", na.strings="NaN", comment.char="");
	
	clone.each$sampleName<-sname;
	
	clones[[index]]<-clone.each
	index<-index+1;
	flush.console();
}
clone.IgG<-clones

#now read IgM for bm
n<-1
files<-list.files(path=here(data.in.dir,tissue[n]),
 pattern=paste0(name_tissue[n],"(.)+R2_IgM.Clones.txt"),full.names=T, recursive=F);
#dirs<-list.dirs(path=".",full.names=T, recursive=F);
##go into each one and run stats
clones<-vector("list", length=length(files));
index<-1
sampleNames<-c()
s<-files[1]
for (s in files)
{
    cat( "Reading sample ", s, ".....\n");
    sname<-basename(sub(paste0("./",tissue[n],"/"),"", s, fixed=T))
	sname<-sub("_.*", "", sname, fixed=F)
    
	sampleNames<-c(sampleNames,sname)
	clone.each<-read.table(s,
				header=T, sep="\t", na.strings="NaN", comment.char="");
	
	clone.each$sampleName<-sname;
	
	clones[[index]]<-clone.each
	index<-index+1;
	flush.console();
}
clone.IgM<-clones
BM.clones<-list(IgM=clone.IgM, IgG=clone.IgG
                    )
########################################
#---------------------now doing SP tissues.
########################################
n<-2
files<-list.files(path=here(data.in.dir,tissue[n]), 
	pattern=paste0(name_tissue[n],"(.)+R2_IgG.Clones.txt"),full.names=T, recursive=F);
#dirs<-list.dirs(path=".",full.names=T, recursive=F);
##go into each one and run stats
clones<-vector("list", length=length(files));
index<-1
sampleNames<-c()
s<-files[1]
for (s in files)
{
    cat( "Reading sample ", s, ".....\n");
    sname<-basename(sub(paste0("./",tissue[n],"/"),"", s, 
    	fixed=T))
	sname<-sub("_.*", "", sname, fixed=F)
    
	sampleNames<-c(sampleNames,sname)
	clone.each<-read.table(s,
				header=T, sep="\t", na.strings="NaN", comment.char="");
	
	clone.each$sampleName<-sname;
	
	clones[[index]]<-clone.each
	index<-index+1;
	flush.console();
}
clone.IgG<-clones

#now read IgM for bm
n<-2
files<-list.files(path=here(data.in.dir,tissue[n]), 
	pattern=paste0(name_tissue[n],"(.)+R2_IgM.Clones.txt"),full.names=T, recursive=F);
#dirs<-list.dirs(path=".",full.names=T, recursive=F);
##go into each one and run stats
clones<-vector("list", length=length(files));
index<-1
sampleNames<-c()
s<-files[1]
for (s in files)
{
    cat( "Reading sample ", s, ".....\n");
    sname<-basename(sub(paste0("./",tissue[n],"/"),"", 
    	s, fixed=T))
	sname<-sub("_.*", "", sname, fixed=F)
    
	sampleNames<-c(sampleNames,sname)
	clone.each<-read.table(s,
				header=T, sep="\t", na.strings="NaN", comment.char="");
	
	clone.each$sampleName<-sname;
	
	clones[[index]]<-clone.each
	index<-index+1;
	flush.console();
}
clone.IgM<-clones
SP.clones<-list(IgM=clone.IgM, IgG=clone.IgG
                    )

#now in this run, we read the conditons from file
conditions<-read.table(file=here(data.in.dir,
		"sampleInfo.txt"), sep="\t", header=T)
                    
save(SP.clones, BM.clones, conditions,
	file=here(data.in.dir,"Figure5/clones.Rdata"))
#load("clones.Rdata")

#turn list in data frame
#BM IgM 
clones<-BM.clones[[1]]
#conditions<-BM.Clones[[2]]
df.clones<-NULL
i<-1
for(i in 1:length(clones))
{
	cat("doing sample ", i, "...........\n")
	temp<-clones[[i]][,c("CloneID", "CDR3Length", "S.Index", "X.Members", "MeanMuFreq","sampleName")]
	#temp$sampleName<-conditions[i,"sampleName"]
	#temp$grouNames<-conditions[i,"groupNames"]
	temp$cloneSize<-temp$"X.Members"/sum(temp$"X.Members")
	
	if(i==1){
		df.clones<-temp;
	} else {
		df.clones<-rbind(df.clones, temp)
	}
	flush.console();
}
df.clones.IgM<-df.clones

#BM IgG
clones<-BM.clones[[2]]
#conditions<-BM.Clones[[2]]
df.clones<-NULL
i<-1
for(i in 1:length(clones))
{
	cat("doing sample ", i, "...........\n")
	temp<-clones[[i]][,c("CloneID", "CDR3Length", "S.Index", "X.Members", "MeanMuFreq","sampleName")]
	#temp$sampleName<-conditions[i,"sampleName"]
	#temp$grouNames<-conditions[i,"groupNames"]
	temp$cloneSize<-temp$"X.Members"/sum(temp$"X.Members")
	
	if(i==1){
		df.clones<-temp;
	} else {
		df.clones<-rbind(df.clones, temp)
	}
	flush.console();
}
df.clones.IgG<-df.clones
BM.df.clones<-list(IgG=df.clones.IgG, IgM=df.clones.IgM)

#SP IgM 
clones<-SP.clones[[1]]
#conditions<-BM.Clones[[2]]
df.clones<-NULL
i<-1
for(i in 1:length(clones))
{
	cat("doing sample ", i, "...........\n")
	temp<-clones[[i]][,c("CloneID", "CDR3Length", "S.Index", "X.Members", "MeanMuFreq","sampleName")]
	#temp$sampleName<-conditions[i,"sampleName"]
	#temp$grouNames<-conditions[i,"groupNames"]
	temp$cloneSize<-temp$"X.Members"/sum(temp$"X.Members")
	
	if(i==1){
		df.clones<-temp;
	} else {
		df.clones<-rbind(df.clones, temp)
	}
	flush.console();
}
df.clones.IgM<-df.clones

#SP IgG
clones<-SP.clones[[2]]
#conditions<-BM.Clones[[2]]
df.clones<-NULL
i<-1
for(i in 1:length(clones))
{
	cat("doing sample ", i, "...........\n")
	temp<-clones[[i]][,c("CloneID", "CDR3Length", "S.Index", "X.Members", "MeanMuFreq","sampleName")]
	#temp$sampleName<-conditions[i,"sampleName"]
	#temp$grouNames<-conditions[i,"groupNames"]
	temp$cloneSize<-temp$"X.Members"/sum(temp$"X.Members")
	
	if(i==1){
		df.clones<-temp;
	} else {
		df.clones<-rbind(df.clones, temp)
	}
	flush.console();
}
df.clones.IgG<-df.clones
SP.df.clones<-list(IgG=df.clones.IgG, IgM=df.clones.IgM)

#now in this run, we read the conditons from file
#conditions<-read.table(file="../sampleConditions.txt", sep="\t", header=T)

#determine the file order
conditions$index.fileOrder=0;
df.clones<-BM.df.clones[[1]]
snames<-unique(df.clones$sampleName);
sample.num<-gsub("MB", "", snames)
index<-order(as.integer(basename(sample.num)))
conditions[conditions$tissue=="Bone Marrow", "index.fileOrder"]<-index;

df.clones<-SP.df.clones[[1]]
snames<-unique(df.clones$sampleName);
sample.num<-gsub("SP", "", snames)
index<-order(as.integer(basename(sample.num)))
conditions[conditions$tissue=="Spleen", "index.fileOrder"]<-index;

#BM.Clones<-list(clones=clones, condition=conditions)
#conditions<-conditions;
 save(SP.df.clones,BM.df.clones, conditions, 
 	file=here(data.in.dir,"Figure5/clones.df.RData"))
  save(SP.df.clones,BM.df.clones, conditions, 
 	file=here(data.in.dir,"Figure6/clones.df.RData"))

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#now we want to read the cloneAssignment file then, we can know each individual+ 
	#sequence information. here more about individual mutation frequency.          +
	#*******************************************************************************
#set working directory first
#setwd("E:\\feng\\LAB\\MSI\\MouseLungProject\\LBSeq8\\IgTotal");

n<-1
files<-list.files(path=here(data.in.dir,tissue[n]), 
	pattern=paste0(name_tissue[n],"(.)+R2_IgG.CloneAssignments.txt"),full.names=T, recursive=F);
#dirs<-list.dirs(path=".",full.names=T, recursive=F);

##go into each one and run stats
clone.Assignment<-vector("list", length=length(files));
sampleNames<-c()
index<-1
for (s in files)
{
	cat( "Reading sample ", s, ".....\n");
	sname<-basename(sub(paste0("./",tissue[n],"/"),"", 
		s, fixed=T))
	sname<-sub("_.*", "", sname, fixed=F)
    
	#if(sname=="Sample05"){
	#	next;
	#}
	sampleNames<-c(sampleNames,sname)
	clone.each<-read.table(s, #file.path(s, "/CloneAssignments.txt"),
				header=T, sep="\t", na.strings="NaN", comment.char="");
	
	clone.each$sampleName<-sname;
	
	clone.Assignment[[index]]<-clone.each
	index<-index+1;
	flush.console();
}

clones<-clone.Assignment
#conditions<-SP.assignment[[2]]
df.clones<-NULL
i<-1
for(i in 1:length(clones))
{
	cat("doing sample ", i, "...........\n")
	temp<-clones[[i]]
	#temp$sampleName<-conditions[i,"sampleName"]
	#temp$grouNames<-conditions[i,"groupNames"]
	#temp$cloneSize<-temp$"X.Members"/sum(temp$"X.Members")
	
	if(i==1){
		df.clones<-temp;
	} else {
		df.clones<-rbind(df.clones, temp)
	}
	flush.console();
}
df.cloneAssign.IgG<-df.clones;
#conditions.assign.sp<-conditions;

####clone assign IgM, BM
n<-1
files<-list.files(path=here(data.in.dir,tissue[n]), 
	pattern=paste0(name_tissue[n],"(.)+R2_IgM.CloneAssignments.txt"),full.names=T, recursive=F);
#dirs<-list.dirs(path=".",full.names=T, recursive=F);

##go into each one and run stats
clone.Assignment<-vector("list", length=length(files));
sampleNames<-c()
index<-1
for (s in files)
{
	cat( "Reading sample ", s, ".....\n");
	sname<-basename(sub(paste0("./",tissue[n],"/"),"", 
		s, fixed=T))
	sname<-sub("_.*", "", sname, fixed=F)
    
	#if(sname=="Sample05"){
	#	next;
	#}
	sampleNames<-c(sampleNames,sname)
	clone.each<-read.table(s, #file.path(s, "/CloneAssignments.txt"),
				header=T, sep="\t", na.strings="NaN", comment.char="");
	
	clone.each$sampleName<-sname;
	
	clone.Assignment[[index]]<-clone.each
	index<-index+1;
	flush.console();
}

clones<-clone.Assignment
#conditions<-SP.assignment[[2]]
df.clones<-NULL
i<-1
for(i in 1:length(clones))
{
	cat("doing sample ", i, "...........\n")
	temp<-clones[[i]]
	#temp$sampleName<-conditions[i,"sampleName"]
	#temp$grouNames<-conditions[i,"groupNames"]
	#temp$cloneSize<-temp$"X.Members"/sum(temp$"X.Members")
	
	if(i==1){
		df.clones<-temp;
	} else {
		df.clones<-rbind(df.clones, temp)
	}
	flush.console();
}
df.cloneAssign.IgM<-df.clones;

BM.df.cloneAssign<-list(IgM=df.cloneAssign.IgM, IgG=df.cloneAssign.IgG)

#######3---Spleen IgG
n<-2
files<-list.files(path=here(data.in.dir,tissue[n]), 
	pattern=paste0(name_tissue[n],"(.)+R2_IgG.CloneAssignments.txt"),full.names=T, recursive=F);
#dirs<-list.dirs(path=".",full.names=T, recursive=F);

##go into each one and run stats
clone.Assignment<-vector("list", length=length(files));
sampleNames<-c()
index<-1
for (s in files)
{
	cat( "Reading sample ", s, ".....\n");
	sname<-basename(sub(paste0("./",tissue[n],"/"),"", 
		s, fixed=T))
	sname<-sub("_.*", "", sname, fixed=F)
    
	#if(sname=="Sample05"){
	#	next;
	#}
	sampleNames<-c(sampleNames,sname)
	clone.each<-read.table(s, #file.path(s, "/CloneAssignments.txt"),
				header=T, sep="\t", na.strings="NaN", comment.char="");
	
	clone.each$sampleName<-sname;
	
	clone.Assignment[[index]]<-clone.each
	index<-index+1;
	flush.console();
}

clones<-clone.Assignment
#conditions<-SP.assignment[[2]]
df.clones<-NULL
i<-1
for(i in 1:length(clones))
{
	cat("doing sample ", i, "...........\n")
	temp<-clones[[i]]
	#temp$sampleName<-conditions[i,"sampleName"]
	#temp$grouNames<-conditions[i,"groupNames"]
	#temp$cloneSize<-temp$"X.Members"/sum(temp$"X.Members")
	
	if(i==1){
		df.clones<-temp;
	} else {
		df.clones<-rbind(df.clones, temp)
	}
	flush.console();
}
df.cloneAssign.IgG<-df.clones;
#conditions.assign.sp<-conditions;

####clone assign IgM, BM
n<-2
files<-list.files(path=here(data.in.dir,tissue[n]), 
	pattern=paste0(name_tissue[n],"(.)+R2_IgM.CloneAssignments.txt"),full.names=T, recursive=F);
#dirs<-list.dirs(path=".",full.names=T, recursive=F);

##go into each one and run stats
clone.Assignment<-vector("list", length=length(files));
sampleNames<-c()
index<-1
for (s in files)
{
	cat( "Reading sample ", s, ".....\n");
	sname<-basename(sub(paste0("./",tissue[n],"/"),
		"", s, fixed=T))
	sname<-sub("_.*", "", sname, fixed=F)
    
	#if(sname=="Sample05"){
	#	next;
	#}
	sampleNames<-c(sampleNames,sname)
	clone.each<-read.table(s, #file.path(s, "/CloneAssignments.txt"),
				header=T, sep="\t", na.strings="NaN", comment.char="");
	
	clone.each$sampleName<-sname;
	
	clone.Assignment[[index]]<-clone.each
	index<-index+1;
	flush.console();
}

clones<-clone.Assignment
#conditions<-SP.assignment[[2]]
df.clones<-NULL
i<-1
for(i in 1:length(clones))
{
	cat("doing sample ", i, "...........\n")
	temp<-clones[[i]]
	#temp$sampleName<-conditions[i,"sampleName"]
	#temp$grouNames<-conditions[i,"groupNames"]
	#temp$cloneSize<-temp$"X.Members"/sum(temp$"X.Members")
	
	if(i==1){
		df.clones<-temp;
	} else {
		df.clones<-rbind(df.clones, temp)
	}
	flush.console();
}
df.cloneAssign.IgM<-df.clones;

SP.df.cloneAssign<-list(IgM=df.cloneAssign.IgM, IgG=df.cloneAssign.IgG)

save(SP.df.cloneAssign, BM.df.cloneAssign,conditions, 
	file=here(data.in.dir,"Figure6/cloneAssignment_df.RData"))


	##########################################
	#reading the mutation file from the VRG analysis
	#########################################
	#set working directory first
n<-1
files<-list.files(path=here(data.in.dir,tissue[n]), 
	pattern=paste0(name_tissue[n],"(.)+R2_IgM.Mutations.txt"),full.names=T, recursive=F);
#dirs<-list.dirs(path=".",full.names=T, recursive=F);

##go into each one and run stats
mutations<-vector("list", length=length(files));
sampleNames<-c()
index<-1
for (s in files)
{
	cat( "Reading sample ", s, ".....\n");
	sname<-basename(sub(paste0("./",tissue[n],"/"),"", 
		s, fixed=T))
	sname<-sub("_.*", "", sname, fixed=F)

	sampleNames<-c(sampleNames,sname)
	mutation.each<-read.table(s,
				header=T, sep="\t", na.strings="NaN", comment.char="");
	
	mutation.each$sampleName<-sname;
	
	mutations[[index]]<-mutation.each
	index<-index+1;
	flush.console();
}

#conditions<-SP.assignment[[2]]
df.mutations<-NULL
i<-1
for(i in 1:length(mutations))
{
	cat("doing sample ", i, "...........\n")
	temp<-mutations[[i]]
	#temp$sampleName<-conditions[i,"sampleName"]
	#temp$grouNames<-conditions[i,"groupNames"]
	#temp$cloneSize<-temp$"X.Members"/sum(temp$"X.Members")
	
	if(i==1){
		df.mutations<-temp;
	} else {
		df.mutations<-rbind(df.mutations, temp)
	}
	flush.console();
}
df.mutations.IgM<-df.mutations;

#----------BM IgG
n<-1
files<-list.files(path=here(data.in.dir,tissue[n]), 
	pattern=paste0(name_tissue[n],"(.)+R2_IgG.Mutations.txt"),full.names=T, recursive=F);
#dirs<-list.dirs(path=".",full.names=T, recursive=F);

##go into each one and run stats
mutations<-vector("list", length=length(files));
sampleNames<-c()
index<-1
for (s in files)
{
	cat( "Reading sample ", s, ".....\n");
	sname<-basename(sub(paste0("./",tissue[n],"/"),"",
	 s, fixed=T))
	sname<-sub("_.*", "", sname, fixed=F)

	sampleNames<-c(sampleNames,sname)
	mutation.each<-read.table(s,
				header=T, sep="\t", na.strings="NaN", comment.char="");
	
	mutation.each$sampleName<-sname;
	
	mutations[[index]]<-mutation.each
	index<-index+1;
	flush.console();
}

#conditions<-SP.assignment[[2]]
df.mutations<-NULL
i<-1
for(i in 1:length(mutations))
{
	cat("doing sample ", i, "...........\n")
	temp<-mutations[[i]]
	#temp$sampleName<-conditions[i,"sampleName"]
	#temp$grouNames<-conditions[i,"groupNames"]
	#temp$cloneSize<-temp$"X.Members"/sum(temp$"X.Members")
	
	if(i==1){
		df.mutations<-temp;
	} else {
		df.mutations<-rbind(df.mutations, temp)
	}
	flush.console();
}
df.mutations.IgG<-df.mutations;
BM.df.mutations<-list(IgG=df.mutations.IgG, IgM=df.mutations.IgM)

#####Spleen IgM
n<-2
files<-list.files(path=here(data.in.dir,tissue[n]), 
	pattern=paste0(name_tissue[n],"(.)+R2_IgM.Mutations.txt"),full.names=T, recursive=F);
#dirs<-list.dirs(path=".",full.names=T, recursive=F);

##go into each one and run stats
mutations<-vector("list", length=length(files));
sampleNames<-c()
index<-1
for (s in files)
{
	cat( "Reading sample ", s, ".....\n");
	sname<-basename(sub(paste0("./",tissue[n],"/"),
		"", s, fixed=T))
	sname<-sub("_.*", "", sname, fixed=F)

	sampleNames<-c(sampleNames,sname)
	mutation.each<-read.table(s,
				header=T, sep="\t", na.strings="NaN", comment.char="");
	
	mutation.each$sampleName<-sname;
	
	mutations[[index]]<-mutation.each
	index<-index+1;
	flush.console();
}

#conditions<-SP.assignment[[2]]
df.mutations<-NULL
i<-1
for(i in 1:length(mutations))
{
	cat("doing sample ", i, "...........\n")
	temp<-mutations[[i]]
	#temp$sampleName<-conditions[i,"sampleName"]
	#temp$grouNames<-conditions[i,"groupNames"]
	#temp$cloneSize<-temp$"X.Members"/sum(temp$"X.Members")
	
	if(i==1){
		df.mutations<-temp;
	} else {
		df.mutations<-rbind(df.mutations, temp)
	}
	flush.console();
}
df.mutations.IgM<-df.mutations;

#----------Spleen  IgG
n<-2
files<-list.files(path=here(data.in.dir,tissue[n]), 
	pattern=paste0(name_tissue[n],"(.)+R2_IgG.Mutations.txt"),full.names=T, recursive=F);
#dirs<-list.dirs(path=".",full.names=T, recursive=F);

##go into each one and run stats
mutations<-vector("list", length=length(files));
sampleNames<-c()
index<-1
for (s in files)
{
	cat( "Reading sample ", s, ".....\n");
	sname<-basename(sub(paste0("./",tissue[n],"/"),
		"", s, fixed=T))
	sname<-sub("_.*", "", sname, fixed=F)

	sampleNames<-c(sampleNames,sname)
	mutation.each<-read.table(s,
				header=T, sep="\t", na.strings="NaN", comment.char="");
	
	mutation.each$sampleName<-sname;
	
	mutations[[index]]<-mutation.each
	index<-index+1;
	flush.console();
}

#conditions<-SP.assignment[[2]]
df.mutations<-NULL
i<-1
for(i in 1:length(mutations))
{
	cat("doing sample ", i, "...........\n")
	temp<-mutations[[i]]
	#temp$sampleName<-conditions[i,"sampleName"]
	#temp$grouNames<-conditions[i,"groupNames"]
	#temp$cloneSize<-temp$"X.Members"/sum(temp$"X.Members")
	
	if(i==1){
		df.mutations<-temp;
	} else {
		df.mutations<-rbind(df.mutations, temp)
	}
	flush.console();
}
df.mutations.IgG<-df.mutations;
SP.df.mutations<-list(IgG=df.mutations.IgG, IgM=df.mutations.IgM)

#save(BM.mutations, SP.mutations, conditions.mutations.sp, conditions.mutations.bm, file="../../mutations.RData");
save(BM.df.mutations, SP.df.mutations, conditions,  
	file=here(data.in.dir,"Figure6/mutations_df.RData"));