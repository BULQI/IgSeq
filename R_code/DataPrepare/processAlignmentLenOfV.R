#
#code to process the alignment v length and prepare the data
#    based on the alignment files (including germline genes and read seqences and alignment)
#   we get the total number of v length and number of sequence 
#     Note (updated 11/11/2020): the cloanaylst has been modified, so 
#           we now have one extra conditions for the alignment file, "NNN" meanng the discontinued part on the reads.
#           "nnn" means n-nucleotides in the alignment seq.
# =======================
#----11/12/2020
#now add the code to extract information from SimpleMarkedUAs.fasta file on the position of discontinued section
#       need to get where it starts relative to the starting of V region
#       and if there is no discontinued region,  using 9999 to indicate.
#============================
#---update 11/11/20202
#       modify the code the process the new alignment file including the "discontinued" region 
#       secondly, the simpleMarkedUAs file is pretty big and we will not download, but will run code on
#       the linux machine (theano.)
#       this file was copied from /LBSeq8 and start doing new pipeline analysis.
#--------------------------------------------------------
#
#-------update 3/18/2020
#add explanation of code in the "*.SimpleMarkedUAs.fasta"
# U :unkown/unsequenced nts, this is the part that has not been sequenced. that leads to a truncated sequence.
# A, B, 1,2,3 belong to v region(does not include the CDR3). still need to figure out the meanings of each code.
# V and n and D and J make up the CDR3???. ?? need to confirm.
# N is only showing on the sequenced read, not on the alignment part.!!! "n" is on alignment part, meaning n-nucleotides at either vd or dj junction.
#-------update 3/12/2020
#copied from WL03R2
#prepare the data for the analysis next
#---------

#here in this code, we go through the *.simpleMarkedUAs.fasta file to 
#determine the VRG results for the v length and matched length.
#we need them for the intraclonal diveristy for pairwised difference.
#

library(Biostrings)
library(stringr)
library(here)

#set working directory first. work directly on linux 
#setwd("E:\\feng\\LAB\\MSI\\MouseLungProject\\LBSeq8\\IgTotal");
#setwd("/projectnb/keplab/feng/Sequencing/130.64.74.72/200121-0384H_Feng_Feng/analysis_2/merge/");

#getting a list of files to be processed.
#now get all the directory names, these will be my sample names too.
#file.names<-list.files(path=".",pattern="*.SimpleMarkedUAs.fasta",full.names=F, recursive=F);

tissue<-c("MB", "SP")#tissue<-c("BM", "SP")
name_tissue<-c("MB", "SP")
data.in.dir<-"Data"
#BM IgG
n<-1
files<-list.files(path=here(data.in.dir,tissue[n]), 
	pattern=paste0(name_tissue[n],"(.)+R2_IgG.SimpleMarkedUAs.fasta"),full.names=T, recursive=F);

##go into each one and run stats
sMarked<-vector("list", length=length(files));
sDiscPosition<-vector("list", length=length(files));
index<-0;
sampleNames<-c()
for(s in files)
{
	cat("reading the file:", s,"\n")
	sname<-sub(paste0("./",tissue[n],"/"),"", s, fixed=T)
	sampleName<-sub("_.*", "", sname, fixed=F)
    
	sampleNames<-c(sampleNames, sampleName)
    s1<-readBStringSet(file.path(s))
	index<-index+1;
	 freq.U<-letterFrequency(s1,c("U", "A", "B", "1","2","3"))
	 
	 freq.U<-freq.U[seq(3,length(s1), by=3),]
	 
	 #also need to processed the sequence names.
	 sname<-names(s1)[seq(3,length(s1), by=3)];
	 #peel off the ig designations
	 sname<-sub("\\|.*", "", sname)
	 dft<-as.data.frame(freq.U)
	 names(dft)<-paste0("freq_", names(dft))
	 dft<-cbind(sname, dft)
	 names(dft)[1]<-"ReadID"
	 
	 #add a few more in there
	 dft$totalVBase<-apply(dft[,-1],1,sum)
	 dft$XVBases<-apply(dft[,-c(1:2, dim(dft)[2])],1,sum)
	 dft$sampleName<-sampleName;
	
    
    #get the discposition information
    #sname<-names(s1)[seq(1,length(s1), by=3)];
	 #peel off the ig designations
	 #sname<-sub("\\|.*", "", sname)
     dis<-as.character(s1[seq(1,length(s1), by=3)])
    dp<-str_locate(dis, "N+")
    colnames(dp)<-c("discPosStart", "discPosEnd")

    dft<-cbind(dft, dp)
    
    sMarked[[index]]<-dft;
	flush.console();
}

#reading in the conditions
#now in this run, we read the conditons from file
#conditions<-read.table(file="../sampleConditions.txt", sep="\t", header=T)


df.vlength<-NULL
i<-1
for(i in 1:length(sMarked))
{
	cat("doing sample ", i, "...........\n")
	temp<-sMarked[[i]]
	#temp$sampleName<-unqiue([i,"sampleName"]
	#temp$grouNames<-conditions[i,"groupNames"]
	#temp$cloneSize<-temp$"X.Members"/sum(temp$"X.Members")
	
	df.vlength<-rbind(df.vlength, temp)
    
	flush.console();
}
#conditions.bm<-conditions;
df.vlength.IgG<-df.vlength

###BM IgM
#BM IgG
n<-1
files<-list.files(paste0("./",tissue[n]), pattern=paste0(name_tissue[n],"(.)+R2_IgM.SimpleMarkedUAs.fasta"),full.names=T, recursive=F);

##go into each one and run stats
sMarked<-vector("list", length=length(files));
sDiscPosition<-vector("list", length=length(files));
index<-0;
sampleNames<-c()
for(s in files)
{
	cat("reading the file:", s,"\n")
	sname<-sub(paste0("./",tissue[n],"/"),"", s, fixed=T)
	sampleName<-sub("_.*", "", sname, fixed=F)
    
	sampleNames<-c(sampleNames, sampleName)
    s1<-readBStringSet(file.path(s))
	index<-index+1;
	 freq.U<-letterFrequency(s1,c("U", "A", "B", "1","2","3"))
	 
	 freq.U<-freq.U[seq(3,length(s1), by=3),]
	 
	 #also need to processed the sequence names.
	 sname<-names(s1)[seq(3,length(s1), by=3)];
	 #peel off the ig designations
	 sname<-sub("\\|.*", "", sname)
	 dft<-as.data.frame(freq.U)
	 names(dft)<-paste0("freq_", names(dft))
	 dft<-cbind(sname, dft)
	 names(dft)[1]<-"ReadID"
	 
	 #add a few more in there
	 dft$totalVBase<-apply(dft[,-1],1,sum)
	 dft$XVBases<-apply(dft[,-c(1:2, dim(dft)[2])],1,sum)
	 dft$sampleName<-sampleName;
     
    #get the discposition information
    #sname<-names(s1)[seq(1,length(s1), by=3)];
	 #peel off the ig designations
	 #sname<-sub("\\|.*", "", sname)
     dis<-as.character(s1[seq(1,length(s1), by=3)])
    dp<-str_locate(dis, "N+")
    colnames(dp)<-c("discPosStart", "discPosEnd")

    dft<-cbind(dft, dp)
    
	sMarked[[index]]<-dft;
	flush.console();
}

#reading in the conditions
#now in this run, we read the conditons from file
#conditions<-read.table(file="../sampleConditions.txt", sep="\t", header=T)


df.vlength<-NULL
i<-1
for(i in 1:length(sMarked))
{
	cat("doing sample ", i, "...........\n")
	temp<-sMarked[[i]]
	#temp$sampleName<-unqiue([i,"sampleName"]
	#temp$grouNames<-conditions[i,"groupNames"]
	#temp$cloneSize<-temp$"X.Members"/sum(temp$"X.Members")
	
	df.vlength<-rbind(df.vlength, temp)
	
	flush.console();
}
#conditions.bm<-conditions;
df.vlength.IgM<-df.vlength
BM.df.vlength<-list(IgG=df.vlength.IgG, IgM=df.vlength.IgM)


####-Spleen IgG
n<-2
files<-list.files(paste0("./",tissue[n]), pattern=paste0(name_tissue[n],"(.)+R2_IgG.SimpleMarkedUAs.fasta"),full.names=T, recursive=F);

##go into each one and run stats
sMarked<-vector("list", length=length(files));
sDiscPosition<-vector("list", length=length(files));
index<-0;
sampleNames<-c()
for(s in files)
{
	cat("reading the file:", s,"\n")
	sname<-sub(paste0("./",tissue[n],"/"),"", s, fixed=T)
	sampleName<-sub("_.*", "", sname, fixed=F)
    
	sampleNames<-c(sampleNames, sampleName)
    s1<-readBStringSet(file.path(s))
	index<-index+1;
	 freq.U<-letterFrequency(s1,c("U", "A", "B", "1","2","3"))
	 
	 freq.U<-freq.U[seq(3,length(s1), by=3),]
	 
	 #also need to processed the sequence names.
	 sname<-names(s1)[seq(3,length(s1), by=3)];
	 #peel off the ig designations
	 sname<-sub("\\|.*", "", sname)
	 dft<-as.data.frame(freq.U)
	 names(dft)<-paste0("freq_", names(dft))
	 dft<-cbind(sname, dft)
	 names(dft)[1]<-"ReadID"
	 
	 #add a few more in there
	 dft$totalVBase<-apply(dft[,-1],1,sum)
	 dft$XVBases<-apply(dft[,-c(1:2, dim(dft)[2])],1,sum)
	 dft$sampleName<-sampleName;
     
    #get the discposition information
    #sname<-names(s1)[seq(1,length(s1), by=3)];
	 #peel off the ig designations
	 #sname<-sub("\\|.*", "", sname)
     dis<-as.character(s1[seq(1,length(s1), by=3)])
    dp<-str_locate(dis, "N+")
    colnames(dp)<-c("discPosStart", "discPosEnd")

    dft<-cbind(dft, dp)
    
	sMarked[[index]]<-dft;
	flush.console();
}

#reading in the conditions
#now in this run, we read the conditons from file
#conditions<-read.table(file="../sampleConditions.txt", sep="\t", header=T)


df.vlength<-NULL
i<-1
for(i in 1:length(sMarked))
{
	cat("doing sample ", i, "...........\n")
	temp<-sMarked[[i]]
	#temp$sampleName<-unqiue([i,"sampleName"]
	#temp$grouNames<-conditions[i,"groupNames"]
	#temp$cloneSize<-temp$"X.Members"/sum(temp$"X.Members")
	
	df.vlength<-rbind(df.vlength, temp)
	
	flush.console();
}
#conditions.bm<-conditions;
df.vlength.IgG<-df.vlength

###SP IgM
n<-2
files<-list.files(paste0("./",tissue[n]), pattern=paste0(name_tissue[n],"(.)+R2_IgM.SimpleMarkedUAs.fasta"),full.names=T, recursive=F);

##go into each one and run stats
sMarked<-vector("list", length=length(files));
sDiscPosition<-vector("list", length=length(files));
index<-0;
sampleNames<-c()
for(s in files)
{
	cat("reading the file:", s,"\n")
	sname<-sub(paste0("./",tissue[n],"/"),"", s, fixed=T)
	sampleName<-sub("_.*", "", sname, fixed=F)
    
	sampleNames<-c(sampleNames, sampleName)
    s1<-readBStringSet(file.path(s))
	index<-index+1;
	 freq.U<-letterFrequency(s1,c("U", "A", "B", "1","2","3"))
	 
	 freq.U<-freq.U[seq(3,length(s1), by=3),]
	 
	 #also need to processed the sequence names.
	 sname<-names(s1)[seq(3,length(s1), by=3)];
	 #peel off the ig designations
	 sname<-sub("\\|.*", "", sname)
	 dft<-as.data.frame(freq.U)
	 names(dft)<-paste0("freq_", names(dft))
	 dft<-cbind(sname, dft)
	 names(dft)[1]<-"ReadID"
	 
	 #add a few more in there
	 dft$totalVBase<-apply(dft[,-1],1,sum)
	 dft$XVBases<-apply(dft[,-c(1:2, dim(dft)[2])],1,sum)
	 dft$sampleName<-sampleName;
     
    #get the discposition information
    #sname<-names(s1)[seq(1,length(s1), by=3)];
	 #peel off the ig designations
	 #sname<-sub("\\|.*", "", sname)
     dis<-as.character(s1[seq(1,length(s1), by=3)])
    dp<-str_locate(dis, "N+")
    colnames(dp)<-c("discPosStart", "discPosEnd")

    dft<-cbind(dft, dp)
    
	sMarked[[index]]<-dft;
	flush.console();
}

#reading in the conditions
#now in this run, we read the conditons from file
#conditions<-read.table(file="../sampleConditions.txt", sep="\t", header=T)


df.vlength<-NULL
i<-1
for(i in 1:length(sMarked))
{
	cat("doing sample ", i, "...........\n")
	temp<-sMarked[[i]]
	#temp$sampleName<-unqiue([i,"sampleName"]
	#temp$grouNames<-conditions[i,"groupNames"]
	#temp$cloneSize<-temp$"X.Members"/sum(temp$"X.Members")
	
	df.vlength<-rbind(df.vlength, temp)
	
	flush.console();
}
#conditions.bm<-conditions;
df.vlength.IgM<-df.vlength
SP.df.vlength<-list(IgG=df.vlength.IgG, IgM=df.vlength.IgM)
#save them
save( BM.df.vlength, SP.df.vlength, file="df.vlength.RData")
