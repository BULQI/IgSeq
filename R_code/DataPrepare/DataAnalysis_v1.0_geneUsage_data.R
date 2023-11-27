#R code to analyze the heavy chain sequencing data
#update 10/5/2021
#  ----writing functions to do reading of recom summary files and also add paramter to 
#           remove zero mutation sequences
#   in this module we do all mutation sequences. 
#           check the other file "DataAnalysis_v1_geneUsage_data_noZeroMu.r" 
#               for no zero mutation  sequences.
#update 12/2/2020
# add code to save the raw recombinationSummaries records (including out of frame ones.)
#---------------------
#update 11/16/2020
# fix a bug for saving sp df ig good recom sum data
#--------------------------
#update 11/11/2020
#add part to save the IG.RecSum as a data frame
#----update 10/5/202
#    add functions to rewrite the data so to save sample.usage.RData
#------update 10/4/2020
# copy this from "windowsD/feng/LAB/MSI/MouseLungProject/LBSeq8/IgTotal/" to do new analyses on new pipeline processed data
#----------- df
#----update 2/22/2020
# 	we also plot the IG reads in this module.
#	copied from WL03R2, IgSeq
#----------------------------
#Feng @BU -----1/30/2020.
# this one is copied from previous work WL03 analysis
#In here we read the data and save them. The analysis is in another file 
# DataAnalysis_v1.0_geneUsage_composition.RData
#-------------the below part is from previous version. read in caution.
#======================================================
#Plan 1) do gene usage to compare between groups
#	  2) 
#------------update 
#--5/28/2019
#	1)read the BM data to put them together with spleen data.
#-----------update
##6/3/2019
#1) new data. new analyzed data set. because the old processing have lots
#-----------update 
##6/11/2019
#1)add code to save the list and data frame of the recombination summary data.
# 
#==========================
library(plyr)
#library(ggplot2)
library(frLib)
library(here)
#first do gene usage.
#Read in the data
#need to figure out the file name

#set working directory first
#do SP first
#setwd("/home/feng/Feng/LAB/Wetzler_PorinB/wl03R2/newPipeline/")

#define the function to read the files for each piece of data.
#'@title reading recombinationSummaries files for processing.
#'@param files character array holding the files (no sub-directories, only the current folder)
#'@param noZeroMuFreq boolean used to get rid of zero mutation sequences. by default this is FALSE meaning
#'          we want to have zero mutation sequences.
#'
readSummaryRecom<-function (files, data.in.dir,tissue, noZeroMuFreq=FALSE)
{
        usages<-vector("list", length=length(files));
        index<-1
        s.summary.goodIG<-vector("list", length=length(files));
        s.summary.all<-vector("list", length=length(files));
        for (s in files)
        {
            cat( "Reading sample from ", s, ".....\n");
            #sname<-sub(paste0("./",tissue[1],"/"),"", s, fixed=T)
            sname<-sub("_.*", "", s, fixed=F)
            #if(sname=="Sample05"){
            #	#next;
            #	cat("go")
            #}
            s.summary<-read.table(here(data.in.dir, tissue,s), #file.path(s,paste0(sname,".RecombinationSummaries.txt")),
                        header=T, sep="\t", na.strings="NaN", comment.char="");
            #clean it up. Get rid of no-good ones.
            #get where we have CDR3  <--------this is not a good criterion.need to talk about it.
            if(noZeroMuFreq)  #<-remove the zero mutation ones
            {
                     s.summary<-s.summary[!is.na(s.summary$MuFreq)&s.summary$MuFreq>0,]
            }
            s.summary.CDR3<-s.summary[nchar(as.character(trimws(s.summary$CDR3)))>0&(trimws(s.summary$CDR3)!="No CDR3"),]
            
            #now get that done. need to do stats, first do the gene usage.
            #get the gene usage.
            s.summary.CDR3$sampleName<-sname;
            s.summary$sampleName<-sname;
            s.summary.goodIG[[index]]<-s.summary.CDR3;
            usage.allele<-s.summary.CDR3[,c("VGene", "DGene", "JGene")];
            usage.allele$sampleName<-sname
            #usage.allele$sampleName<-sub("_.+","",s);
            #usage.allele$sampleName<-sub("^./","",usage.allele$sampleName);
            #need to figure out the unique ones and ran stats
            usages[[index]]<-usage.allele
            s.summary.all[[index]]<-s.summary
            index<-index+1;
            flush.console();
        }
        #geneUsage.IgM<-usages
        #summary.goodIG.IgM.BM<-s.summary.goodIG
        #summary.all.IgM.BM<-s.summary.all
        list(usages=usages, goodIG=s.summary.goodIG, all=s.summary.all)
}


#read BM data first 
#go to BM to get recombinationSummaries
tissue<-c("BM", "SP")
data.in.dir<-"Data"
data.out.dir<-"Data"
#setwd("/home/feng/Feng/LAB/Wetzler_PorinB/wl03R2/newPipeline/BM")
files<-list.files(
        path=here(data.in.dir,"BM"),
        pattern="MB(.)+_IgM.RecombinationSummaries.txt",
        full.names=F, recursive=F);

#dirs<-list.dirs(path=".",full.names=T, recursive=F);
sname<-sub(paste0("./",tissue[1],"/"),"", files, fixed=T)
x<-readSummaryRecom(sname, data.in.dir, tissue[1])

geneUsage.IgM<-x[["usages"]]
summary.goodIG.IgM.BM<-x[["goodIG"]]
summary.all.IgM.BM<-x[["all"]]

###reading IgG
#setwd("/home/feng/Feng/LAB/Wetzler_PorinB/wl03R2/newPipeline/BM/")
files<-list.files(path=here(data.in.dir,"BM"), 
        pattern="MB(.)+_IgG.RecombinationSummaries.txt",full.names=F, recursive=F);
#dirs<-list.dirs(path=".",full.names=T, recursive=F);
x<-readSummaryRecom(files, data.in.dir, tissue[1])
geneUsage.IgG<-x[["usages"]]
summary.goodIG.IgG.BM<-x[["goodIG"]]
summary.all.IgG.BM<-x[["all"]]
BM<-list(IgM=geneUsage.IgM, IgG=geneUsage.IgG #, sgood.IgG=summary.goodIG.IgG, sgood.IgM=summary.goodIG.IgM)
                    )
#=========================================Spleen=================
#now doing Spleen
#read SP data first 
#go to SP to get recombinationSummaries
#setwd("/home/feng/Feng/LAB/Wetzler_PorinB/wl03R2/newPipeline/SP/")
files<-list.files(path=here(data.in.dir,"SP"), 
    pattern="SP(.)+_IgM.RecombinationSummaries.txt",full.names=F, recursive=F);
#dirs<-list.dirs(path=".",full.names=T, recursive=F);
x<-readSummaryRecom(files, data.in.dir, tissue[2])

geneUsage.IgM<-x[["usages"]]
summary.goodIG.IgM.SP<-x[["goodIG"]]
summary.all.IgM.SP<-x[["all"]]

###reading IgG
files<-list.files(path=here(data.in.dir,tissue[2])
    , pattern="SP(.)+_IgG.RecombinationSummaries.txt",
    full.names=F, recursive=F);
#dirs<-list.dirs(path=".",full.names=T, recursive=F);
x<-readSummaryRecom(files,data.in.dir,tissue[2])
##go into each one and run stats

geneUsage.IgG<-x[["usages"]]
summary.goodIG.IgG.SP<-x[["goodIG"]]
summary.all.IgG.SP<-x[["all"]]
SP<-list(IgM=geneUsage.IgM, IgG=geneUsage.IgG #, sgood.IgG=summary.goodIG.IgG, sgood.IgM=summary.goodIG.IgM)
                )
save(BM, SP ,file=here(data.out.dir,"geneUsage.Rdata"))
#load(file="geneUsage.Rdata");
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##updated on 5/09/2018. add the total number of good Igs
totalNum.BM.IgG<-rep(0, length(BM[["IgG"]]))
totalNum.BM.IgM<-rep(0, length(BM[["IgM"]]))
totalNum.SP.IgG<-rep(0, length(SP[["IgG"]]))
totalNum.SP.IgM<-rep(0, length(SP[["IgM"]]))
totalNum.BM.snames<-rep("", length(BM[["IgG"]]))
totalNum.SP.snames<-rep("", length(SP[["IgG"]]))
for(i in 1:length(totalNum.BM.IgG))
{
	totalNum.BM.IgG[i]<-dim(BM[["IgG"]][[i]])[1]
    totalNum.BM.IgM[i]<-dim(BM[["IgM"]][[i]])[1]
	totalNum.BM.snames[i]<-unique(BM[["IgG"]][[i]]$sampleName)
}
names(totalNum.BM.IgG)<-totalNum.BM.snames
names(totalNum.BM.IgM)<-totalNum.BM.snames
#SP
for(i in 1:length(totalNum.SP.IgG))
{
	totalNum.SP.IgG[i]<-dim(SP[["IgG"]][[i]])[1]
    totalNum.SP.IgM[i]<-dim(SP[["IgM"]][[i]])[1]
	totalNum.SP.snames[i]<-unique(SP[["IgG"]][[i]]$sampleName)
}
names(totalNum.SP.IgG)<-totalNum.SP.snames
names(totalNum.SP.IgM)<-totalNum.SP.snames

#now in this run, we read the conditons from file
#setwd("/home/feng/Feng/LAB/Wetzler_PorinB/wl03R2/newPipeline/")
    conditions<-read.table(file=here(data.in.dir,"sampleInfo.txt")
            , sep="\t", header=T)
    rownames(conditions)<-conditions$ID
    conditions$index.fileOrder<-0;
    conditions<-conditions[order(conditions$ID),]
    conditions[conditions$tissue=="Bone Marrow", "index.fileOrder"]<-order(conditions[conditions$tissue=="Bone Marrow", "snum"])
    conditions[conditions$tissue=="Spleen", "index.fileOrder"]<-order(conditions[conditions$tissue=="Spleen", "snum"])
#usages=BM[["IgG"]]; totalNum.Igs=totalNum.BM.IgG; s.summary.goodIG<-summary.goodIG.IgG.BM
makeUsageTable<-function(usages,    conditions,   totalNum.Igs,   s.summary.goodIG, s.summary.all )
{
    #now that we reading the gene usage data, do the stats
    #for vgene first
    #make a table 
    #we will run stats at gene level not allele level, so strip off
    usage.VGene<-NULL
    usage.DGene<-NULL
    usage.JGene<-NULL
    usage.VJGene<-NULL
    usage.VDJGene<-NULL
    IG.RecSum<-NULL;
    all.RecSum<-NULL;
    #vgenes<-NULL
    for(i in 1:length(usages))
    {
        vallele<-usages[[i]]$VGene
        dallele<-usages[[i]]$DGene
        jallele<-usages[[i]]$JGene
        #v_jallele<-usages[[i]]
        v<-sub("[*].+","",vallele) #<----get rid of allele number, V1-3.1*01, "*01" is the allele numer
        d<-sub("[*].+","",dallele)
        j<-sub("[*].+","",jallele)
        
        #vgenes<-unique(c(vgenes,unique(v)))
        tempv<-data.frame("sampleName"=usages[[i]]$sampleName,"VGene"=v)
        tempv<-count(tempv, vars="VGene")
        colnames(tempv)<-c("VGene",usages[[i]]$sampleName[1]);
        
        tempd<-data.frame("sampleName"=usages[[i]]$sampleName,"DGene"=d)
        tempd<-count(tempd, vars="DGene")
        colnames(tempd)<-c("DGene",usages[[i]]$sampleName[1]);
        
        tempj<-data.frame("sampleName"=usages[[i]]$sampleName,"JGene"=j)
        tempj<-count(tempj, vars="JGene")
        colnames(tempj)<-c("JGene",usages[[i]]$sampleName[1]);
        
        tempvj<-data.frame("sampleName"=usages[[i]]$sampleName,"VGene"=v,"JGene"=j )
        tempvj<-count(tempvj, vars=c("VGene","JGene"))
        colnames(tempvj)<-c("VGene","JGene",usages[[i]]$sampleName[1]);
        
        tempvdj<-data.frame("sampleName"=usages[[i]]$sampleName,"VGene"=v,"DGene"=d, "JGene"=j )
        tempvdj<-count(tempvdj, vars=c("VGene","DGene","JGene"))
        colnames(tempvdj)<-c("VGene","DGene","JGene",usages[[i]]$sampleName[1]);
        
        #build the data frame
        if(i==1)
        {
            usage.VGene<-tempv
            usage.DGene<-tempd
            usage.JGene<-tempj
            usage.VJGene<-tempvj;
            usage.VDJGene<-tempvdj;
        } else {
            usage.VGene<-merge(usage.VGene, tempv, all.x=T, all.y=T)
            usage.DGene<-merge(usage.DGene, tempd, all.x=T, all.y=T)
            usage.JGene<-merge(usage.JGene, tempj, all.x=T, all.y=T)
            usage.VJGene<-merge(usage.VJGene, tempvj, all.x=T, all.y=T)
            usage.VDJGene<-merge(usage.VDJGene, tempvdj, all.x=T, all.y=T);
        }
        
        #build a dataframe for good IGs for future usa.
        temp.IG<-s.summary.goodIG[[i]]
        temp.s<-s.summary.all[[i]]
        #temp.IG$sampleName<-usages[[i]]$sampleName
        IG.RecSum<-rbind(IG.RecSum, temp.IG);
        all.RecSum<-rbind(all.RecSum, temp.s)
        cat(i, "Roung, dim is ", dim(usage.VDJGene), "\n");
        flush.console();
    }

    #replace NA with zeros in order to do stats and math
    usage.noNA<-usage.VGene
    for(i in 1:dim(usage.noNA)[1])
    {
        usage.noNA[i,is.na(usage.noNA[i,])]<-0;
    }
    usage.VGene<-usage.noNA

    usage.noNA<-usage.VJGene
    for(i in 1:dim(usage.noNA)[1])
    {
        usage.noNA[i,is.na(usage.noNA[i,])]<-0;
    }
    usage.VJGene<-usage.noNA

    usage.noNA<-usage.VDJGene
    for(i in 1:dim(usage.noNA)[1])
    {
        usage.noNA[i,is.na(usage.noNA[i,])]<-0;
    }
    usage.VDJGene<-usage.noNA

    
    #determine the file order
    snames<-sort(colnames(usage.VGene)[-1]);
    #sample.num<-sub("[a-zA-Z]{2}", "", snames)
    #index<-order(as.integer(sample.num))
    #conditions[snames,]$index.fileOrder<-index;
    #now sort the samples to put them into groups
    #groupNames<-c("PBS","Ova", 
    #				"Ova+PorB","Ova+CpG",
    #				"Ova+Alum");

    #sampleName.group<-rep(groupNames,c(3,2,3,3,3));				
    #conditions<-data.frame("group.num"=c(1,1,1,2,2,3,3,3,4,4,4,5,5,5), 
    #			"sampleName"=totalNum.snames[c(1,8:14,2:7)],
    #			"groupNames"=c("PBS","PBS","PBS",
    #										"Ova", "Ova",
    #										"Ova+PorB","Ova+PorB","Ova+PorB",
    #										"Ova+CpG","Ova+CpG","Ova+CpG",
    #										"Ova+Alum","Ova+Alum","Ova+Alum"),
    #			index.fileOrder=c(1,8:14,2:7)
    #			);

    #IG.RecSum.sp<-IG.RecSum;
    #conditions.sp<-conditions
    #s.summary.goodIG.sp<-s.summary.goodIG;

    #save(IG.RecSum)

    #determine the index for the grouping the samples together
    #index<-c();

    #snames<-substr(snames, start=1, stop=nchar(snames)-5);
    #for( g.num in unique(conditions$group.num))
    #{
    #	#get the sample names to search through the s names
    #	index.cur<-which(conditions$group.num==g.num)
    #	temp.name<-conditions$sampleName[index.cur]
    #	conditions$index.fileOrder[index.cur]<-(match(temp.name, snames))
    #}



    totalNum.Igs<-totalNum.Igs[conditions[snames,]$index.fileOrder]
    #names(totalNum.Igs)<-totalNum.snames[conditions$index.fileOrder]
    #save(totalNum.Igs, file="WL03_totalNum.Igs.BySummaryGoodCR3.RData");			
    Gene.Seg.name.ind<-vector("list", 5);
    Gene.Seg.name.ind[[1]]<-usage.VGene[,1];
    Gene.Seg.name.ind[[2]]<-usage.DGene[,1];
    Gene.Seg.name.ind[[3]]<-usage.JGene[,1];
    Gene.Seg.name.ind[[4]]<-usage.VJGene[,c(1,2)];
    Gene.Seg.name.ind[[5]]<-usage.VDJGene[,c(1,2,3)];

    #sort the colnames in each output
    usage.VGene<-usage.VGene[, c(1,conditions[snames,]$index.fileOrder+1)]
    usage.DGene<-usage.DGene[, c(1,conditions[snames,]$index.fileOrder+1)]
    usage.JGene<-usage.JGene[, c(1,conditions[snames,]$index.fileOrder+1)]
    usage.VJGene<-usage.VJGene[,c(1,2,conditions[snames,]$index.fileOrder+2)]
    usage.VDJGene<-usage.VDJGene[,c(1,2,3,conditions[snames,]$index.fileOrder+3)]


    #now saving the data for the latter usage
    sample.usage<-list(VGene=usage.VGene,
                    DGene=usage.DGene,
                    JGene=usage.JGene,
                    VJGene=usage.VJGene,
                    VDJGene=usage.VDJGene,
                    conditions=conditions,
                    totalNum.Igs=totalNum.Igs,
                    Gene.Seg.name=Gene.Seg.name.ind);
    return (list(sample.usage, IG.RecSum, all.RecSum))
} # end of make usage table 

#now call the function to do job
#1) BM IgG
 #usages=BM[["IgG"]]; totalNum.Igs=totalNum.BM.IgG; s.summary.goodIG<-summary.goodIG.IgG.BM
 x<-makeUsageTable(usages=BM[["IgG"]], conditions=conditions, 
                totalNum.Igs=totalNum.BM.IgG, s.summary.goodIG=summary.goodIG.IgG.BM, 
                    s.summary.all=summary.all.IgG.BM)
      sample.usage.IgG<-x[[1]]
      RecSum.df.IgG<-x[[2]]
      RecSum.all.IgG<-x[[3]]
 x<-makeUsageTable(usages=BM[["IgM"]], conditions=conditions, 
                totalNum.Igs=totalNum.BM.IgM, s.summary.goodIG=summary.goodIG.IgM.BM, summary.all.IgM.BM)
    sample.usage.IgM<-x[[1]]
    RecSum.df.IgM<-x[[2]]
    RecSum.all.IgM<-x[[3]]
BM.susage<-list(IgG=sample.usage.IgG, IgM=sample.usage.IgM)
IG.RecSum.bm<-list(IgG=RecSum.df.IgG, IgM=RecSum.df.IgM)    
All.RecSum.bm<-list(IgG=RecSum.all.IgG, IgM=RecSum.all.IgM)

x <-makeUsageTable(usages=SP[["IgG"]], conditions=conditions, 
                totalNum.Igs=totalNum.SP.IgG, s.summary.goodIG=summary.goodIG.IgG.SP, summary.all.IgG.SP)
         sample.usage.IgG<-x[[1]]    
         RecSum.df.IgG<-x[[2]]   
         RecSum.all.IgG<-x[[3]]
x <-makeUsageTable(usages=SP[["IgM"]], conditions=conditions, 
                totalNum.Igs=totalNum.SP.IgM, s.summary.goodIG=summary.goodIG.IgM.SP, summary.all.IgM.SP)
                sample.usage.IgM<-x[[1]]
                RecSum.df.IgM<-x[[2]]
                RecSum.all.IgM<-x[[3]]
SP.susage<-list(IgG=sample.usage.IgG, IgM=sample.usage.IgM)
IG.RecSum.sp<-list(IgG=RecSum.df.IgG, IgM=RecSum.df.IgM)
All.RecSum.sp<-list(IgG=RecSum.all.IgG, IgM=RecSum.all.IgM)          

#save the data to the correct location.      
save(SP.susage, BM.susage, conditions, 
    file=here(data.in.dir,"Figure2/sample.geneUsage.RData"))

#save(s.summary.goodIG,  conditions,  file="./Ig.good.RecSum.RData");

save(IG.RecSum.bm, IG.RecSum.sp, conditions, 
    file=here(data.in.dir,"Figure4/Ig.good.RecSum.df.RData"));
save(IG.RecSum.bm, IG.RecSum.sp, conditions, 
    file=here(data.in.dir,"Figure6/Ig.good.RecSum.df.RData"));

save(All.RecSum.bm, All.RecSum.sp, conditions, 
    file=here(data.in.dir,"Figure4/All.RecSum.df.RData"))



