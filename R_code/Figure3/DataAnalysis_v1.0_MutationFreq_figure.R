#R code to analyze the heavy chain sequencing data, mutation frequency
#the mutation analysis for manuscript Figure 3.
#
## Make sure you have read ./ReadMe.txt for instructions
# of running scripts and preparing data.
#==========================
library(plyr)
library(frLib)
library(emmeans)
library(multcompView)
library(multcomp) #for doing coxbox transformation testing
library(ggplot2)
library(ggpubr)
library(here)
library(car)

output.dir<-"R_code/Figure3"
data.dir<-"Data"
data.figure3.dir<-"Data/Figure3"

#read BM data first 
#go to BM to get recombinationSummaries
tissue<-c("SP", "BM")
isotype<-c("IgG", "IgM")
#files<-list.files(paste0("./",tissue[1]), pattern=paste0("MB(.)+_","isotype[1]",".RecombinationSummaries.txt"),full.names=T, recursive=F);
#sname<-sub(paste0("./",tissue[1],"/"),"", files, fixed=T)
# sname<-sub("_.*", "", sname, fixed=F)

readSummaryFile<-function(files, sname)
{##go into each one and run stats
        MuFreq<-vector("list", length=length(files));
        index<-1
        for (s in files)
        {
            cat( "Reading sample from ", s, ".....\n");
            #sname<-sub(paste0("./",tissue[1],"/"),"", s, fixed=T)
            #sname<-sub("_.*", "", sname, fixed=F)
            #if(sname=="Sample05"){
            # #next;
            # cat("go")
            #}
            s.summary<-read.table(s, #file.path(s,paste0(sname,".RecombinationSummaries.txt")),
                        header=T, sep="\t", na.strings="NaN", comment.char="");
            #s.summary<-read.table(s,
            #     header=T, sep="\t", na.strings="NaN", comment.char="");
            #clean it up. Get rid of no-good ones.
            #get where we have CDR3  <--------this is not a good criterion.need to talk about it.
            s.summary.CDR3<-s.summary[nchar(as.character(trimws(s.summary$CDR3)))>0&(trimws(s.summary$CDR3)!="No CDR3"),]
            
            #add another criterion 
            s.summary.CDR3<-s.summary.CDR3[s.summary.CDR3$X.VBases*1.1>s.summary.CDR3$DiscntLen,]
            #now get that done. need to do stats, first do the gene usage.
            #get the gene usage.
            #sampleName<-sub("_.+","",s);
            #sampleName<-sub("^./","",sampleName);
            sampleName<-sname[index]
            if(dim(s.summary.CDR3)[1]==0)
            {
                MuFreq.temp<-s.summary
                MuFreq.temp$sampleName<-sampleName
                MuFreq.temp<-MuFreq.temp[c(),c("MuFreq","sampleName")]
                
            } else {
                MuFreq.temp<-data.frame("MuFreq"=s.summary.CDR3[,c("MuFreq")], "sampleName"=sampleName);
            }
            #MuFreq.temp$sampleName<-s;
            #need to figure out the unique ones and ran stats
            MuFreq[[index]]<-MuFreq.temp
            index<-index+1;
            flush.console();    
        }
    return(MuFreq)
}

files<-list.files(path=here(data.dir,tissue[1]), pattern=paste0(tissue[1], "(.)+_", isotype[1],".RecombinationSummaries.txt"),full.names=T, recursive=F);
sname<-basename(files)#sub(paste0("./",tissue[1],"/"),"", files, fixed=T)
  sname<-sub("_.*", "", sname, fixed=F)
MuFreq.SP.IgG<-readSummaryFile(files=files, sname=sname)

files<-list.files(path=here(data.dir,tissue[1]), pattern=paste0(tissue[1], "(.)+_", isotype[2],".RecombinationSummaries.txt"),full.names=T, recursive=F);
sname<-basename(files)#sub(paste0("./",tissue[1],"/"),"", files, fixed=T)
  sname<-sub("_.*", "", sname, fixed=F)
MuFreq.SP.IgM<-readSummaryFile(files=files, sname=sname)

files<-list.files(path=here(data.dir,tissue[2]), pattern=paste0("MB", "(.)+_", isotype[1],".RecombinationSummaries.txt"),full.names=T, recursive=F);
sname<-basename(files)#sub(paste0("./",tissue[2],"/"),"", files, fixed=T)
  sname<-sub("_.*", "", sname, fixed=F)
MuFreq.BM.IgG<-readSummaryFile(files=files, sname=sname)

files<-list.files(path=here(data.dir,tissue[2]), pattern=paste0("MB", "(.)+_", isotype[2],".RecombinationSummaries.txt"),full.names=T, recursive=F);
sname<-basename(files) #sub(paste0("./",tissue[2],"/"),"", files, fixed=T)
  sname<-sub("_.*", "", sname, fixed=F)
MuFreq.BM.IgM<-readSummaryFile(files=files, sname=sname)

conditions<-read.table(file=here(data.dir,"sampleInfo.txt"), sep="\t", header=T)
    rownames(conditions)<-conditions$ID
    conditions$index.fileOrder<-0;
    conditions<-conditions[order(conditions$ID),]
    conditions[conditions$tissue=="Bone Marrow", "index.fileOrder"]<-order(conditions[conditions$tissue=="Bone Marrow", "snum"])
    conditions[conditions$tissue=="Bone Marrow", ]<-conditions[order(conditions[conditions$tissue=="Bone Marrow", "snum"]),]
    conditions[conditions$tissue=="Spleen", "index.fileOrder"]<-order(conditions[conditions$tissue=="Spleen", "snum"])
    conditions[conditions$tissue=="Spleen", ]<-conditions[conditions$tissue=="Spleen", ][order(conditions[conditions$tissue=="Spleen", "snum"]),]
    rownames(conditions)<-conditions$ID
    conditions$treatment<-factor(conditions$treatment, levels=c("PBS", "OVA", "OVA+PorB", "OVA+CpG", "OVA+Alum"))
#save(MuFreq.BM.IgG,MuFreq.BM.IgM,MuFreq.SP.IgG,MuFreq.SP.IgM, conditions,file="MuFreq.Rdata")
#load(file="MuFreq.Rdata")


#############+++++++++++
### BM IgM  and IgG together                                    |
########################
MuFreq.group<-NULL
MuFreq<-MuFreq.BM.IgM
#vgenes<-NULL
for(i in 1:length(MuFreq))
{
  mf<-MuFreq[[i]]$MuFreq
  
  temp<-data.frame("sampleName"=MuFreq[[i]]$sampleName,"MuFreq"=mf)
  
  #build the data frame
  if(i==1)
  {
    MuFreq.group<-temp
  } else {
    MuFreq.group<-rbind(MuFreq.group, temp)#, all.x=T, all.y=T)
  }
  cat(i, "Roung, dim is ", dim(MuFreq.group), "\n");
  flush.console();
}
###somehow it is possible that we have a good CDR3, but not mutation frequency, why???
##don't know
#MuFreq.group<-MuFreq.group[!is.na(MuFreq.group[,2])&MuFreq.group[,2]!=0,] #<--------no zero 
#MuFreq.group<-MuFreq.group[!is.na(MuFreq.group[,2]),]
MuFreq.group.IgM<-MuFreq.group
MuFreq.group.IgM$isotype<-"IgM"
MuFreq.group.IgM$tissue<-"Bone Marrow"
MuFreq.group.IgM$sampleName<-as.character(MuFreq.group.IgM$sampleName)
MuFreq.group.IgM$treatment<-conditions[MuFreq.group.IgM$sampleName,"treatment"]
###IgG
MuFreq.group<-NULL
MuFreq<-MuFreq.BM.IgG
#vgenes<-NULL
for(i in 1:length(MuFreq))
{
  mf<-MuFreq[[i]]$MuFreq
  
  temp<-data.frame("sampleName"=MuFreq[[i]]$sampleName,"MuFreq"=mf)
  
  #build the data frame
  if(i==1)
  {
    MuFreq.group<-temp
  } else {
    MuFreq.group<-rbind(MuFreq.group, temp)#, all.x=T, all.y=T)
  }
  cat(i, "Roung, dim is ", dim(MuFreq.group), "\n");
  flush.console();
}
###somehow it is possible that we have a good CDR3, but not mutation frequency, why???
##don't know
#MuFreq.group<-MuFreq.group[!is.na(MuFreq.group[,2])&MuFreq.group[,2]!=0,] #<--------no zero 
#MuFreq.group<-MuFreq.group[!is.na(MuFreq.group[,2]),]
MuFreq.group.IgG<-MuFreq.group

MuFreq.group.IgG$isotype<-"IgG"
MuFreq.group.IgG$tissue<-"Bone Marrow"
MuFreq.group.IgG$sampleName<-as.character(MuFreq.group.IgG$sampleName)
MuFreq.group.IgG$treatment<-conditions[MuFreq.group.IgG$sampleName,"treatment"]



#sort the colnames in each output
#MuFreq.table<-MuFreq.table[, c(1,conditions$index.fileOrder+1)]

MuFreq.table<-rbind(MuFreq.group.IgG, MuFreq.group.IgM)
MuFreq.table.BM<-MuFreq.table

#############+++++++++++
### SP IgM  and IgG together                                    |
########################
MuFreq.group<-NULL
MuFreq<-MuFreq.SP.IgM
#vgenes<-NULL
for(i in 1:length(MuFreq))
{
  mf<-MuFreq[[i]]$MuFreq
  
  temp<-data.frame("sampleName"=MuFreq[[i]]$sampleName,"MuFreq"=mf)
  
  #build the data frame
  if(i==1)
  {
    MuFreq.group<-temp
  } else {
    MuFreq.group<-rbind(MuFreq.group, temp)#, all.x=T, all.y=T)
  }
  cat(i, "Roung, dim is ", dim(MuFreq.group), "\n");
  flush.console();
}
###somehow it is possible that we have a good CDR3, but not mutation frequency, why???
##don't know
#MuFreq.group<-MuFreq.group[!is.na(MuFreq.group[,2])&MuFreq.group[,2]!=0,] #<--------no zero 
#MuFreq.group<-MuFreq.group[!is.na(MuFreq.group[,2]),]
MuFreq.group.IgM<-MuFreq.group
MuFreq.group.IgM$isotype<-"IgM"
MuFreq.group.IgM$tissue<-"Spleen"
MuFreq.group.IgM$sampleName<-as.character(MuFreq.group.IgM$sampleName)
MuFreq.group.IgM$treatment<-conditions[MuFreq.group.IgM$sampleName,"treatment"]

###IgG
MuFreq.group<-NULL
MuFreq<-MuFreq.SP.IgG
#vgenes<-NULL
for(i in 1:length(MuFreq))
{
  mf<-MuFreq[[i]]$MuFreq
  
  temp<-data.frame("sampleName"=MuFreq[[i]]$sampleName,"MuFreq"=mf)
  
  #build the data frame
  if(i==1)
  {
    MuFreq.group<-temp
  } else {
    MuFreq.group<-rbind(MuFreq.group, temp)#, all.x=T, all.y=T)
  }
  cat(i, "Roung, dim is ", dim(MuFreq.group), "\n");
  flush.console();
}
###somehow it is possible that we have a good CDR3, but not mutation frequency, why???
##don't know
#MuFreq.group<-MuFreq.group[!is.na(MuFreq.group[,2])&MuFreq.group[,2]!=0,] #<--------no zero 
#MuFreq.group<-MuFreq.group[!is.na(MuFreq.group[,2]),]
MuFreq.group.IgG<-MuFreq.group

MuFreq.group.IgG$isotype<-"IgG"
MuFreq.group.IgG$tissue<-"Spleen"
MuFreq.group.IgG$sampleName<-as.character(MuFreq.group.IgG$sampleName)
MuFreq.group.IgG$treatment<-conditions[MuFreq.group.IgG$sampleName,"treatment"]


#sort the colnames in each output
MuFreq.table<-rbind(MuFreq.group.IgG, MuFreq.group.IgM)
MuFreq.table.SP<-MuFreq.table

#############################################################################3
#####now do three way  analysis
###############################

#first carry the 

#MuFreq.table.SP<-MuFreq.table
MuFreq.table<-rbind(MuFreq.table.SP, MuFreq.table.BM)
 
save( file=here(data.figure3.dir,"muFreq_table_all.RData"),MuFreq.table)#muFreq.table, contains all the data includ

MuFreq.table<-MuFreq.table[!is.na(MuFreq.table[,2]),]#&MuFreq.table[,2]!=0,]    <---------------YES, with zero.
#MuFreq.table<-MuFreq.table[!is.na(MuFreq.table[,2])&MuFreq.table[,2]!=0,]

mufreq.zero.ratio<-aggregate(MuFreq.table$MuFreq, by=list(tissue=MuFreq.table$tissue, 
                                            treatment=MuFreq.table$treatment, isotype=MuFreq.table$isotype
                                                     , sampleName=MuFreq.table$sampleName), FUN=function(x){sum(x>0)/length(x)})
temp<-aggregate(MuFreq.table$MuFreq, by=list(tissue=MuFreq.table$tissue, 
                                            treatment=MuFreq.table$treatment, isotype=MuFreq.table$isotype
                                                     , sampleName=MuFreq.table$sampleName), FUN=function(x){length(x)})                                                     
#now get mean summarized data
mu.means<-aggregate(MuFreq.table$MuFreq, by=list(tissue=MuFreq.table$tissue, 
                                            treatment=MuFreq.table$treatment, isotype=MuFreq.table$isotype
                                                     , sampleName=MuFreq.table$sampleName), mean)
mu.means$treatment<-factor(mu.means$treatment, levels=c("PBS","OVA","OVA+PorB", "OVA+CpG", "OVA+Alum"))
mu.means$isotype<-factor(mu.means$isotype, levels=c("IgM", "IgG"))
mu.means$tissue<-factor(mu.means$tissue, levels=c("Spleen","Bone Marrow"))

mu.means.np<-mu.means[mu.means$treatment!="OVA+PorB",]
#mu.means.np<-mu.means.np[-20,]
#mu.means.np$trans<-log(mu.means.np$x/(1-mu.means.np$x))
#now we need to do testing
options(contrasts = c("contr.sum", "contr.poly"))
model_ls<-lm(x~tissue*treatment*isotype, data=mu.means.np)

bcx<-boxcox(model_ls)
lamda<-bcx$x[which.max(bcx$y)]
mu.means.np$trs<-(mu.means.np$x^lamda-1)/lamda

#testing for equal variance
bartlett.test(x~tissue, data=mu.means.np)
bartlett.test(x~isotype, data=mu.means.np)
bartlett.test(x~treatment, data=mu.means.np)

model_ls<-lm(x~tissue*treatment*isotype, data=mu.means.np)
#normality
shapiro.test(model_ls$residual)
qqnorm(model_ls$residual)
qqline(model_ls$residual)


Anova(model_ls, type=2)


#model_ls<-lm(x~tissue*treatment*isotype, data=mu.means.np)

model_ls.t<-lm(x~tissue*treatment*isotype, data=mu.means.np)


#em_tr<-emmeans(mod, ~tissue|treatment*isotype)
em_tr<-emmeans(model_ls, ~treatment|isotype*tissue)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
    "OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), "OVA - PBS"=c(-1, 1,0,0)
  #  "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#    "Bone Marrow - Spleen" = c(1,-1)
  )
  

ec1<-emmeans::contrast(em_tr, Set1, adjust='FDR')

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
    "OVA/Alum - OVA" = c(0,1,0,-1), "OVA/CpG - OVA"=c(0,-1,1,0), "OVA - PBS"=c(-1, 1,0,0)
  #  "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#    "Bone Marrow - Spleen" = c(1,-1)
  )
  

ec1.OVA<-emmeans::contrast(em_tr, Set1, adjust='FDR')


em_tr<-emmeans(model_ls, ~isotype|treatment)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
   # "OVA/Alum - PBS" = c(-1,0,0,0,1), "OVA/CpG - PBS"=c(-1,0,0,1,0), "OVA - PBS"=c(-1, 1,0,0,0)
    "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#    "Bone Marrow - Spleen" = c(1,-1)
  )

ec2<-emmeans::contrast(em_tr, Set1, adjust='FDR')

em_tr<-emmeans(model_ls, ~isotype|treatment*tissue)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
   # "OVA/Alum - PBS" = c(-1,0,0,0,1), "OVA/CpG - PBS"=c(-1,0,0,1,0), "OVA - PBS"=c(-1, 1,0,0,0)
    "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#    "Bone Marrow - Spleen" = c(1,-1)
  )

ec2iso<-emmeans::contrast(em_tr, Set1, adjust='FDR')

em_tr<-emmeans(model_ls, ~tissue|isotype*treatment)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
   # "OVA/Alum - PBS" = c(-1,0,0,0,1), "OVA/CpG - PBS"=c(-1,0,0,1,0), "OVA - PBS"=c(-1, 1,0,0,0)
   # "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
    "Bone Marrow - Spleen" = c(1,-1)
  )

ec3<-emmeans::contrast(em_tr, Set1, adjust='FDR')

####three way interaction for treatment
model_ls3<-lm(x~tissue*treatment*isotype, data=mu.means.np)
#em_tr<-emmeans(mod, ~tissue|treatment*isotype)
em_tr<-emmeans(model_ls3, ~treatment|isotype*tissue)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
    "OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), "OVA - PBS"=c(-1, 1,0,0)
  #  "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#    "Bone Marrow - Spleen" = c(1,-1)
  )

ect3<-emmeans::contrast(em_tr, Set1, adjust='FDR')


#do summarizing
m.ti<-aggregate(mu.means.np$x, by=list(tissue=mu.means.np$tissue), mean)
v.ti<-aggregate(mu.means.np$x, by=list(tissue=mu.means.np$tissue), FUN=function(x){sd(x)/sqrt(length(x))})
mv<-data.frame(mean=m.ti$x, var=v.ti$x, tissue=m.ti$tissue)
m1<-ggplot(data=mv, aes(x=tissue, y=mean, fill=tissue))+
    geom_bar(stat="identity")+theme_bw(base_size=15)+ylab("Mu Freq")+
    theme(legend.position="none", axis.title.x=element_blank())+#ylim(c(0,0.02))+
    geom_point(data=mu.means.np, aes(y=x, x=tissue),colour="grey", position=position_dodge(width = .9),size=5)+
    geom_errorbar(aes(ymin=mean-var, ymax=mean+var), width=.2,
                 position=position_dodge(.9))+
    geom_segment(x=1,xend=2, y=0.0125, yend=0.0125)+annotate(geom="text",x=1.5,y=0.013, label="p<0.0001")

   #save the ggplot2 figure for the future.
   save(mv, m1, mu.means.np, 
      file=here(data.figure3.dir,"mutFreq_tissueOnly_fig.RData"))

m.ti<-aggregate(mu.means.np$x, by=list(treatment=mu.means.np$treatment, isotype=mu.means.np$isotype), mean)
v.ti<-aggregate(mu.means.np$x, by=list(treatment=mu.means.np$treatment, isotype=mu.means.np$isotype), FUN=function(x){sd(x)/sqrt(length(x))})
mv<-data.frame(mean=m.ti$x, var=v.ti$x, treatment=m.ti$treatment, isotype=m.ti$isotype)   
m2<-ggplot(data=mv, aes(x=treatment, y=mean, fill=isotype))+
    geom_bar(stat="identity", position = position_dodge(width = 0.9))+theme_bw(base_size=15)+ylab("Mu Freq")+
    ylim(c(0,0.0185))+
    theme(legend.position=c(0.08,0.93), axis.title.x=element_blank(), legend.text=element_text(size=8),
           legend.title = element_blank() )+
     geom_point(data=mu.means.np, aes(y=x), colour="grey",position=position_dodge(width = .9),size=5)+
    geom_errorbar(aes(ymin=mean-var, ymax=mean+var), width=.2,
                 position=position_dodge(.9))+
     geom_segment(x=1.75,xend=2.25, y=0.0165, yend=0.0165)+annotate(geom="text",x=2,y=0.017, label="p<0.0001")+
     geom_segment(x=2.75,xend=3.25, y=0.0175, yend=0.0175)+ annotate(geom="text",x=3,y=0.018, label="p<0.0001")+
     geom_segment(x=3.75,xend=4.25, y=0.0170, yend=0.0170)+ annotate(geom="text",x=4,y=0.0175, label="p<0.0001")

m3<-ggplot(data=mv, aes(x=isotype, y=mean, fill=treatment))+
    geom_bar(stat="identity", position = position_dodge(width = 0.9))+theme_bw(base_size=15)+ylab("Mu Freq")+
    ylim(c(0,0.021))+
    theme(legend.position="right", axis.title.x=element_blank())+
    geom_point(data=mu.means.np, aes(y=x), colour="grey",position=position_dodge(width = .9),size=5)+
    geom_errorbar(aes(ymin=mean-var, ymax=mean+var), width=.2,
                 position=position_dodge(.9))    +
     geom_segment(x=0.65,xend=0.9, y=0.017, yend=0.017)+annotate(geom="text",x=0.775,y=0.0175, label="p<0.05")+
     geom_segment(x=0.65,xend=1.1, y=0.0185, yend=0.0185)+ annotate(geom="text",x=0.9,y=0.019, label="p<0.05")+
     geom_segment(x=0.65,xend=1.3, y=0.020, yend=0.020)+ annotate(geom="text",x=1,y=0.0205, label="p<0.05")+
     geom_segment(x=1.65,xend=1.9, y=0.0145, yend=0.0145)+annotate(geom="text",x=1.775,y=0.015, label="p<0.05")+
     geom_segment(x=1.65,xend=2.1, y=0.016, yend=0.016) +annotate(geom="text",x=1.9,y=0.0165, label="p<0.05")+
     geom_segment(x=1.65,xend=2.3, y=0.0175, yend=0.0175)    +annotate(geom="text",x=2,y=0.018, label="p<0.05")

       save(mv, m3, m2,mu.means.np, 
          file=here(data.figure3.dir,"mutFreq_TrByIso_fig.RData"))
tiff(file=here(output.dir,"treatment_pooledTissue.tiff"), 
      width=700, height=650)
m3
dev.off()     
#####do three way interaction
#-----------

m.3way<-aggregate(mu.means.np$x, by=list(treatment=mu.means.np$treatment, isotype=mu.means.np$isotype, tissue=mu.means.np$tissue), mean)
v.3way<-aggregate(mu.means.np$x, by=list(treatment=mu.means.np$treatment, isotype=mu.means.np$isotype, tissue=mu.means.np$tissue), FUN=function(x){sd(x)/sqrt(length(x))})
mv<-data.frame(mean=m.3way$x, var=v.3way$x, treatment=m.3way$treatment, isotype=m.3way$isotype, tissue=m.3way$tissue)   
    
    ######Note: here we are trying to plotting using the estimated pool equal variance!!!
    #we revert back to use the true sd for plotting.<--
#mv<-data.frame(em_tr)
#colnames(mv)[c(4,5)]<-c("mean", "var")
m3way<-ggplot(data=mv, aes(x=isotype, y=mean, fill=treatment))+
    geom_bar(stat="identity", position = position_dodge(width = 0.9))+theme_bw(base_size=15)+ylab("Mu Freq")+
    ylim(c(0,0.0135))+scale_fill_discrete(name=element_blank())+
    facet_grid(.~tissue)+
    theme(legend.position=c(0.15,0.8), axis.title.x=element_blank())+
    #theme(legend.position="right", axis.title.x=element_blank())+
     #geom_point(data=mu.means.np, aes(y=x), colour="grey",position=position_dodge(width = .9),size=5)+
    geom_errorbar(aes(ymin=mean-var, ymax=mean+var), width=.2,
                 position=position_dodge(.9))+
     geom_segment(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                            x=1.65,xend=2.1, y=0.012, yend=0.012)+
       geom_text(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                        x=1.9,y=0.0125, label="p=0.04")+
     geom_segment(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                        x=1.65,xend=2.3, y=0.013, yend=0.013)+
        geom_text(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                        x=2.0,y=0.0135, label="p=0.07")+

     geom_segment(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgM"),
                                            x=0.65,xend=0.9, y=0.0100, yend=0.0100)+
         geom_text(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgM"),
                                        x=0.8,y=0.0105, label="p=0.01")+
     geom_segment(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                            x=0.65,xend=1.1, y=0.011, yend=0.011)+
         geom_text(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                        x=0.9,y=0.0115, label="p<0.01")+
         geom_segment(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                            x=0.65,xend=1.3, y=0.012, yend=0.012)+
         geom_text(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                        x=1,y=0.0125, label="p<0.05")+
    geom_segment(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgG"),
                                            x=1.65,xend=1.9, y=0.0085, yend=0.0085)+
         geom_text(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgG"),
                                        x=1.8,y=0.009, label="p<0.05")+
     geom_segment(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgG"),
                                            x=1.85,xend=2.325, y=0.0095, yend=0.0095)+
         geom_text(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgG"),
                                        x=2.05,y=0.01, label="p<0.05")#+
tiff(file=here(output.dir,"muFreq3way_treatment.tiff"), 
    width=800, height=500)
m3way
dev.off();

save(mv, m3way,mu.means.np, 
    file=here(data.figure3.dir,"m3way_treatment.RData" ))

tiff(file=here(output.dir,"mufreq.tiff"),width=850, height=650)       
#   library(ggpubr)     
        ggarrange(
                        ggarrange(m1,m2,
                                nrow=1, ncol=2, widths=c(1,2.5), labels=c("A","B")
                        ), m3way,
                        nrow=2, ncol=1, labels=c("","C")
        
        )
        dev.off()

tiff(file=here(output.dir,"mufreq_treatReduce.tiff"),
      width=850, height=650)       
#   library(ggpubr)     
        ggarrange(
                        ggarrange(m1,m2,
                                nrow=1, ncol=2, widths=c(1,2.5), labels=c("A","B")
                        ), m3,
                        nrow=2, ncol=1, labels=c("","C")
        
        )
dev.off()

########    isotype  3 way
m.3way<-aggregate(mu.means.np$x, by=list(treatment=mu.means.np$treatment, isotype=mu.means.np$isotype, tissue=mu.means.np$tissue), mean)
v.3way<-aggregate(mu.means.np$x, by=list(treatment=mu.means.np$treatment, isotype=mu.means.np$isotype, tissue=mu.means.np$tissue), FUN=function(x){sd(x)/sqrt(length(x))})
mv<-data.frame(mean=m.3way$x, var=v.3way$x, treatment=m.3way$treatment, isotype=m.3way$isotype, tissue=m.3way$tissue)   

mv.split<-mv
#mv.split<-data.frame(em_tr)
# colnames(mv.split)[c(4,5)]<-c("mean","var")
 mv.split$section<-"Immunized"
 mv.split[mv.split$treatment=="PBS","section"]<-"PBS"

m3wayIso<-ggplot(data=mv.split, aes(x=treatment, y=mean, fill=isotype))+
    geom_bar(stat="identity", position = position_dodge(width = 0.9))+theme_bw(base_size=15)+ylab("Mu Freq")+
    ylim(c(0,0.0135))+
    facet_grid(.~tissue)+
    theme(legend.position="right", axis.title.x=element_blank())+
   #  geom_dotplot(data=mu.means.np, aes(y=x, x=treatment), colour="black",
   #                         binaxis='y', stackdir='center', position="dodge" ,binwidth=0.0004)+
    geom_errorbar(aes(ymin=mean-var, ymax=mean+var), width=.2,
                 position=position_dodge(.9))+
     geom_segment(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                            x=1.75,xend=2.25, y=0.0115, yend=0.0115)+
       geom_text(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                        x=2.05,y=0.012, label="p<0.001")+
     geom_segment(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                        x=2.75,xend=3.25, y=0.012, yend=0.012)+
        geom_text(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                        x=3.05,y=0.0125, label="p<0.001")+
     geom_segment(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                            x=3.75,xend=4.25, y=0.0125, yend=0.0125)+
         geom_text(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                        x=4.05,y=0.013, label="p<0.001")+
     geom_segment(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgG"),
                                            x=1.75,xend=2.25, y=0.0085, yend=0.0085)+
         geom_text(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgG"),
                                        x=2.05,y=0.009, label="p<0.001")+
    geom_segment(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgG"),
                                            x=2.75,xend=3.25, y=0.009, yend=0.009)+
         geom_text(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgG"),
                                        x=3.05,y=0.0095, label="p<0.001")+
     geom_segment(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgG"),
                                            x=3.75,xend=4.25, y=0.0075, yend=0.0075)+
         geom_text(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgG"),
                                        x=4.05,y=0.008, label="p<0.001")#+
 
 m3wayIso.pbs<-ggplot(data=mv.split[mv.split$section=="PBS",], aes(x=tissue, y=mean, fill=isotype))+
    geom_bar(stat="identity", position = position_dodge(width = 0.9))+
            theme_bw(base_size=15)+ylab("Mu Freq")+
    ylim(c(0,0.0135))+
    facet_grid(.~section)+
    theme(legend.position=c(0.3,0.8), axis.title.x=element_blank())+
    geom_errorbar(aes(ymin=mean-var, ymax=mean+var), width=.2,
                 position=position_dodge(.9))
                 
  m3wayIso.imm<-ggplot(data=mv.split[mv.split$section=="Immunized",], aes(x=treatment, y=mean, fill=isotype))+
    geom_bar(stat="identity", position = position_dodge(width = 0.9))+
            theme_bw(base_size=15)+ylab("Mu Freq")+
    ylim(c(0,0.0135))+
    facet_grid(.~tissue)+
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank()
                , axis.text.y=element_blank())+
    geom_errorbar(aes(ymin=mean-var, ymax=mean+var), width=.2,
                 position=position_dodge(.9))+
                 geom_segment(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                            x=0.75,xend=1.25, y=0.0115, yend=0.0115)+
       geom_text(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                        x=1.05,y=0.012, label="p<0.001")+
     geom_segment(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                        x=1.75,xend=2.25, y=0.012, yend=0.012)+
        geom_text(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                        x=2.05,y=0.0125, label="p<0.001")+
     geom_segment(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                            x=2.75,xend=3.25, y=0.0125, yend=0.0125)+
         geom_text(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                        x=3.05,y=0.013, label="p<0.001")+
     geom_segment(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgG"),
                                            x=0.75,xend=1.25, y=0.0085, yend=0.0085)+
         geom_text(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgG"),
                                        x=1.05,y=0.009, label="p<0.005")+
    geom_segment(data=data.frame(tissue="Spleen", treatment="OVA+CpG", isotype="IgG"),
                                            x=1.75,xend=2.25, y=0.009, yend=0.009)+
         geom_text(data=data.frame(tissue="Spleen", treatment="OVA+CpG", isotype="IgG"),
                                        x=2.05,y=0.0095, label="p<0.005")+
     geom_segment(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgG"),
                                            x=2.75,xend=3.25, y=0.0075, yend=0.0075)+
         geom_text(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgG"),
                                        x=3.05,y=0.008, label="p<0.01")
 
 ggarrange(
                        m3wayIso.pbs, m3wayIso.imm,
                        nrow=1, ncol=2,widths=c(1,1.8)
        
        )     
save(m3wayIso.pbs, m3wayIso.imm, mv.split, 
    file=here(data.figure3.dir,"m3wayIso_split.RData"))        

tiff(file=here(output.dir,"muFreq3way_isotype.tiff"), width=800, height=500)
m3wayIso
dev.off();

########################tissue
mv.split.tis<-mv
#mv.split.tis<-data.frame(em_tr)
# colnames(mv.split.tis)[c(4,5)]<-c("mean","var")
 mv.split.tis$section<-"Immunized"
 mv.split.tis[mv.split.tis$treatment=="PBS","section"]<-"PBS"
 m3wayTis.pbs<-ggplot(data=mv.split.tis[mv.split.tis$section=="PBS",], aes(x=isotype, y=mean, fill=tissue))+
    geom_bar(stat="identity", position = position_dodge(width = 0.9))+
            theme_bw(base_size=15)+ylab("Mu Freq")+
    ylim(c(0,0.0135))+
    facet_grid(.~section)+
    theme(legend.position=c(0.3,0.88), axis.title.x=element_blank())+
    geom_errorbar(aes(ymin=mean-var, ymax=mean+var), width=.2,
                 position=position_dodge(.9))+
    geom_segment(data=data.frame(tissue="Spleen", treatment="PBS", isotype="IgM"),
                                            x=1.75,xend=2.25, y=0.0095, yend=0.0095)+
         geom_text(data=data.frame(tissue="Spleen", treatment="PBS", isotype="IgM"),
                                        x=2.0,y=0.010, label="p<0.005")+
      geom_segment(data=data.frame(tissue="Spleen", treatment="PBS", isotype="IgM"),
                                            x=0.75,xend=1.25, y=0.0085, yend=0.0085)+
         geom_text(data=data.frame(tissue="Spleen", treatment="PBS", isotype="IgM"),
                                        x=1.0,y=0.009, label="p<0.005")
                 

m3wayTis<-ggplot(data=mv.split.tis, aes(x=treatment, y=mean, fill=tissue))+
    geom_bar(stat="identity", position = position_dodge(width = 0.9))+theme_bw(base_size=15)+ylab("Mu Freq")+
    ylim(c(0,0.013))+
    facet_grid(.~isotype)+
    theme(legend.position="right", axis.title.x=element_blank())+
     geom_dotplot(data=mu.means.np, aes(y=x, x=treatment), colour="black",
                            binaxis='y', stackdir='center', position="dodge" ,binwidth=0.0003)+
    geom_errorbar(aes(ymin=mean-var, ymax=mean+var), width=.2,
                 position=position_dodge(.9))+
     geom_segment(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                            x=0.75,xend=1.25, y=0.01, yend=0.010)+
       geom_text(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                        x=1.05,y=0.0105, label="p=0.02")+
     geom_segment(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                        x=2.75,xend=3.25, y=0.0115, yend=0.0115)+
        geom_text(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                        x=3.05,y=0.012, label="p=0.03")+
     geom_segment(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                            x=3.75,xend=4.25, y=0.012, yend=0.012)+
         geom_text(data=data.frame(tissue="Bone Marrow", treatment="OVA", isotype="IgG"),
                                        x=4.05,y=0.0125, label="p<0.01")+
                                        
     geom_segment(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgM"),
                                            x=0.75,xend=1.25, y=0.01, yend=0.01)+
         geom_text(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgM"),
                                        x=1.05,y=0.0105, label="p<0.015")+
    geom_segment(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgM"),
                                            x=2.75,xend=3.25, y=0.007, yend=0.007)+
         geom_text(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgM"),
                                        x=3.05,y=0.0075, label="p=0.08")+
     geom_segment(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgM"),
                                            x=3.75,xend=4.25, y=0.008, yend=0.008)+
         geom_text(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgM"),
                                        x=4.05,y=0.0085, label="p<0.005")#+
#     geom_segment(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgG"),
 #                                           x=1.65,xend=2.3, y=0.0115, yend=0.0115)+
 #        geom_text(data=data.frame(tissue="Spleen", treatment="OVA", isotype="IgG"),
 #                                       x=2.0,y=0.0119, label="p<0.08")


m3wayTis2<-ggplot(data=mv.split.tis, aes(x=isotype, y=mean, fill=tissue))+
    geom_bar(stat="identity", position = position_dodge(width = 0.9))+theme_bw(base_size=15)+ylab("Mu Freq")+
    ylim(c(0,0.013))+
    facet_grid(.~treatment)+
    theme(legend.position=c(0.35,0.88), axis.title.x=element_blank())+
    # geom_dotplot(data=mu.means.np, aes(y=x, x=isotype), colour="black",
    #                        binaxis='y', stackdir='center', position="dodge" ,binwidth=0.0003)+
    geom_errorbar(aes(ymin=mean-var, ymax=mean+var), width=.2,
                 position=position_dodge(.9))+
      geom_segment(data=data.frame(tissue="Spleen", treatment="PBS", isotype="IgM"),
                                            x=1.75,xend=2.25, y=0.0095, yend=0.0095)+
         geom_text(data=data.frame(tissue="Spleen", treatment="PBS", isotype="IgM"),
                                        x=2.0,y=0.010, label="p<0.005")+
      geom_segment(data=data.frame(tissue="Spleen", treatment="PBS", isotype="IgM"),
                                            x=0.75,xend=1.25, y=0.0085, yend=0.0085)+
         geom_text(data=data.frame(tissue="Spleen", treatment="PBS", isotype="IgM"),
                                        x=1.0,y=0.009, label="p<0.005")+
          geom_segment(data=data.frame(tissue="Spleen", treatment="OVA+CpG", isotype="IgM"),
                                            x=1.75,xend=2.25, y=0.011, yend=0.011)+
         geom_text(data=data.frame(tissue="Spleen", treatment="OVA+CpG", isotype="IgM"),
                                        x=2.0,y=0.0115, label="p<0.005")+
      geom_segment(data=data.frame(tissue="Spleen", treatment="OVA+Alum", isotype="IgM"),
                                            x=1.75,xend=2.25, y=0.012, yend=0.012)+
         geom_text(data=data.frame(tissue="Spleen", treatment="OVA+Alum", isotype="IgG"),
                                        x=2.0,y=0.0125, label="p<0.0005")

tiff(file=here(output.dir,"muFreq3way_tissue.tiff"), width=800, height=500)
m3wayTis2
dev.off();                                        

#now don't save the data frame, since it is save in treatment plot
save(m3wayTis2, m3wayIso, m3wayTis.pbs, mv.split.tis,
    file=here(data.figure3.dir,"m3way_tisIso.RData"))

model_ls<-lm(x~tissue*treatment*isotype, data=mu.means.np)
#model_ls<-lm(x~tissue*treatment, data=mu.means.np[mu.means.np$isotype=="IgM",])
#em_tr<-emmeans(mod, ~tissue|treatment*isotype)
#em_tr<-emmeans(model_ls, ~treatment|tissue)
em_tr<-emmeans(model_ls, ~treatment|isotype*tissue)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
    "OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), "OVA - PBS"=c(-1, 1,0,0)
  #  "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#    "Bone Marrow - Spleen" = c(1,-1)
  )

ec1<-emmeans::contrast(em_tr, Set1, adjust='FDR')

em_tr<-emmeans(model_ls, ~isotype|treatment*tissue)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
   # "OVA/Alum - PBS" = c(-1,0,0,0,1), "OVA/CpG - PBS"=c(-1,0,0,1,0), "OVA - PBS"=c(-1, 1,0,0,0)
    "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#    "Bone Marrow - Spleen" = c(1,-1)
  )

ec2<-emmeans::contrast(em_tr, Set1, adjust='FDR')


em_tr<-emmeans(model_ls, ~tissue|isotype*treatment)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
   # "OVA/Alum - PBS" = c(-1,0,0,0,1), "OVA/CpG - PBS"=c(-1,0,0,1,0), "OVA - PBS"=c(-1, 1,0,0,0)
   # "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
    "Bone Marrow - Spleen" = c(1,-1)
  )

ec3<-emmeans::contrast(em_tr, Set1, adjust='FDR')


#difference of difference
model_ls3<-lm(x~tissue*treatment*isotype, data=mu.means.np)
#model_ls3<-lm(x~tissue*treatment, data=mu.means.np[mu.means.np$isotype=="IgG",])
#now doing constrast for doing comparison of difference of BM spleen
#em_tr<-emmeans(model_ls3, ~treatment*tissue|isotype)
em_tr<-emmeans(model_ls3, ~treatment*tissue)
#em_tr %>% test(joint=F)

Set1 <- list(  #order PBS BM, OVA BM, CpG BM, Alum BM, PBS spleen, OVA spleen, CpG spleen and Alum BM
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
    "BM-Spleen | OVA vs PBS" = c(-1,1,0,0, 1,-1,0,0), "BM-Spleen | CpG vs OVA" = c(0,-1,1,0, 0,1,-1,0),
    "BM-Spleen | Alum vs. OVA" = c(0,-1,0,1, 0,1,0,-1),"BM-Spleen | Alum vs PBS" = c(-1,0,0,1,1,0,0,-1)
    ,"BM-Spleen | Alum vs CpG" = c(0,0,-1,1,0,0,1,-1)
  #  "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#    "Bone Marrow - Spleen" = c(1,-1)
  )
  df.ec<-emmeans::contrast(em_tr, Set1, adjust='FDR')
  
  Set1 <- list(  #order PBS BM, OVA BM, CpG BM, Alum BM, PBS spleen, OVA spleen, CpG spleen and Alum BM
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
    #"BM-Spleen | OVA vs PBS" = c(-1,1,0,0, 1,-1,0,0), 
    "BM-Spleen | CpG vs OVA" = c(0,-1,1,0, 0,1,-1,0),
    "BM-Spleen | Alum vs. OVA" = c(0,-1,0,1, 0,1,0,-1),#"BM-Spleen | Alum vs PBS" = c(-1,0,0,1,1,0,0,-1)
    "BM-Spleen | Alum vs CpG" = c(0,0,-1,1,0,0,1,-1)
  #  "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#    "Bone Marrow - Spleen" = c(1,-1)
  )
  df.ec2<-emmeans::contrast(em_tr, Set1, adjust='FDR')
#start plotting. the 


#now we do

em_tr<-emmeans(model_ls, ~tissue|isotype*treatment)
#em_tr<-emmeans(model_ls, ~tissue|treatment)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
   # "OVA/Alum - PBS" = c(-1,0,0,0,1), "OVA/CpG - PBS"=c(-1,0,0,1,0), "OVA - PBS"=c(-1, 1,0,0,0)
   # "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
    "Bone Marrow - Spleen" = c(-1,1)
  )

ec.df2<-emmeans::contrast(em_tr, Set1, adjust='FDR')
x<-plot(ec.df2, horizontal=F, by=NULL, plotit=F)
#x$isotype<-factor(x$isotype, levels=c("IgM", "IgG"))
w=0.1
x$min_y<-x$the.emmean-x$SE*0.6
x$max_y<-x$the.emmean+x$SE*0.6
md.line<-ggplot(data=x, aes(x=treatment,y=the.emmean, colour=isotype, group=isotype  ))+
                                        ylab("Mu Freq Difference\n (Bone Marrow - Spleen)")+geom_point(aes(), size=5)+
                                       geom_line(size=0.25, linetype="dotted")+scale_x_discrete(expand=expansion(-0.0,0.0))+
                                       theme (text=element_text(size=16),axis.title.x = element_blank())+
                                      geom_errorbar(aes(ymin=min_y, ymax=max_y), width=.0)+
                                       geom_segment(aes(x=0.8, xend=4.2, y=0.0, yend=0),size=0.8, colour="grey", linetype=2)+
                                       geom_text(#data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0.0005, treatment="OVA", isotype="IgG"),
                                        aes(x=1.1, y=0.0002,
                                            label=c("Bone Marrow")), hjust=0,size=5,colour="red")+
                                     geom_text(#data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0.0005, treatment="OVA", isotype="IgG"),
                                        aes(x=1.1, y=-0.0002,
                                            label=c("Spleen")), hjust=0,size=5, colour="red")+
                                    geom_segment(aes(x=1.0, xend=1.0, y=0.0001, yend=0.0005),
                                                                                    arrow = arrow(length = unit(0.2, "cm"),type="closed"), size=1, colour="red", linetype=1)+
                                    geom_segment(aes(x=1.0, xend=1.0, y=-0.0001, yend=-0.0005),
                                                                                    arrow = arrow(length = unit(0.2, "cm"),type="closed"), size=1, colour="red", linetype=1)#+
                                     #geom_segment(#data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0, treatment="OVA", isotype=x[1,2]),
                                    #            aes(x=1, xend=2, y=0.0043, yend=0.0043),
                                    #                    size=1,  linetype=1)+
                                    # geom_segment(data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0, treatment="OVA", isotype="IgM"),
                                    #            aes(x=2, xend=4, y=0.005, yend=0.005),
                                    #                    size=1,  linetype=1)+
                                   # geom_text(data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0, treatment="OVA", isotype=c("IgM","IgM")),
                                    #    aes(y=c(0.0045, 0.0052),x=c(1.5,3), 
                                     #       label=c("p<0.05","p<0.05")),size=6)+theme(legend.position = c(0.8, 0.3))

png(file=here(output.dir,"tissueDiffLine.png"), width=800, height=500)
  md.line
dev.off(); 

save(x, w, md.line, 
    file=here(data.figure3.dir,"mutationDifference_fig.RData" ))                                           

em_tr<-emmeans(model_ls, ~tissue|treatment)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
   # "OVA/Alum - PBS" = c(-1,0,0,0,1), "OVA/CpG - PBS"=c(-1,0,0,1,0), "OVA - PBS"=c(-1, 1,0,0,0)
   # "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
    "Bone Marrow - Spleen" = c(1,-1)
  )

ec.df2<-emmeans::contrast(em_tr, Set1, adjust='FDR')
x<-plot(ec.df2, horizontal=F, by=NULL, plotit=F)
#x$isotype<-factor(x$isotype, levels=c("IgM", "IgG"))
w=0.1
                                            
md<-ggplot(data=x, aes(x=treatment,y=the.emmean ))+
                                        ylab("Mu Freq Difference\n (Bone Marrow - Spleen)")+geom_point(aes(), size=3, colour="red")+
                                       # xlab("")+
                                   geom_rect(aes(  xmin=as.integer(x$treatment)-w,xmax=as.integer(x$treatment)+w,
                                        ymin=lower.CL, ymax=upper.CL ),colour="dark green", size=0.1, fill="light blue")+
                                   theme_bw(base_size=12)+theme(axis.title.x=element_blank())+
                                   geom_point(aes(), size=3, colour="red")+
                                   geom_segment(aes(x=0, xend=5, y=0.0, yend=0),size=0.8, colour="grey", linetype=2)+
                                   geom_segment(data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0, treatment="OVA", isotype=x[1,2]),
                                                aes(x=0, xend=5, y=x[1,3], yend=x[1,3]),size=1, colour="red", linetype=2)+
                                    geom_segment(#data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0, treatment="OVA", isotype=x[1,2]),
                                                aes(x=1, xend=2, y=x$upper.CL[1]+0.0005, yend=x$upper.CL[1]+0.0005),
                                                        size=1, colour="black", linetype=1)+
                                    geom_segment(#data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0, treatment="OVA", isotype=x[1,2]),
                                                aes(x=2, xend=4, y=x$upper.CL[4]+0.0005, yend=x$upper.CL[4]+0.0005),
                                                        size=1, colour="black", linetype=1)+
                                   geom_segment(aes(x=0.3, xend=0.3, y=0.0009, yend=0.0018),
                                                                                    arrow = arrow(length = unit(0.2, "cm"),type="closed"), size=1, colour="red", linetype=1)+
                                    geom_segment(aes(x=0.3, xend=0.3, y=-0.0009, yend=-0.0018),
                                                                                    arrow = arrow(length = unit(0.2, "cm"),type="closed"), size=1, colour="red", linetype=1)+                                                
                                   geom_text(#data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0, treatment="OVA", isotype="IgG"),
                                        aes(y=upper.CL+0.00075,x=c(1.5, 1,2,3), 
                                            label=c("p=0.05",  "","", "p<0.01")))+
                                    geom_text(#data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0.0005, treatment="OVA", isotype="IgG"),
                                        aes(x=0.45, y=0.0005,
                                            label=c("Bone Marrow")),size=4)+
                                     geom_text(#data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0.0005, treatment="OVA", isotype="IgG"),
                                        aes(x=0.45, y=-0.0005,
                                            label=c("Spleen")),size=4)
png(file=here(output.dir,"tissueDiffBar.png"), width=800, height=500)
md
dev.off(); 


x.2<-plot(ec.df2, horizontal=F, by=NULL, plotit=F)
x.2<-x.2[x.2$treatment!="PBS",]
x.2$treatment<-factor(x.2$treatment, levels=c("OVA", "OVA+CpG", "OVA+Alum"))
#x$isotype<-factor(x$isotype, levels=c("IgM", "IgG"))
w=0.1
md2<-ggplot(data=x.2, aes(x=treatment,y=the.emmean ))+
                                        ylab("Mu Freq Difference\n (Bone Marrow - Spleen)")+geom_point(aes(), size=3, colour="red")+
                                       # xlab("")+
                                   geom_rect(aes(  xmin=as.integer(x.2$treatment)-w,xmax=as.integer(x.2$treatment)+w,
                                        ymin=lower.CL, ymax=upper.CL ),colour="dark green", size=0.1, fill="light blue")+
                                   theme_bw(base_size=12)+theme(axis.title.x=element_blank())+
                                   geom_point(aes(), size=3, colour="red")+
                                   geom_segment(aes(x=0, xend=3.4, y=0.0, yend=0),size=0.8, colour="grey", linetype=2)+
                                   #geom_segment(data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0, treatment="OVA", isotype=x[1,2]),
                                   #            aes(x=0, xend=5, y=x[1,3], yend=x[1,3]),size=1, colour="red", linetype=2)+
                                    geom_segment(#data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0, treatment="OVA", isotype=x[1,2]),
                                                aes(x=1, xend=2, y=x.2$upper.CL[2]+0.0005, yend=x.2$upper.CL[2]+0.0005),
                                                        size=1, colour="black", linetype=1)+
                                    geom_segment(#data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0, treatment="OVA", isotype=x[1,2]),
                                                aes(x=1, xend=3, y=x.2$upper.CL[3]+0.0005, yend=x.2$upper.CL[3]+0.0005),
                                                        size=1, colour="black", linetype=1)+
                                   geom_segment(aes(x=0.3, xend=0.3, y=0.0009, yend=0.0018),
                                                                                    arrow = arrow(length = unit(0.2, "cm"),type="closed"), size=1, colour="red", linetype=1)+
                                    geom_segment(aes(x=0.3, xend=0.3, y=-0.0009, yend=-0.0018),
                                                                                    arrow = arrow(length = unit(0.2, "cm"),type="closed"), size=1, colour="red", linetype=1)+                                                
                                   geom_text(#data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0, treatment="OVA", isotype="IgG"),
                                        aes(y=upper.CL+0.00075,x=c(1.5, 1.5,2), 
                                            label=c("","p=0.05",  "p=0.01")))+
                                    geom_text(#data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0.0005, treatment="OVA", isotype="IgG"),
                                        aes(x=0.45, y=0.0005,
                                            label=c("Bone Marrow")),size=4)+
                                     geom_text(#data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0.0005, treatment="OVA", isotype="IgG"),
                                        aes(x=0.45, y=-0.0005,
                                            label=c("Spleen")),size=4)#+
                                  # facet_grid(.~isotype)

save(x.2, w, md2, 
    file=here(data.figure3.dir,"mutationDifferenceWOPBS_fig.RData" ))

###now want to do tissue difference manually
# get individual differences.
mu.mean.np.mouse<-mu.means.np
mu.mean.np.mouse$mouseNum<-substr(mu.mean.np.mouse$sampleName, 3,10)
mvm<-aggregate(mu.mean.np.mouse$x, by=list(mu.mean.np.mouse$mouseNum, mu.mean.np.mouse$isotype), 
                        FUN=function(x){x[1]-x[2]})
mvm$treatment<-mu.mean.np.mouse[mu.mean.np.mouse$isotype=="IgG"&mu.mean.np.mouse$tissue=="Bone Marrow","treatment"]
colnames(mvm)<-c("mouse", "isotype", "muDiff", "treatment")
mvm<-mvm[mvm$mouse!=5,]
lm.mvm<-lm(data=mvm, muDiff~isotype*treatment)
#library(multcomp)
bcx<-boxcox(lm.mvm)
lamda<-bcx$x[which.max(bcx$y)]
mvm$trs<-(mvm$muDiff^lamda-1)/lamda
lm.mvmtrs<-lm(data=mvm, trs~isotype*treatment)
qqnorm(lm.mvmtrs$residual)
qqline(lm.mvmtrs$residual)
Anova(lm.mvmtrs)
#Response: trs
#                     Sum Sq Df F value Pr(>F)  
#isotype           0.0067768  1  4.7437 0.0470 *
#treatment         0.0039471  3  0.9210 0.4561  
#isotype:treatment 0.0030191  3  0.7044 0.5650  
#Residuals         0.0200003 14  

em_tr<-emmeans(lm.mvmtrs, ~treatment|isotype)
#em_tr %>% test(joint=F)

Set1 <- list(
    "OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), "OVA - PBS"=c(-1, 1,0,0)
    , "OVA/Alum - OVA" = c(0,-1,0,1)
   # "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
  #  "Bone Marrow - Spleen" = c(1,-1)
  )

ec.d<-emmeans::contrast(em_tr, Set1, adjust='none')

#########now start doing transformation before getting mean
#now let's do compositions
#library(compositions)
#
mu.means<-aggregate(MuFreq.table$MuFreq, by=list(tissue=MuFreq.table$tissue, 
                                            isotype=MuFreq.table$isotype
                                                     , sampleName=MuFreq.table$sampleName), FUN=function(x){
                                                                    xmin=min(x[x>0]);
                                                                   xmin*0.33
                                                     })
names(mu.means)[4]="MuDL"
mu.means<-merge(mu.means, MuFreq.table)
mu.means[mu.means$MuFreq==0,"MuFreq"]<-mu.means[mu.means$MuFreq==0,"MuDL"]                                                     
mu.means<-aggregate((MuFreq.table$MuFreq), by=list(tissue=MuFreq.table$tissue, 
                                            treatment=MuFreq.table$treatment, isotype=MuFreq.table$isotype
                                                     , sampleName=MuFreq.table$sampleName), mean)
#mu.means<-aggregate(log(mu.means$MuFreq/(1-mu.means$MuFreq)), by=list(tissue=MuFreq.table$tissue, 
#                                            treatment=MuFreq.table$treatment, isotype=MuFreq.table$isotype
#                                                     , sampleName=MuFreq.table$sampleName), median)
mu.means$treatment<-factor(mu.means$treatment, levels=c("PBS","OVA","OVA+PorB", "OVA+CpG", "OVA+Alum"))
mu.means$isotype<-factor(mu.means$isotype, levels=c("IgM", "IgG"))
mu.means$tissue<-factor(mu.means$tissue, levels=c("Spleen","Bone Marrow"))

mu.means.np<-mu.means[mu.means$treatment!="OVA+PorB",]

mvm$muDiff<-mvm$trs
mvm.mean<-aggregate(mvm$muDiff, by=list(mvm$isotype,mvm$treatment), mean)
mvm.se<-as.data.frame(em_tr)#aggregate(mvm$muDiff, by=list(mvm$isotype,mvm$treatment), FUN=function(x){sd(x)/sqrt(length(x))})
colnames(mvm.mean)<-c("isotype", "treatment", "x")
#colnames(mvm.se)<-c("isotype", "treatment", "x")
mvm.mean$se<-mvm.se$SE
ggplot(data=mvm.mean, aes(x=treatment,y=x,group=isotype, color=isotype))+
    geom_line()+geom_point()+
    geom_errorbar(aes(ymin=x-se, ymax=x+se), width=.2)#,
               #  position=position_dodge(.8))

#no we need to do the difference from the beginning for each individual sequences

