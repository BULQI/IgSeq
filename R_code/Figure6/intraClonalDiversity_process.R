#R code to process intraclonal diversity results
#
# Make sure you have read ./ReadMe.txt for instructions
# of running scripts and preparing data.
# 

#library
library(plyr)
library(ggplot2)
library(frLib)
library(rcPkg)
library(vegan)
library(compositions)
library(here)
library(ggalt)
library(ggpubr)
library(MASS)

 #setwd("/home/feng/Feng/LAB/Wetzler_PorinB/wl03R2/newPipeline")
 data.dir<-"Data/Figure6"

 load(here(data.dir,"clones.df.RData"))
 topn<-40
 #get the top40 by size the 
 #first IgG BM
 temp<-BM.df.clones[["IgG"]]
 snames<-unique(temp$sampleName)
 dto<-data.frame()
 for(i in 1:length(snames))
 {
    tp<-temp[temp$sampleName==snames[i],]
    tp<-tp[order(tp$X.Members, decreasing=T),]
    tp<-tp[1:topn,]
    dto<-rbind(dto, tp)
    dto<-dto[!is.na(dto$X.Members),]
 }
 meanMu.BM.IgG.top40<-dto
 
 temp<-BM.df.clones[["IgM"]]
 snames<-unique(temp$sampleName)
 dto<-data.frame()
 for(i in 1:length(snames))
 {
    tp<-temp[temp$sampleName==snames[i],]
    tp<-tp[order(tp$X.Members, decreasing=T),]
    tp<-tp[1:topn,]
    dto<-rbind(dto, tp)
    dto<-dto[!is.na(dto$X.Members),]
 }
 meanMu.BM.IgM.top40<-dto
 
 temp<-SP.df.clones[["IgG"]]
 snames<-unique(temp$sampleName)
 dto<-data.frame()
 for(i in 1:length(snames))
 {
    tp<-temp[temp$sampleName==snames[i],]
    tp<-tp[order(tp$X.Members, decreasing=T),]
    tp<-tp[1:topn,]
    dto<-rbind(dto, tp)
    dto<-dto[!is.na(dto$X.Members),]
 }
 meanMu.SP.IgG.top40<-dto
 
 temp<-SP.df.clones[["IgM"]]
 snames<-unique(temp$sampleName)
 dto<-data.frame()
 for(i in 1:length(snames))
 {
    tp<-temp[temp$sampleName==snames[i],]
    tp<-tp[order(tp$X.Members, decreasing=T),]
    tp<-tp[1:topn,]
    dto<-rbind(dto, tp)
    dto<-dto[!is.na(dto$X.Members),]
 }
 meanMu.SP.IgM.top40<-dto
 
 #
 load(here(data.dir,"intraClonalDiversity_BM.RData"))  #idis.BM.IgM, idis.BM.IgG
#generated in ./intraClonal_batch.R, 
 load(here(data.dir,"intraClonalDiversity_SP.RData"))  #idis.SP.IgM, idis.SP.IgG
#generated in ./intraClonal_batch.R  

#get the top 40 clones for testing

top40<-idis.BM.IgG
dto<-data.frame()
for( i in 1:length(unique(top40[,1])))
{
    sname<-unique(top40[,1])[i]
    temp<-top40[top40$sampleName==sname,]
    temp<-temp[order(temp$xMembers, decreasing=T),]
    temp<-temp[c(1:topn),]
    dto<-rbind(dto, temp)
}
idis.BM.IgG.top40<-dto
idis.BM.IgG.top40<-idis.BM.IgG.top40[!is.na(idis.BM.IgG.top40[,1]),]

top40<-idis.BM.IgM
dto<-data.frame()
for( i in 1:length(unique(top40[,1])))
{
    sname<-unique(top40[,1])[i]
    temp<-top40[top40$sampleName==sname,]
    temp<-temp[order(temp$xMembers, decreasing=T),]
    temp<-temp[c(1:topn),]
    dto<-rbind(dto, temp)
}
idis.BM.IgM.top40<-dto
idis.BM.IgM.top40<-idis.BM.IgM.top40[!is.na(idis.BM.IgM.top40[,1]),]

top40<-idis.SP.IgM
dto<-data.frame()
for( i in 1:length(unique(top40[,1])))
{
    sname<-unique(top40[,1])[i]
    temp<-top40[top40$sampleName==sname,]
    temp<-temp[order(temp$xMembers, decreasing=T),]
    temp<-temp[c(1:topn),]
    dto<-rbind(dto, temp)
}
idis.SP.IgM.top40<-dto

top40<-idis.SP.IgG
dto<-data.frame()
for( i in 1:length(unique(top40[,1])))
{
    sname<-unique(top40[,1])[i]
    temp<-top40[top40$sampleName==sname,]
    temp<-temp[order(temp$xMembers, decreasing=T),]
    temp<-temp[c(1:topn),]
    dto<-rbind(dto, temp)
}
idis.SP.IgG.top40<-dto


conditions<-read.table(here("Data","sampleInfo.txt"), sep="\t", header=T)
rownames(conditions)<-conditions$ID
conditions$treatment<-factor(conditions$treatment, levels=c("PBS","OVA", "OVA+PorB", "OVA+CpG", "OVA+Alum" ))

#merge meanmutation and idis
mi.BM.IgM.top40<-merge(idis.BM.IgM.top40, meanMu.BM.IgM.top40, by.x=c("sampleName", "cloneID"), by.y=c("sampleName", "CloneID"))
mi.BM.IgG.top40<-merge(idis.BM.IgG.top40, meanMu.BM.IgG.top40, by.x=c("sampleName", "cloneID"), by.y=c("sampleName", "CloneID"))
mi.SP.IgM.top40<-merge(idis.SP.IgM.top40, meanMu.SP.IgM.top40, by.x=c("sampleName", "cloneID"), by.y=c("sampleName", "CloneID"))
mi.SP.IgG.top40<-merge(idis.SP.IgG.top40, meanMu.SP.IgG.top40, by.x=c("sampleName", "cloneID"), by.y=c("sampleName", "CloneID"))

#by sample

data.array<-list(mi.BM.IgM.top40, mi.BM.IgG.top40, mi.SP.IgM.top40, mi.SP.IgG.top40)
name.array<-c("Bone Marrow (IgM)", "Bone Marrow (IgG)", "Spleen (IgM)", "Spleen (IgG)");
#save(data.array, conditions, name.array, file="intraClonalDiversityDataArray_top15.RData")
#save(data.array, conditions, name.array, file="intraClonalDiversityDataArray_top30.RData")
save(data.array, conditions, name.array, file=here(data.dir,"intraClonalDiversityDataArray_top40.RData"))

#first do linear regression, 


mi.temp<-data.array[[1]]   #<-------------build threshold based on 
x<-rlm(MeanMuFreq~idi, data=mi.temp,psi=psi.huber)
plot(mi.temp$idi, mi.temp$MeanMuFreq, type="p")
cx<-seq(0, 0.05, by=0.001)
cy<-cx*x$coefficients[2]+x$coefficients[1]
lines(cx,cy, col=2, lwd=2, lty=2)
cy2<-cx*x$coefficients[2]+0.01
lines(cx,cy2, col=3, lwd=2, lty=3)
#pdf("intraClonalDiversity.pdf")
op<-par(mfrow=c(2,2), mar=c(3,2.5,2,0.2), mgp=c(1.5,0.4,0.0))
for (i in 1:4)
{
    plot(x=c(0.00, 0.034), y=c(0,0.04), type="n", xlab="Intra-Clonal Dissimilarity", 
                                                    ylab="Mean Mutation Frequency",main=name.array[i])
    lty=2

    mi.temp<-data.array[[i]]#mi.SP.IgM.top20
            #actually by sample
    sname<-unique(mi.temp[,1])
    for(j in 1:length(sname))
    {
        #get treatment 
        points(mi.temp[mi.temp$sampleName==sname[j],3], mi.temp[mi.temp$sampleName==sname[j],8],
                        col=as.integer(conditions[sname[j],"treatment"]),
                    pch=as.integer(conditions[sname[j],"treatment"]) )
    }
    cx<-seq(0, 0.05, by=0.001)
    cy<-cx*x$coefficients[2]+x$coefficients[1]
    lines(cx,cy, col=2, lwd=2, lty=2)
    cy2<-cx*x$coefficients[2]+0.01
    lines(cx,cy2, col=3, lwd=2, lty=3)
    legend(x=0.025, y=0.015, legend=unique(as.character(conditions$treatment)), col=unique(as.integer(conditions$treatment))
                    , pch=unique(as.integer(conditions$treatment)), cex=0.8)
}
par(op)

#now let's do ggplot2 with polygon range around the points
#op<-par(mfrow=c(2,2), mar=c(3,2.5,2,0.2), mgp=c(1.5,0.4,0.0))
tiff(here("R_code/Figure6","intraClonalDiversity_threshold.tiff"), width=800, height=700)
gp<-vector("list", length=4);
for (i in 1:4)
{
    
    mi.temp<-data.array[[i]][,c("sampleName", "idi", "MeanMuFreq")]
    mi.temp<-cbind(mi.temp, conditions[mi.temp$sampleName, c("tissue", "treatment") ])
    gp[[i]]<-ggplot(data=mi.temp, aes(x=idi, y=MeanMuFreq, color=treatment, group=treatment))+
        geom_point(aes(shape=treatment))+lims(x=c(0.00, 0.034),y=c(0,0.04))+theme(legend.position = c(0.79, 0.18))+
        #geom_encircle(expand=0)+#stat_ellipse(level=0.5) +  ###add the circles 
        geom_abline(slope=x$coefficients[2], intercept=x$coefficients[1], color="red")+  #<-add lines
        geom_abline(slope=x$coefficients[2], intercept=0.01, color="green", linetype="dashed")+ #<-add lines
        labs(x="Intra-Clonal Dissimilarity", y="Mean Mutation Freq", title=name.array[i])
   
}
figure<-ggarrange(gp[[1]], gp[[2]], gp[[3]], gp[[4]], ncol=2, nrow=2)
#par(op)
figure
dev.off()

#now based on threshold, we categorize the points.
categorizeCloneIntraDiversity<-function(dat, slope, intercept)
{
    x<-dat[,1]
    y<-dat[,2]
    return(c("-"=sum(slope*x+intercept>y)/length(x), "+"=1-sum(slope*x+intercept>y)/length(x)))
}

categorizeCloneIntraDiversity(mi.temp[,c("idi", "MeanMuFreq")], x$coefficients[2],0.01)
#categorizeCloneIntraDiversity_unselected(mi.temp[,c("idi", "MeanMuFreq")], slope=x$coefficients[2],0.01, 0.01)
#prepare the data for each tissue, isotype and treatment 
intraD<-data.frame()
tissues<-c("Bone Marrow", "Bone Marrow", "Spleen", "Spleen");
isotypes<-c("IgM", "IgG", "IgM", "IgG")
for(i in 1:4)
{
    mi.temp<-data.array[[i]]#mi.SP.IgM.top20
            #actually by sample
    sname<-unique(mi.temp[,1])
    temp<-data.frame()
    for(j in 1:length(sname))
    {
        #get treatment 
        id<-categorizeCloneIntraDiversity(mi.temp[mi.temp$sampleName==sname[j],c("idi", "MeanMuFreq")], x$coefficients[2],0.01)
        
        temp<-rbind(temp,data.frame(neg=id["-"], pos=id["+"], sampleName=sname[j]))
    }
    #temp$sampleName<-sname
    temp$tissue<-tissues[i]
    temp$isotype<-isotypes[i]
    intraD<-rbind(intraD, temp)
}
#add the treatment factors
intraD$treatment<-conditions[as.character(intraD$sampleName), "treatment"]
rownames(intraD)<-paste0(intraD$sampleName, intraD$isotype)
intraD$logRatio<-log(intraD$pos/intraD$neg)
intraD.comp<-acomp(as.matrix(intraD[,1:2]))
intraD.comp<-zeroreplace(intraD.comp, d=c(1/topn, 1/topn), a=1/3)
bmatrix<-matrix(0, nrow=1, ncol=2)
bmatrix[1,]=c( -1,1)
#bmatrix[2,]=c(rep(-1, 3), 1,rep(0,1))

#intraD.ilr<-ilr(intraD.comp)
intraD.ilr<-ilr(intraD.comp, V=buildBalanceBase(bmatrix))
intraD.clr<-clr(intraD.comp)
intraD$isotype<-factor(intraD$isotype, levels=c("IgM", "IgG"))
intraD$tissue<-factor(intraD$tissue, levels=c( "Spleen","Bone Marrow"))
intraD$treatment<-factor(intraD$treatment, levels=c("PBS", "OVA", "OVA+PorB", "OVA+CpG", "OVA+Alum"))
intraD$ilr<-intraD.ilr
save(intraD, intraD.ilr, file=here(data.dir,"intraClonalDiversity_intraD_top40.RData"));
