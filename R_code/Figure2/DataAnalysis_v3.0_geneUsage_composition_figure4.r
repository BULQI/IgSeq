#R code to analyze the heavy chain sequencing data
#
##=====-----make-up notes: 8/11/2023; this is for make up the previous run
#    this is the latest one that run on all sequences (not just zero or no zero mutation sequences.)
#
##=======================-------------------------------------
#10/26/2021   ---- updated
#  change different conditions for testing:
#    zero mutation sequences.
#   remove the IGHV reads with too few (originally were <250, now it is <700 and also try <1000 total reads.)
#   Scaled and CENTERED, before doing PCA !!!
#   we also spent much time testing the best way to replace zero values since we need to do scale and ratio --a (see line 114 and line 120.)
# final paramter: no scale of zero replace value and ratio of 2/3
#           remove <700 total number of reads IGHV
#           centered, but not scaled for pca and clr data.
##===========================
#  4/4/2021 ----
## modify the code to plot loading and pcs separately and also 
#       change to VC1 
#   output the top 10 genes both directions 
#
#######
#3/12/21, add code to plot the individual lines of top 10 from both end showing the treand.
# Feng
#====
#3/5/2021. copied from the v1.0 version. in this version we will do heatmap on the "raw" data for each
# IGHV to confirm the PC effects.
#
#=====================
#now this is the one doing version, to show the individual effect with raw data like figure 4e, f.  
#meaning to confirm the PC effects by the raw data.(raw percent data)
# next version v2.0, try to show the pc effect by heatmap (centered and scaled effects for each IgV genes)
#====
#updatecd 2/12/2021
# we also add code to remove zoro counts and also we need to do contrast for doing type III anova.
# type I and II we don't have to do that.
# need to set up factors to do post hoc correctly with contrast. 
#-================
#copied from /newPipeline/DataAnalysis_v1_geneUsage_composition.R2
# read the saved data 
# plot for figure 4.
#-----------------------------------
#copied from LBSeq8
#  doing the gene usage analysis for wl03R2 new pipeline.
#      10/5/2020
#----------------------------- 
#LBSeq 8, copied from WL03 R2 isotype folder.
#   doing gene usage analysis.
#   In here we read the data from previous save RData. 
#  DataAnalysis_v1.0_geneUsageData.R (previously DataAnalysis_v1.0_geneUsage.R)
#     -------2/18/2020
#-----------------
# WL03, IgSeq
#Feng @BU ------1/30/2020.
# this one is copied from previous work WL03 analysis
#In here we read the data from previous save RData. 
#  DataAnalysis_v1.0_geneUsageData.R (previously DataAnalysis_v1.0_geneUsage.R)
#then we do analysis using composition.
#
#--------------the below part is from the previous version. Read them in caution.---
#Plan 1) do gene usage to compare between groups
#	  2) 
#-------
##5/28/2019
#  1) using composition to do the analysis
#		i) linear model /anova for difference
#		ii) tSNE for overa pattern
#		iii)PCA for pattern, each component for different pairwise different  
#		iv) cluster for grouping.
#  2) reading the data from each group and each tissue.
#		this is done in file "DataAnalysis_v1.0_geneUsage.R"
#		two lists BM.usage and SP.usage are saved. They are two
#			lists contains the v, d, j gene usage and conditoins and totalIgs
##------------
##6/3/2019
#1) new data. new analyzed data set. because the old processing have lots 
#==========================
library(plyr)
library(compositions)
library(ggplot2)
library(MASS)
library(here)
library(ggrepel)
library(dplyr)
library(emmeans)
library(ggpubr)
library(car)
library(ggfortify)

#set working directory first
#do SP first
#setwd("E://feng//LAB//MSI//MouseLungProject//LBSeq8//IgTotal")
#setwd("/home/feng/Windows/windowsD/feng/LAB/MSI/MouseLungProject/LBSeq8/IgTotal")
#setwd("/home/feng/Feng/LAB/Wetzler_PorinB/wl03R2/newPipeline/")
data.dir<-"Data/Figure2"
output.dir<-"R_code/Figure2"
load(here(data.dir, "sample.geneUsage.RData"))  # save in newPipeline/DataAnalysis_v1_geneUsage_data.r  <-all sequences
                                                                               #and DataAnalysis_v1_geneUsage_data_noZeroMu. <-no zero sequences
#load("sample.geneUsage_noZeroMu.RData")                                                                               
#so far 4 lists loaded, BM/SP IgG/IgM and conditions

#============================
#   IgG/IgM, BM/Spleen, Treatment
#++++++++++++++++++++++++++++
BM.IgG<-BM.susage[["IgG"]][["VGene"]]
BM.IgM<-BM.susage[["IgM"]][["VGene"]]
BM<-merge(BM.IgG, BM.IgM, by="VGene", all=T, suffixes=c("_IgG", "_IgM"))

SP.IgG<-SP.susage[["IgG"]][["VGene"]]
SP.IgM<-SP.susage[["IgM"]][["VGene"]]
SP<-merge(SP.IgG, SP.IgM, by="VGene", all=T, suffixes=c("_IgG", "_IgM"))

#merge
isotype<-merge(SP, BM, by="VGene", all=T, suffixes=c("_SP", "_BM"))
#check to replace NA with zero
isotype.noNA<-isotype
    for(i in 1:dim(isotype.noNA)[1])
    {
        isotype.noNA[i,is.na(isotype.noNA[i,])]<-0;
    }
usage.VGene<-isotype.noNA
usage.VGene.clean<- usage.VGene[-(which(apply(usage.VGene[,-1], 1, sum)<750)),]
#usage.VGene.clean<- usage.VGene[-(which(apply(usage.VGene[,-1], 1, sum)<1000)),]
usage.VGene.t<-t(usage.VGene.clean[,-1])
colnames(usage.VGene.t)<-(usage.VGene.clean[,1])
dl<-apply(usage.VGene.clean[,-1],2, sum)
dl<-1/dl
dl<-rep(dl,dim(usage.VGene.t)[2])
dl<-matrix(dl, ncol=dim(usage.VGene.t)[2], nrow=dim(usage.VGene.t)[1], byrow=F)
min_dl<-log(min(dl))
max_dl<-log(max(dl))
scale_dl<- log(max(dl)/min(dl))#log(2) #log(10)#
range_dl<-max_dl-min_dl

dl<-(log(dl)-min_dl)/range_dl *scale_dl+min_dl 
dl<-exp(dl)
usage.VGene.comp<-acomp(usage.VGene.t)
usage.VGene.comp<-zeroreplace(usage.VGene.comp, d=dl, a=2/3)

usage.VGene.clr<-clr(usage.VGene.comp)
#usage.VGene.ilr<-ilr(usage.VGene.comp)
dt.clr<-data.frame(usage.VGene.clr)
#prepare conditions to combine
temp<-conditions[conditions$tissue=="Spleen",]
temp$isotype<-"IgG"
temp<-temp[order(temp$snum),]
conditions.all<-temp
temp<-conditions[conditions$tissue=="Spleen",]
temp$isotype<-"IgM"
temp<-temp[order(temp$snum),]
conditions.all<-rbind(conditions.all, temp)

temp<-conditions[conditions$tissue=="Bone Marrow",]
temp$isotype<-"IgG"
temp<-temp[order(temp$snum),]
conditions.all<-rbind(conditions.all, temp)
temp<-conditions[conditions$tissue=="Bone Marrow",]
temp$isotype<-"IgM"
temp<-temp[order(temp$snum),]
conditions.all<-rbind(conditions.all, temp)

dt.clr<-cbind(dt.clr, conditions.all)

    ##########################
    #     new notes:: here!!! 
    #     do we need to save the data here???
    #     need to check later
    ############################
#//now save the data, so they can be used by clustering and MFA in figure 6.
#setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/figure4/")
#   save(file="geneUsage_all_clr.RData", dt.clr)

##end of checking for saving section 




#now start doing the visualization
#observation is sample/subject
#variable is the vgene family

  #------------Vgene	

   dt.clr.data<-dt.clr[dt.clr$treatment!="OVA+PorB",c(1:(dim(dt.clr)[2]-6))]
   #replace name
    colnames(dt.clr.data)<-sub("^IGHV", "V",colnames(dt.clr.data))
    colnames(dt.clr.data)<-sub(".", "-",colnames(dt.clr.data),fix=T)
   dt.clr.all<-dt.clr[dt.clr$treatment!="OVA+PorB",]
   dt.clr.all$treatment<-factor(dt.clr.all$treatment, levels=c("PBS", "OVA", "OVA+CpG", "OVA+Alum"))
	pca.VGene.clrd<-prcomp(dt.clr.data, scale=F, center=T)


  save(pca.VGene.clrd, dt.clr.data,  
      file=here(data.dir,"pcaData_cronbachAlapha.RData"))#save the data so to do crobach alapha for stability.
                ###==== see the R file PC_consistency_CronbachA.R
    #get scaled data ready for doing trend for confirming 
	dt.clr.data.scaled<-scale(dt.clr.data, scale=F, center=T)	
   #pca.VGene.ilr<-pcaCoDa(data.frame(usage.VGene.ilr))
   #pca.VGene.ilrd<-prcomp(data.frame(usage.VGene.ilr), scale=T)
#   setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/figure4/")
   png(file=here(output.dir,"pc_AllPooled_dlScaleNo.png"), width=500, height=500)
   num.pc<-30
   plot(seq(1,28,1),(pca.VGene.clrd$sdev^2/sum(pca.VGene.clrd$sdev^2))[seq(1,28,1)],type="b",
		main="PC contribution", xlab="PCs", ylab="portion",ylim=c(0,1.05), col=2, pch=15
	)
   lines(seq(1,num.pc,1),summary(pca.VGene.clrd)$importance[3,seq(1,num.pc,1)],type="b")
   lines(c(0,num.pc),c(0.05,0.05), lty=2, lwd=2,col=3)
   lines(c(6,6),c(0,1.9), lty=2, lwd=2,col=2)
   legend(x=15, y=0.75, legend=c("individual","accumative" ), pch=c(15,1), col=c(2,1), lty=c(1,1), cex=1.4)
   text(x=7.65, y=0.85,label="PC 6", cex=1.6)
   text(x=22, y=0.10, label="5% contribution", cex=1.5)
   dev.off()
   
   #plot the above again using ggplot2
   dt1<-data.frame(con.ac=summary(pca.VGene.clrd)$importance[3,seq(1,num.pc,1)], ind=seq(1,num.pc,1), 
                                                    con=(pca.VGene.clrd$sdev^2/sum(pca.VGene.clrd$sdev^2))[seq(1,num.pc,1)])
                                                    
   cn<- ggplot(data=dt1, aes(x=ind, y=con.ac))+
            geom_bar(data=dt1, aes(x=ind,y=con), stat="identity", fill="turquoise2")+
            scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))+
            geom_point()+ylab("Contribution (%)")+xlab("PCs")+ggtitle("PC contributions")+theme_bw(base_size=14)+
            geom_line(linetype=2)+
            geom_segment(x=6,xend=6, y=0,yend=1.1, colour="red", linetype=3 )+
            geom_segment(x=0,xend=50, y=0.05,yend=0.05, colour="red", linetype=3 )+
            annotate(geom="text", label=c("Accumulative", "Individual"), x=c(25,25),y=c(0.9, 0.1), colour="red",size=7 )+
            theme(plot.title= element_text(margin = margin(b = -20), hjust=0.01))
	#   lines(c(0,64),c(0.015,0.015), lty=2, lwd=2,col=3)

g<-autoplot(pca.VGene.clrd,x=1,y=2,data=dt.clr.all, #colour="isotype",shape="tissue",
				label=F ,frame=F,  size=-1,#frame.type="t"
			loadings=T, loadings.label=T, loadings.label.size=4, loadings.colour="orange", loadings.label.colour="red"
		) +
		#scale_colour_manual(name="isotype", values= c("forest green", "red3"))   +# c("forest green", "red3", "dark blue"))+
  ggtitle("PCA VGene Usage (loadings)")+labs(y="PC2", x="PC1")+
  theme_bw(base_size=13) + xlim(-0.2, 0.2)+ylim(-0.2, 0.2)+
  theme(legend.position = c(0.95,0.05),
                                    legend.justification=c(1,0),plot.title = element_text(margin = margin(b = -20)))
#ggbiplot(pca.VGene.clrd)
g2<-autoplot(pca.VGene.clrd,x=1,y=2,data=dt.clr.all, colour="isotype",shape="tissue",
				label=F,frame=F,  size=4,#frame.type="t"
			loadings=F, loadings.label=F, loadings.label.size=4, loadings.colour="orange", loadings.label.colour="red"
		) +
		#scale_colour_manual(name="isotype", values= c("forest green", "red3"))   +# c("forest green", "red3", "dark blue"))+
  ggtitle("PCA VGene Usage")+labs(y="PC2", x="PC1")+
  theme_bw(base_size=13) + #xlim(-0.32, 0.3)+ylim(-0.32, 0.3)+
  theme(legend.position = c(0.01,0.3),legend.background = element_rect(color = "lightblue", linetype = "solid",size=0.5),
                                    legend.justification=c(0,1),plot.title = element_text(margin = margin(b = -20)))
#ggbiplot(pca.VGene.clrd)
h<-autoplot(pca.VGene.clrd,x=1,y=2,data=dt.clr.all, colour="treatment",shape="isotype",
				label=F ,frame=F, size=4,#frame.type="t", size=5
			loadings=F, loadings.label=F, loadings.label.size=4, loadings.colour="orange", loadings.label.colour="orange"
		) + scale_colour_discrete(name="Immunization")+
		#scale_colour_manual(name="treatment", values= c("black","pink","forest green", "red3" ))   +# c("forest green", "red3", "dark blue"))+
  ggtitle("PCA VGene Usage")+labs(y="PC2", x="PC1")+
  theme_bw(base_size=13) +
  theme(legend.position = c(0.01,0.30),legend.justification=c(0,1),plot.title = element_text(margin = margin(b = -20))) 

#ggbiplot(pca.VGene.clrd)
J<-autoplot(pca.VGene.clrd,x=5,y=6,data=dt.clr.all, colour="treatment",
				label=F ,frame=F, size=4,#frame.type="t", size=5
			loadings=T, loadings.label=T, loadings.label.size=4, loadings.colour="orange", loadings.label.colour="orange"
		) +
		#scale_colour_manual(name="treatment", values= c("black","pink","forest green", "red3" ))   +# c("forest green", "red3", "dark blue"))+
  ggtitle("PCA VGene Usage")+labs(y="PC6", x="PC5")+
  theme_bw(base_size=13) +
  theme(legend.position = c(0.95,0.95),legend.justification=c(1,1),plot.title = element_text(margin = margin(b = -20))) 

#save for plotting in the future.
save(file=here(data.dir,"pca_loading.RData"), g, g2, h, J, pca.VGene.clrd, cn, dt1 , dt.clr.all)
  #now do linear regressoin/anova
  png(file=here(output.dir,"pca_represent_isotypeTissue.png"),width=500, height=500)
  g2
  dev.off()
  ####now let's do anova 
	sample.conditions<-conditions.all[conditions.all$treatment!="OVA+PorB",]#sample.usage[["conditions"]]

	#now we do linear regression and also we don't add the factor of mouse number, since
	#it results in an unbalanced design (some groups of combination don't exist; for example,
	#there is the group with the treatment of saline and mouse #1, but no the one with
	#IVA and mouse #1). This might lead to a singular fitting.
	pc.uplim<-8
  
#-----this following line is to set up the correct contrast 
#---- to do ANOVA for using library(car), please be careful.
	options(contrasts = c("contr.sum", "contr.poly"))

#now we need to do individual ANOVAs on each PCs 
pc.uplim<-6
i<-1
model.ind.anova<-NULL

xx<-cbind(pca.VGene.clrd$x[,i], sample.conditions[,c("treatment","isotype","tissue")])
colnames(xx)[1]="y"
model.ind<-lm(y~treatment*isotype*tissue, data=xx)
Anova(model.ind, type=2)


for(i in 1:pc.uplim)
{
	model.ind<-lm(pca.VGene.clrd$x[,i]~sample.conditions$treatment*sample.conditions$isotype*sample.conditions$tissue);#+sample.conditions$Tissue*sample.conditions$Memory);
	temp<-cbind(Anova(model.ind, type=2)[c(1:7),],PCs=i, effect=rownames(Anova(model.ind, type=2))[1:7])
    #temp<-cbind(Anova(model.ind, type=3)[c(1:8),],PCs=i, effect=rownames(Anova(model.ind, type=3))[1:8])
	rownames(temp)<-NULL
	temp$effect<-gsub("sample.conditions\\$","",temp$effect)
	temp$p.adj.ind<-p.adjust(temp$"Pr(>F)", method="BH")
	model.ind.anova<-rbind(model.ind.anova,temp)
	
}	
model.ind.anova$p.adj<-p.adjust(model.ind.anova$"Pr(>F)", method="fdr")
#sort(model.ind.anova$p.adj)

#now we start doing the Post hoc tests
#                              for PC1       <-------------------

pcs<-1
dt.mod<-data.frame(ilr=pca.VGene.clrd$x[,pcs], treatment=sample.conditions$treatment, isotype=sample.conditions$isotype, tissue=sample.conditions$tissue)
dt.mod$treatment<-factor(dt.mod$treatment, levels=c("PBS","OVA", "OVA+CpG", "OVA+Alum"))
mod<-lm(ilr~treatment*isotype*tissue, data=dt.mod)
em_tr<-emmeans(mod, ~tissue|treatment*isotype)
#em_tr<-emmeans(mod, ~isotype|treatment*tissue)
em_tr<-emmeans(mod, ~treatment|tissue*isotype)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
    "OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), "OVA - PBS"=c(-1, 1,0,0)
 # "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#   "Bone Marrow - Spleen" = c(1,-1)
  )

ec<-emmeans::contrast(em_tr, Set1, adjust='FDR')


Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
    "OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0),
    "OVA/Alum - OVA" = c(0,-1,0,1), "OVA/CpG - OVA"=c(0,-1,1,0)
 # "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#   "Bone Marrow - Spleen" = c(1,-1)
  )

ec.OVA<-emmeans::contrast(em_tr, Set1, adjust='FDR')

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
    "OVA - PBS" = c(-1,1,0,0), #"OVA/CpG - PBS"=c(-1,0,1,0),
    "OVA/Alum - OVA" = c(0,-1,0,1), "OVA/CpG - OVA"=c(0,-1,1,0)
 # "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#   "Bone Marrow - Spleen" = c(1,-1)
  )

ec.OVA2<-emmeans::contrast(em_tr, Set1, adjust='FDR')


#plot the difference
pc1.tr<-ggplot(data=dt.mod, aes(y=-1*ilr, x=tissue, color=treatment))+
        #geom_point(aes(shape=treatment))+#lims(x=c(0.00, 0.034),y=c(0,0.04))+theme(legend.position = c(0.79, 0.18))+
        geom_boxplot(aes(shape=tissue), notch=F, notchwidth=0.1)+
        #lims(x=c(0.00, 0.034),y=c(0,0.04))+
        guides(shape="none")+
        labs( y="PC1 (IGHV usage)")+#, title=name.array[i]) +
       facet_grid(.~isotype)+ 
       #geom_hline(data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"), aes(yintercept=yintercept),
       #                     linetype="dashed", color="red")
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=1.7, xend=1.9, y=10.5, yend=10.5), colour="black",size=0.9) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
#         geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Bone Marrow"),
#                    aes(x=0.7, xend=1.3, y=1.5, yend=1.5),color="black", size=0.9) +   
         geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                   aes(x=0.7, xend=1.1, y=8.5, yend=8.5),color="black", size=0.9) + 
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                   aes(x=0.7, xend=1.3, y=9.5, yend=9.5),color="black", size=0.9) + 
         #geom_point(data=data.frame(treatment=c(2.2, 2.3,2.5,2.6, 1.3 ),ilr=c(0.47,0.47, 0.83,0.83, 0.2), isotype="IgG", tissue="Bone Marrow"),
         #           color="black", size=2, shape=8) +   
         geom_text(data=data.frame(treatment=c(1.0,1.1,1.0),ilr=c(-9, -10.0, -11.0), isotype=c("IgG",
         "IgG", "IgG"), tissue=c(0.9,1.0,1.8)),
                    color="black", size=4, aes(label=c("p<0.01", "p<0.005", "p<0.01")))+

        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=1.9, xend=2.1, y=11.5, yend=11.5), colour="black",size=0.9) +
         geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                   aes(x=1.9, xend=2.3, y=12.5, yend=12.5),color="black", size=0.9) + 
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                   aes(x=0.9, xend=1.1, y=10.5, yend=10.5),color="black", size=0.9) + 
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                   aes(x=0.9, xend=1.3, y=11.5, yend=11.5),color="black", size=0.9) +
         geom_text(data=data.frame(treatment=c(1.0,1.1,1.0,1.1),ilr=c(-11.0, -12.0, -12.0,-13), isotype=c("IgG",
         "IgG", "IgG","IgG"), tissue=c(1.0,1.1,2.0, 2.1)),
                    color="black", size=4, aes(label=c("p<0.05", "p<0.01", "p<0.001","p<0.05")))+
        theme(legend.position = c(0.8, 0.78),text= element_text(size=18),axis.title.x=element_blank())
pdf(file=here(output.dir,"pc1_treatment_dlnoScale_dataScaleCentered.pdf"))
    pc1.tr
dev.off()
save(file=here(data.dir,"figure4_pc1_draw.RData"), dt.mod, pc1.tr)


#now plot the contributions
#two kinds of contribution
# PC loadings/rotation matrix

#pcs<-1

#contributions of different variables to PC1, square the egein vector (without corrected by eigen value)
contr<-abs(pca.VGene.clrd$rotation[,pcs])*pca.VGene.clrd$rotation[,pcs]/sum(pca.VGene.clrd$rotation[,pcs]^2)
contr<-sort(contr)
contr.names<-names(contr)
rotation<-pca.VGene.clrd$rotation[,pcs]
save(file=here(data.dir,"gene_contribution_ordered_PC1.RData"), 
    contr, rotation)

dt.contr<-data.frame(contr=contr, IGHVs=contr.names, index=seq(1,length(contr)))
dt.contr.pos<-dt.contr[dt.contr$contr< -0.015 | dt.contr$contr> 0.015,]

contr.LD1<-ggplot(data=dt.contr, aes(x=index,y=contr))+
    geom_point()+scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))+ theme_bw(base_size=15)+
    geom_segment(x=-10, xend=140, y=0, yend=0, colour="blue", linetype=2,size=0.65 )+
    geom_segment(x=-10, xend=140, y=0.02, yend=0.02, colour="blue", linetype=3,size=0.4 )+
    geom_segment(x=-10, xend=140, y=-0.02, yend=-0.02, colour="blue", linetype=3,size=0.4 )+
    geom_text_repel(data=filter(dt.contr, contr< -0.020 | contr >0.020), aes(label=IGHVs),max.overlaps = Inf)+
    ggtitle(paste0("IGHV contribution to PC", pcs))+ xlab("")+ ylab("% contribution")+
    theme(plot.title = element_text(margin = margin(b = -20),size=15))
  
###contribution of variable to each PC
#ig<-c("V14-1" , "V15-2", "V5-4" ,    "V1-81")#,"V2-2")
ig<-c("V15-2" , "V5-9", "V2-2" ,    "V5-4")#,"V2-2")

contr<-(pca.VGene.clrd$rotation[ig,])*matrix(rep(pca.VGene.clrd$sde, length(ig)), ncol=length(pca.VGene.clrd$sde), nrow=length(ig),byrow=T)#/sum(pca.VGene.clrd$rotation[,pcs]^2)
contr<-abs(contr)*contr/apply(contr^2, 1, sum)
#contr<-sort(abs(contr))
#contr.names<-names(contr[order(contr)])
cb<-vector("list", length=length(ig))
for(i in 1:length(ig)){
dt.contr<-data.frame(contr=contr[i,], IGHVs=colnames(contr), index=seq(1,dim(contr)[2]))
#dt.contr.pos<-dt.contr[dt.contr$contr< -0.015 | dt.contr$contr> 0.015,]

cb[[i]]<-ggplot(data=dt.contr, aes(x=index,y=contr))+
    geom_bar(stat="identity", fill="turquoise2")+#ylim(-0.85, 0.1)+
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))+
    #geom_segment(x=-10, xend=140, y=0.02, yend=0.02, colour="blue", linetype=2,size=0.35 )+
    geom_segment(x=-10, xend=140, y=-0.00, yend=-0.0, colour="blue", linetype=2,size=0.35 )+theme_bw(base_size=15)+
    #geom_segment(x=-10, xend=140, y=-0.02, yend=-0.02, colour="blue", linetype=3,size=0.4 )+
    geom_text_repel(data=filter(dt.contr, contr< -0.25 | contr >0.25), aes(label=IGHVs),max.overlaps = Inf,size=5)+
    ggtitle(paste0(ig[i]," contribution to PC", pcs))+ xlab("PCs")+ ylab("contribution (%)")+
        
    theme(plot.title = element_text(margin = margin(b = 0), hjust=0.5))
    #cb<-c(cb, temp)
}

tiff(file=here(output.dir,"PC1contributions_figure4.tiff"), 
    width=900, height=630)
ggarrange(contr.LD1, 
                ggarrange(cb[[1]],cb[[2]],cb[[3]], cb[[4]],  
                          ncol = 2, nrow = 2,
                          heights = c(2, 2)),  
          ncol = 1, nrow = 2, labels=c("A","B"),
          heights = c(2, 3))
dev.off()


        #put genes in order
        genes.contr2<-c("V5-9","V14-1" , "V1-63", "V3-6", "V1-20", 
                    "V5-4","V9-1", "V8-8",    "V14-3" ,"V2-2")
#        eff<-t(dt.clr.data.scaled[, genes.contr2])
        #now build the column order
        #BM first 
        conditions$treatment<-factor(conditions$treatment, levels=c("PBS", "OVA", "OVA+PorB", "OVA+CpG","OVA+Alum"))
#        conditions.noPorB<-conditions[conditions$treatment!="PorB",]
        bm.con<-conditions[conditions$tissue=="Bone Marrow",]
        bm.con<-bm.con[bm.con$index.fileOrder,]
        bm.con<-bm.con[bm.con$treatment!="OVA+PorB",]
        bm.IgG.cnames<-paste0(bm.con$ID, "_IgG" )

        #Spleen
        sp.con<-conditions[conditions$tissue=="Spleen",]
        sp.con<-sp.con[sp.con$index.fileOrder,]
        sp.con<-sp.con[sp.con$treatment!="OVA+PorB",]
        sp.IgG.cnames<-paste0(sp.con$ID, "_IgG" )
        
        bm.IgM.cnames<-paste0(bm.con$ID, "_IgM" )
        #Spleen
        sp.IgM.cnames<-paste0(sp.con$ID, "_IgM" )
        
        
      #plot lines
      #plot(c(1,5), c(-2,2))
      #'@title make data read for pictures for line pattern (single tissue with one isotype)
      #'@ description to do for only tissue isotype combination, like doing Bone marrow IgG.
      #'            in this function we taking the data frame of all samples/IGHVs and pick the ones with 
      #'            a correct order so to plot only the top a few to show the pattern 
        #'          It assume for each tissue and isotype type 
      #'@param eff  this is the data from holding the clr transformed data to show pattern
      #'         it assume rows of samples and cols of IgHV. It is for one  isotype and tissue combination
       #'       so that we can summarize it in regards to the treatment.
      #'@param orderByName a character array containing the name of IGHV names in order 
      #'            the order is sorted according to the contribution to each specific PC. 
      #'@param numOfElement a scalar specifying the number of elements to pick. we pick from both the top
      #'            and bottom, and therefore we will pick 2x numOfElement columns in total.
      #'@param tissueName name of the tissue of the output data frame. 
      #'@param isotypeName name of the isotye of the output data frame
      #'@param conds condition table holding all the necessary meta information. 
      #'     in this case it at least holding treatment column that is compatible with the input 
      dataLinePattern_single<-function(d, orderByName, numOfElement=10, 
                tissueName, isotypeName, conds, signs_flipped=F )
      {
                ###----- Bone Marrow IgG
       dtt<- d[,orderByName]
       dta<-data.frame()
       for(i in 1:length(orderByName[1:numOfElement])){
                   dtt.sum<-aggregate(dtt[,i], by=list(conds$treatment), mean)
                   colnames(dtt.sum)[1]<-"treatment"
                   dtt.sum$VGene<-orderByName[i]
                  dta<-rbind(dta, dtt.sum)
                  #lines(dtt.sum[,1], dtt.sum$x)
              }
              dta$tissue<-tissueName
              if(signs_flipped){
                dta$sign<-"positive"
              }
              else{
                dta$sign<-"negative"
              }
              #dta$isotype<-isotypeName
                temp<-data.frame()
      for(i in length(orderByName):(length(orderByName)-numOfElement+1)){
                   dtt.sum<-aggregate(dtt[,i], by=list(conds$treatment), mean)
                  colnames(dtt.sum)[1]<-"treatment"
                  dtt.sum$VGene<-orderByName[i]
                  temp<-rbind(temp, dtt.sum)
                  #lines(dtt.sum[,1], dtt.sum$x, col=2)
              }
              temp$tissue<-tissueName
              temp$sign<-"positive"
              if(signs_flipped){
                temp$sign<-"negative"
              }
              else{
                temp$sign<-"positive"
              }
              
              dta<-rbind(dta, temp)
              dta$isotype<-isotypeName
              dta$sign<-factor(dta$sign, levels=c("positive", "negative"))
              dta$isotype<-factor(dta$isotype, levels=c("IgM","IgG"))
              dta$tissue<-factor(dta$tissue, levels=c("Bone Marrow","Spleen"))
                return(dta)
      }      
      
      #Bone marrow IgG
      nElement<-9;
      #dt.clr.data.scaled<-dt.clr.data
      eff<-dt.clr.data.scaled[bm.IgG.cnames,]
      xt<-dataLinePattern_single(d=eff, orderByName=contr.names, numOfElement=nElement,
                    tissueName="Bone Marrow", isotypeName="IgG", conds=bm.con, signs_flipped=T)
        #Bone marrow IgM
      eff<-dt.clr.data.scaled[bm.IgM.cnames,]      
      temp<-dataLinePattern_single(d=eff, orderByName=contr.names, numOfElement=nElement,
                    tissueName="Bone Marrow", isotypeName="IgM", conds=bm.con, signs_flipped=T)
                    
         xt<-rbind(xt,temp)
         xt2<-xt
             xt2$geneIsotype<-paste0(xt$VGene, xt$isotype)
             
             trendline<-aggregate(xt2$x, by=list(treatment=xt2$treatment, isotype=xt2$isotype, sign=xt2$sign, tissue= xt2$tissue), mean)        
            trendline$geneIsotype=trendline$isotype
        #plot 
        pc1.line<-ggplot(data=xt2, aes(x=treatment, y=x, group=geneIsotype,colour=isotype))+
            geom_line(linetype=2)+
            geom_point()+
            geom_line(data=trendline, aes(y=x,x=treatment,colour=isotype), linetype=1,size=1.5)+
            facet_grid(sign~.)+theme_bw(base_size=15)+theme(axis.title.x=element_blank())
     
     
     #       pdf("pc4_lines_dlScaledLog2.pdf")
            pc1.line
      #      dev.off()
       #save data for the future
save(file=here(data.dir,"figure4_pc1_lines.RData"), xt2 , pc1.line , trendline);       
            #spleen IgG
      eff<-dt.clr.data.scaled[sp.IgG.cnames,]
      temp<-dataLinePattern_single(d=eff, orderByName=contr.names, numOfElement=nElement,
                    tissueName="Spleen", isotypeName="IgG", conds=sp.con)
      xt<-rbind(xt, temp)
        #Bone marrow IgM
      eff<-dt.clr.data.scaled[sp.IgM.cnames,]      
      temp<-dataLinePattern_single(d=eff, orderByName=contr.names, numOfElement=nElement,
                    tissueName="Spleen", isotypeName="IgM", conds=sp.con)
                    
     xt<-rbind(xt, temp)
     xt$geneIsotype<-paste0(xt$VGene, xt$isotype)
     trendline<-aggregate(xt$x, by=list(treatment=xt$treatment, isotype=xt$isotype, sign=xt$sign, tissue= xt$tissue), mean)        
trendline$geneIsotype=trendline$isotype
          
          pc1.line.all<-ggplot(data=xt, aes(x=treatment, y=x, group=geneIsotype,colour=isotype))+
            geom_line(linetype=2 )+
            geom_point()+
            geom_line(data=trendline, aes(y=x,x=treatment,colour=isotype), linetype=1,size=1.5)+
            facet_grid(sign~tissue)+theme_bw(base_size=15)+theme(axis.title.x=element_blank())
            
         
            tiff(file=here(output.dir,"pc1_lines_trend.tiff"), 
              width=800, height=650)
            pc1.line.all
            dev.off()
        
    #---------------------------------------------------
    ######     END                                                     ---
    #++++++++++++++++++++++++++++++++++++++
 ########################################
#              ----------------PC2  ---------------------  isotype effects     ##
#########################################
#now we start doing the Post hoc tests
pcs<-3
dt.mod.pc2<-data.frame(ilr=pca.VGene.clrd$x[,pcs], treatment=sample.conditions$treatment, isotype=sample.conditions$isotype, tissue=sample.conditions$tissue)
dt.mod.pc2$treatment<-factor(dt.mod.pc2$treatment, levels=c("PBS", "OVA", "OVA+CpG", "OVA+Alum"))
mod<-lm(ilr~treatment*isotype*tissue, data=dt.mod.pc2)
#em_tr<-emmeans(mod, ~tissue|treatment)
em_tr<-emmeans(mod, ~isotype|treatment*tissue)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
#    "OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), "OVA - PBS"=c(-1, 1,0,0)
    "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#    "Bone Marrow - Spleen" = c(1,-1)
  )

ec<-emmeans::contrast(em_tr, Set1, adjust='n')

em_tr<-emmeans(mod, ~treatment|isotype*tissue)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
    "OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), "OVA - PBS"=c(-1, 1,0,0)
#    "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#    "Bone Marrow - Spleen" = c(1,-1)
  )

ec.tr<-emmeans::contrast(em_tr, Set1, adjust='F')

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
    "OVA/CpG- OVA" = c(0,-1,1,0), "OVA/CpG - OVA/Alum"=c(0, 0,1,-1), "OVA/CpG - PBS"=c(-1, 0,1,0)
#    "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#    "Bone Marrow - Spleen" = c(1,-1)
  )

ec.tr.OVA<-emmeans::contrast(em_tr, Set1, adjust='F')



em_tr<-emmeans(mod, ~tissue|treatment*isotype)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
    #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), "OVA - PBS"=c(-1, 1,0,0)
   #"IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
 "Bone Marrow - Spleen" = c(1,-1)
  )

ec.tis<-emmeans::contrast(em_tr, Set1, adjust='n')

tr.pc2<-ggplot(data=dt.mod.pc2, aes(y=ilr, x=tissue, color=treatment))+
        #geom_point(aes(shape=treatment))+#lims(x=c(0.00, 0.034),y=c(0,0.04))+theme(legend.position = c(0.79, 0.18))+
        geom_boxplot(aes(shape=treatment), notch=F, notchwidth=0.1)+
        #lims(x=c(0.00, 0.034),y=c(0,0.04))+
        guides(shape="none")+
        labs( y="PC3 (IGHV usage)")+#, title=name.array[i]) +
       facet_grid(.~isotype)+
       geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Spleen"),
                               aes(x=1.75, xend=2.1, y=10.5, yend=10.5),color="black", size=0.9) + 
                     #geom_point(data=data.frame(treatment=c(2.2, 2.3,2.5,2.6, 1.3 ),ilr=c(0.47,0.47, 0.83,0.83, 0.2), isotype="IgG", tissue="Bone Marrow"),
                     #           color="black", size=2, shape=8) +   
        geom_text(data=data.frame(treatment=1.9,ilr=c(11.0 ), isotype=c("IgG"), tissue=rep("Spleen",1)),
                                color="black", size=4, aes(label=c("p<0.0005")))+
         geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Spleen"),
                               aes(x=1.9, xend=2.1, y=11.7, yend=11.7),color="black", size=0.9) + 
                     #geom_point(data=data.frame(treatment=c(2.2, 2.3,2.5,2.6, 1.3 ),ilr=c(0.47,0.47, 0.83,0.83, 0.2), isotype="IgG", tissue="Bone Marrow"),
                     #           color="black", size=2, shape=8) +   
        geom_text(data=data.frame(treatment=2.0,ilr=c(12.2 ), isotype=c("IgG"), tissue=rep("Spleen",1)),
                                color="black", size=4, aes(label=c("p=0.07")))+                        
        
        theme(legend.position = c(0.89, 0.75),text= element_text(size=18), axis.title.x=element_blank(),
                                    axis.text.x = element_text(angle = 0))
                                    
    save(file=here(data.dir,"figure4_pc2_draw.RData"), 
        dt.mod.pc2, tr.pc2)


#now plot the contributions
#two kinds of contribution
# PC loadings/rotation matrix
library(ggrepel)

#contributions of different variables to PC1, square the egein vector (without corrected by eigen value)
contr<-abs(pca.VGene.clrd$rotation[,pcs])*pca.VGene.clrd$rotation[,pcs]/sum(pca.VGene.clrd$rotation[,pcs]^2)
contr<-sort(contr)
contr.names<-names(contr[order(contr)])
rotation<-pca.VGene.clrd$rotation[,pcs]
save(file=here(data.dir,"gene_contribution_ordered_PC2.RData"), 
    contr, rotation)

dt.contr<-data.frame(contr=contr, IGHVs=contr.names, index=seq(1,length(contr)))
dt.contr.pos<-dt.contr[dt.contr$contr< -0.015 | dt.contr$contr> 0.015,]

contr.LD2<-ggplot(data=dt.contr, aes(x=index,y=contr))+
    geom_point()+scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))+ theme_bw(base_size=15)+
    geom_segment(x=-10, xend=140, y=0, yend=0, colour="blue", linetype=2,size=0.65 )+
    geom_segment(x=-10, xend=140, y=0.02, yend=0.02, colour="blue", linetype=3,size=0.4 )+
    geom_segment(x=-10, xend=140, y=-0.02, yend=-0.02, colour="blue", linetype=3,size=0.4 )+
    geom_text_repel(data=filter(dt.contr, contr< -0.020 | contr >0.020), aes(label=IGHVs),max.overlaps = Inf)+
    ggtitle(paste0("IGHV contribution to PC", pcs))+ xlab("")+ ylab("% contribution")+
    theme(plot.title = element_text(margin = margin(b = -20),size=15))

            ###contribution of variable to each PC
            #ig<-c("V9-2" , "V1-76", "V11-2" ,      "V12-3")
            ig<-c("V1-34", "V5-6" , "V5-9.1" ,      "V1-5")

            contr<-(pca.VGene.clrd$rotation[ig,])*matrix(rep(pca.VGene.clrd$sde, length(ig)), ncol=length(pca.VGene.clrd$sde), nrow=length(ig),byrow=T)#/sum(pca.VGene.clrd$rotation[,pcs]^2)

            contr<-abs(contr)*contr/apply(contr^2, 1, sum)
            #contr<-sort(abs(contr))
            #contr.names<-names(contr[order(contr)])
            cb<-vector("list", length=length(ig))
            for(i in 1:length(ig)){
              dt.contr<-data.frame(contr=contr[i,], IGHVs=colnames(contr), index=seq(1,dim(contr)[2]))
              #dt.contr.pos<-dt.contr[dt.contr$contr< -0.015 | dt.contr$contr> 0.015,]

              cb[[i]]<-ggplot(data=dt.contr, aes(x=index,y=contr))+
                  geom_bar(stat="identity", fill="turquoise2")+#ylim(-0.85, 0.1)+
                  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))+
                  #geom_segment(x=-10, xend=140, y=0.02, yend=0.02, colour="blue", linetype=2,size=0.35 )+
                  geom_segment(x=-10, xend=140, y=-0.00, yend=-0.0, colour="blue", linetype=2,size=0.35 )+theme_bw(base_size=15)+
                  #geom_segment(x=-10, xend=140, y=-0.02, yend=-0.02, colour="blue", linetype=3,size=0.4 )+
                  geom_text_repel(data=filter(dt.contr, contr< -0.25 | contr >0.25), aes(label=IGHVs),max.overlaps = Inf,size=5)+
                  ggtitle(paste0(ig[i]," contribution to PC", pcs))+ xlab("PCs")+ ylab("contribution (%)")+
                      
                  theme(plot.title = element_text(margin = margin(b = 0), hjust=0.5))
                  #cb<-c(cb, temp)
            }


tiff(here(output.dir,"PC2contributions_figure4.tiff"), 
    width=900, height=630)
ggarrange(contr.LD2, 
                ggarrange(cb[[1]],cb[[2]],cb[[3]], cb[[4]],  
                          ncol = 2, nrow = 2,
                          heights = c(2, 2)),  
          ncol = 1, nrow = 2, labels=c("A","B"),
          heights = c(2, 3))
dev.off()

 
    ###################################
    ##                                   plot lines                    ###
    ##____________________________###
        #original.names<-contr.names
      nElement<-9;
     #contr.names<-original.names
      #contr.names<-c("V5-16","V5-16")
      #dt.clr.data.scaled<-dt.clr.data
      eff<-dt.clr.data.scaled[sp.IgG.cnames,]
      xt.pc2<-dataLinePattern_single(d=eff, orderByName=contr.names, numOfElement=nElement,
                    tissueName="Spleen", isotypeName="IgG", conds=sp.con)
        #Bone marrow IgM
      eff<-dt.clr.data.scaled[sp.IgM.cnames,]      
      temp<-dataLinePattern_single(d=eff, orderByName=contr.names, numOfElement=nElement,
                    tissueName="Spleen", isotypeName="IgM", conds=sp.con)
                    
         xt.pc2<-rbind(xt.pc2,temp)
             xt.pc2$geneIsotype<-paste0(xt.pc2$VGene, xt.pc2$isotype)
             
             trendline.pc2<-aggregate(xt.pc2$x, by=list(treatment=xt.pc2$treatment, isotype=xt.pc2$isotype, 
                                                    sign=xt.pc2$sign, tissue= xt.pc2$tissue), mean)        
            trendline.pc2$geneIsotype=trendline.pc2$isotype
        #plot 
        pc2.line<-ggplot(data=xt.pc2, aes(x=treatment, y=x, group=geneIsotype,colour=isotype))+
            geom_line(linetype=2)+
            geom_point()+
            geom_line(data=trendline.pc2, aes(y=x,x=treatment,colour=isotype), linetype=1,size=1.5)+
            facet_grid(sign~.)+theme_bw(base_size=15)+theme(axis.title.x=element_blank())
     
     
     #       pdf("pc4_lines_dlScaledLog2.pdf")
            pc2.line
save(file=here(data.dir,"figure4_pc2_lines.RData"), 
    pc2.line, xt.pc2, trendline.pc2)
      #      dev.off()



########################################
#              ----------------PC3  ---------------------  actuall PC2     ##
#########################################
#now we start doing the Post hoc tests
pcs<-2
dt.mod<-data.frame(ilr=pca.VGene.clrd$x[,pcs], treatment=sample.conditions$treatment, isotype=sample.conditions$isotype, tissue=sample.conditions$tissue)
dt.mod$treatment<-factor(dt.mod$treatment, levels=c("PBS","OVA","OVA+CpG", "OVA+Alum"))
dt.mod$isotype<-factor(dt.mod$isotype, levels=c("IgM", "IgG"))
mod<-lm(ilr~treatment*isotype*tissue, data=dt.mod)
#em_tr<-emmeans(mod, ~tissue|treatment)
em_tr<-emmeans(mod, ~treatment|tissue*isotype)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
    "OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), "OVA - PBS"=c(-1, 1,0,0)
  #  "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#    "Bone Marrow - Spleen" = c(1,-1)
  )

ec<-emmeans::contrast(em_tr, Set1, adjust='FDR')

em_tr<-emmeans(mod, ~isotype)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
    #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), "OVA - PBS"=c(-1, 1,0,0)
   "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#    "Bone Marrow - Spleen" = c(1,-1)
  )

ec.iso<-emmeans::contrast(em_tr, Set1, adjust='FDR')

em_tr<-emmeans(mod, ~isotype|tissue*treatment)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
    #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), "OVA - PBS"=c(-1, 1,0,0)
   "IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
#    "Bone Marrow - Spleen" = c(1,-1)
  )

ec.iso2<-emmeans::contrast(em_tr, Set1, adjust='FDR')

em_tr<-emmeans(mod, ~tissue|isotype*treatment)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
    #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), "OVA - PBS"=c(-1, 1,0,0)
   #"IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
   "Bone Marrow - Spleen" = c(1,-1)
  )

ec.tissue<-emmeans::contrast(em_tr, Set1, adjust='FDR')

em_tr<-emmeans(mod, ~tissue|isotype)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
   #             #"OVA/PorB - PBS"=c(-1,0,1,0,0), 
   #             "OVA - PBS"=c(-1, 1,0,0)
    #"OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), "OVA - PBS"=c(-1, 1,0,0)
   #"IgM - IgG" = c(-1,1)#, "OVA/CpG - PBS"=c(-1,1)
   "Bone Marrow - Spleen" = c(1,-1)
  )

ec.tissue2<-emmeans::contrast(em_tr, Set1, adjust='FDR')


#plot the difference
iso.pc3<-ggplot(data=dt.mod, aes(y=ilr, x=isotype,  color=isotype))+
        #geom_boxplot(aes(shape=isotype), notch=F, notchwidth=0.1, position=position_dodge(0.8))+
        geom_dotplot( binaxis='y', position=position_dodge(0.8),stackdir='center',# colour="grey",
              stackratio=1.0, dotsize=0.5)+#lims(x=c(0.00, 0.034),y=c(0,0.04))+theme(legend.position = c(0.79, 0.18))+
               #lims(x=c(0.00, 0.034),y=c(0,0.04))+
        guides(shape="none")+
        labs(y="PC2 (IGHV usage)")+#, title=name.array[i]) 
        #facet_grid(.~tissue)+
         geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Bone Marrow"),
                    aes(x=1, xend=2, y=10.1, yend=10.1),color="black", size=0.9) +   
         #geom_point(data=data.frame(treatment=c(2.2, 2.3,2.5,2.6, 1.3 ),ilr=c(0.47,0.47, 0.83,0.83, 0.2), isotype="IgG", tissue="Bone Marrow"),
         #           color="black", size=2, shape=8) +   
         geom_text(data=data.frame(treatment=c(0.8),ilr=c(10.5), isotype=c(1.5), tissue=c(0.8)),
                    color="black", size=4, aes(label=c( "p<0.0001")))+
        theme(legend.position = c(0.89, 0.68),text= element_text(size=18), axis.title.x=element_blank())
tiff(file=here(output.dir,"figureS4_PC3_isotype.tiff"))
iso.pc3
dev.off()

iso2.pc3<-ggplot(data=dt.mod, aes(y=ilr, x=treatment,  color=isotype))+
        geom_boxplot(aes(shape=isotype), notch=F, notchwidth=0.1, position=position_dodge(0.8))+
        #geom_dotplot( binaxis='y', position=position_dodge(0.8),stackdir='center',# colour="grey",
        #      stackratio=1.0, dotsize=0.5)+#lims(x=c(0.00, 0.034),y=c(0,0.04))+theme(legend.position = c(0.79, 0.18))+
               lims(y=c(-5.5,11))+
        guides(shape="none")+
        labs(y="PC2 (IGHV usage)")+#, title=name.array[i]) 
        facet_grid(.~tissue)+
         geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Bone Marrow"),
                    aes(x=0.8, xend=1.2, y=2.0, yend=2.0),color="black", size=0.9) +   
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Bone Marrow"),
                    aes(x=1.8, xend=2.2, y=2.6, yend=2.6),color="black", size=0.9) +   
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Bone Marrow"),
                    aes(x=2.8, xend=3.2, y=4.8, yend=4.8),color="black", size=0.9) +   
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Bone Marrow"),
                    aes(x=3.8, xend=4.2, y=4.9, yend=4.9),color="black", size=0.9) +               
         geom_text(data=data.frame(treatment=c(1,2,3,4),ilr=c(2.2,2.8,5,5.1), isotype=c(1.5), tissue="Bone Marrow"),
                    color="black", size=6, aes(label=c( "*","***","****","****")))+
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Spleen"),
                    aes(x=0.8, xend=1.2, y=9.8, yend=9.8),color="black", size=0.9) +   
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Spleen"),
                    aes(x=1.8, xend=2.2, y=3.5, yend=3.5),color="black", size=0.9) +   
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Spleen"),
                    aes(x=2.8, xend=3.2, y=8, yend=8),color="black", size=0.9) +   
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Spleen"),
                    aes(x=3.8, xend=4.2, y=2.5, yend=2.5),color="black", size=0.9) +               
         geom_text(data=data.frame(treatment=c(1,2,3,4),ilr=c(10,3.7,8.2,2.7), isotype=c(1.5), tissue="Spleen"),
                    color="black", size=6, aes(label=c( "****","***","****","***")))+
        theme(legend.position = c(0.85, 0.11),text= element_text(size=18), axis.title.x=element_blank())
tiff(file=here(output.dir,"pc3_iso_figure4.tiff"), 
    width=600, height=500)
    iso2.pc3
dev.off()

tr.pc3<-ggplot(data=dt.mod, aes(y=ilr, x=isotype, color=treatment))+
        #geom_point(aes(shape=treatment))+#lims(x=c(0.00, 0.034),y=c(0,0.04))+theme(legend.position = c(0.79, 0.18))+
        geom_boxplot(aes(shape=treatment), notch=F, notchwidth=0.1)+
        #lims(x=c(0.00, 0.034),y=c(0,0.04))+
        guides(shape="none")+
        labs( y="PC2 (IGHV usage)")+#, title=name.array[i]) +
       facet_grid(.~tissue)+
       geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Bone Marrow"),
                    aes(x=1.7, xend=2.1, y=5, yend=5),color="black", size=0.9) +   
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Bone Marrow"),
                    aes(x=1.7, xend=2.3, y=6, yend=6),color="black", size=0.9) +
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Bone Marrow"),
                    aes(x=1.9, xend=2.1, y=7, yend=7),color="black", size=0.9) +   
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Bone Marrow"),
                    aes(x=1.9, xend=2.3, y=8, yend=8),color="black", size=0.9)+
         geom_text(data=data.frame(treatment=c(0.8),ilr=c(5.2,6.2,7.2,8.2), isotype=c(1.9,2.0,2.0,2.1), tissue=c("Bone Marrow")),
                    color="black", size=6, aes(label=c( "****","****","**","**")))+
                    
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Spleen"),
                    aes(x=1.7, xend=1.9, y=10, yend=10),color="black", size=0.9) +   
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Spleen"),
                    aes(x=1.7, xend=2.3, y=11, yend=11),color="black", size=0.9) +
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Spleen"),
                    aes(x=1.9, xend=2.1, y=12, yend=12),color="black", size=0.9) +   
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Spleen"),
                    aes(x=2.1, xend=2.3, y=13, yend=13),color="black", size=0.9)+
         geom_text(data=data.frame(treatment=c(0.8),ilr=c(10.2,11.2,12.2,13.2), isotype=c(1.8,2.0,2.0,2.2), tissue=c("Spleen")),
                    color="black", size=6, aes(label=c( "****","*****","***","***")))+
        theme(legend.position = c(0.12, 0.7),
                    text= element_text(size=18), axis.title.x=element_blank())
       
tiff(file=here(output.dir,"figureS4_PC3_treatment.tiff"), 
    width=600, height=500)
tr.pc3
dev.off()

#now plot the contributions
#two kinds of contribution
# PC loadings/rotation matrix


#contributions of different variables to PC1, square the egein vector (without corrected by eigen value)
contr<-abs(pca.VGene.clrd$rotation[,pcs])*pca.VGene.clrd$rotation[,pcs]/sum(pca.VGene.clrd$rotation[,pcs]^2)
contr<-sort(contr)
contr.names<-names(contr)
rotation<-pca.VGene.clrd$rotation[,pcs]
save(file=here(data.dir,"gene_contribution_ordered_PC3.RData"), 
    contr,rotation)

dt.contr<-data.frame(contr=contr, IGHVs=contr.names, index=seq(1,length(contr)))
dt.contr.pos<-dt.contr[dt.contr$contr< -0.015 | dt.contr$contr> 0.015,]

contr.LD4<-ggplot(data=dt.contr, aes(x=index,y=contr))+
    geom_point()+scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))+ theme_bw(base_size=15)+
    geom_segment(x=-10, xend=140, y=0, yend=0, colour="blue", linetype=2,size=0.65 )+
    geom_segment(x=-10, xend=140, y=0.02, yend=0.02, colour="blue", linetype=3,size=0.4 )+
    geom_segment(x=-10, xend=140, y=-0.02, yend=-0.02, colour="blue", linetype=3,size=0.4 )+
    geom_text_repel(data=filter(dt.contr, contr< -0.020 | contr >0.020), aes(label=IGHVs),max.overlaps = Inf)+
    ggtitle(paste0("IGHV contribution to PC",pcs))+ xlab("")+ ylab("% contribution")+
    theme(plot.title = element_text(margin = margin(b = -20),size=15))
 
###contribution of variable to each PC
ig<-c("V11-2","V6-6" ,  "V1-76" ,      "V9-2")

contr<-(pca.VGene.clrd$rotation[ig,])*matrix(rep(pca.VGene.clrd$sde, length(ig)), ncol=length(pca.VGene.clrd$sde), nrow=length(ig),byrow=T)#/sum(pca.VGene.clrd$rotation[,pcs]^2)

contr<-abs(contr)*contr/apply(contr^2, 1, sum)
#contr<-sort(abs(contr))
#contr.names<-names(contr[order(contr)])
cb<-vector("list", length=length(ig))
for(i in 1:length(ig))
{
  dt.contr<-data.frame(contr=contr[i,], IGHVs=colnames(contr), index=seq(1,dim(contr)[2]))
  #dt.contr.pos<-dt.contr[dt.contr$contr< -0.015 | dt.contr$contr> 0.015,]

  cb[[i]]<-ggplot(data=dt.contr, aes(x=index,y=contr))+
      geom_bar(stat="identity", fill="turquoise2")+#ylim(-0.85, 0.1)+
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))+
      #geom_segment(x=-10, xend=140, y=0.02, yend=0.02, colour="blue", linetype=2,size=0.35 )+
      geom_segment(x=-10, xend=140, y=-0.00, yend=-0.0, colour="blue", linetype=2,size=0.35 )+theme_bw(base_size=15)+
      #geom_segment(x=-10, xend=140, y=-0.02, yend=-0.02, colour="blue", linetype=3,size=0.4 )+
      geom_text_repel(data=filter(dt.contr, contr< -0.25 | contr >0.25), aes(label=IGHVs),max.overlaps = Inf,size=5)+
      ggtitle(paste0(ig[i]," contribution to PC", pcs))+ xlab("PCs")+ ylab("contribution (%)")+
          
      theme(plot.title = element_text(margin = margin(b = 0), hjust=0.5))
    #cb<-c(cb, temp)
}

tiff(file=here(output.dir,"PC3contributions_figure4.tiff"), 
    width=800, height=600)
ggarrange(contr.LD4, 
                ggarrange(cb[[1]],cb[[2]],cb[[3]], cb[[4]],  
                          ncol = 2, nrow = 2,
                          heights = c(2, 2)),  
          ncol = 1, nrow = 2, labels=c("A","B"),
          heights = c(2, 3))
dev.off()



    ###################################
    ##                                   plot lines                    ###
    ##____________________________###
     
     nElement<-10;
     #contr.names<-original.names
      #contr.names<-c("V5-16","V5-16")
      #dt.clr.data.scaled<-dt.clr.data
      eff<-dt.clr.data.scaled[sp.IgG.cnames,]
      xt<-dataLinePattern_single(d=eff, orderByName=contr.names, numOfElement=nElement,
                    tissueName="Spleen", isotypeName="IgG", conds=sp.con)
        #Bone marrow IgM
      eff<-dt.clr.data.scaled[sp.IgM.cnames,]      
      temp<-dataLinePattern_single(d=eff, orderByName=contr.names, numOfElement=nElement,
                    tissueName="Spleen", isotypeName="IgM", conds=sp.con)
                    
         xt<-rbind(xt,temp)
         
          eff<-dt.clr.data.scaled[bm.IgG.cnames,]
      temp<-dataLinePattern_single(d=eff, orderByName=contr.names, numOfElement=nElement,
                    tissueName="Bone Marrow", isotypeName="IgG", conds=bm.con)
      xt<-rbind(xt, temp)
        #Bone marrow IgM
      eff<-dt.clr.data.scaled[bm.IgM.cnames,]      
      temp<-dataLinePattern_single(d=eff, orderByName=contr.names, numOfElement=nElement,
                    tissueName="Bone Marrow", isotypeName="IgM", conds=bm.con)
                    
     xt<-rbind(xt, temp)
#     xt$geneIsotype<-paste0(xt$VGene, xt$isotype)
     trendline<-aggregate(xt$x, by=list( isotype=xt$isotype, sign=xt$sign), mean)        
trendline$geneIsotype=trendline$sign
#         trendline$VGene=trendline$isotype
         
             xt$geneIsotype<-paste0(xt$VGene,xt$tissue, xt$treatment)
            
        #plot 
        xt.s<-xt[xt$treatment=="PBS",]
        pc3.line<-ggplot(data=xt.s, aes(x=isotype, y=x, group=geneIsotype))+
            geom_line(linetype=2)+
            geom_point()+
            geom_line(data=trendline, aes(y=x,x=isotype), linetype=1,size=1.5,colour="blue")+
            facet_grid(sign~tissue)+
            theme_bw(base_size=15)+theme(axis.title.x=element_blank())
     
     
            tiff(file=here(output.dir,"pc3_lines_isotype.tiff"))
            pc3.line
           dev.off()
           
           
           
           




        ##################################
        ####now do the plotting for output   ###
        ###################################
#figure 4 in the main text  <-- now this is obsoleted, 
#   please go to the R code, figure4_plotting_v1.0.R

tiff(here(output.dir,"figure4_2.tiff"), 
    width=1200, height=2000)
ggarrange(
#level 1 for 
     ggarrange(  cn, g,  g2,
                          ncol = 3, nrow = 1, labels=c("A","B","C")
                          ),  
#level 2 for PC1 treatment
     ggarrange(  pc1.tr,  pc1.line,
                            #ggarrange(effect.ind.pc1[[3]], effect.ind.pc1[[1]],
                            #        ncol = 1, nrow = 2, labels=c("E","F")
                            #    ),
                         ncol = 2, nrow = 1, labels=c("D","E")
                         ),  
#level 3 for PC1 tissue
#     ggarrange(  pc1.ti, 
#                                ggarrange(effect.ind.pc1_2[[1]], effect.ind.pc1_2[[3]],
#                                       ncol = 1, nrow = 2, labels=c("H","I")
#                                ),
#                          ncol = 2, nrow = 1,  labels=c("G")
#                          ),                            
#level 4 for PC2 tissue
     ggarrange( tr.pc2, pc2.line,
                            #ggarrange(effect.ind.pc2[[1]], effect.ind.pc2[[3]],
                            #     ncol = 1, nrow = 2, labels=c("H","I")
                            #),
                          ncol = 2, nrow = 1,  labels=c("F","G")
                          ),     
#level 5 for PC4 tissue
#  ggarrange(  treat.pc4,  pc4.line,
                            #ggarrange(effect.ind.tr.pc4[[4]], effect.ind.tr.pc4[[2]],
                            #       ncol = 1, nrow = 2, labels=c("K","L")
                            #),
#                          ncol = 2, nrow = 1, labels=c("H","I")
 #                         ),                           
    ncol = 1, nrow = 3, heights=c(2,3,3)#, labels=c("A","D","F","H","J"),
    )
dev.off()




