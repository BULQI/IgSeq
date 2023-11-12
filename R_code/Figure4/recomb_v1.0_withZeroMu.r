#####R modules to analyze the recombinations events.
#    For Figure 4 in manuscript
## Make sure you have read ./ReadMe.txt for instructions
# of running scripts and preparing data.
##=======================

library(ggplot2)
library(here)
library(car)
library(ggpubr)
library(emmeans)


data.dir<-"Data/Figure4"
output.dir<-"R_code/Figure4"


#read the data from previously saved
load(here(data.dir,"All.RecSum.df.RData"))# list of  lists All.RecSum.bm, All.RecSum.sp; each of them are list of data frames
                                                            # All.RecSum.bm[["IgG"]], All.RecSum.bm[["IgM"]], All.RecSum.sp[["IgM"]] and All.RecSum.sp[["IgG"]]
                                                            # data frames 
                                                            #note this is the raw data conditions all the records/sequences.
load(here(data.dir,"Ig.good.RecSum.df.RData"))

isotype<-c("IgG", "IgM")
tissue<-c("Bone Marrow", "Spleen")


#now we need to get by sample.
#plot by samples, 
conditions$treatment<-factor(conditions$treatment, levels=c("PBS", "OVA", "OVA+PorB", "OVA+CpG", "OVA+Alum"))

########################now start doing the plotting of CDR3 length
rs<-IG.RecSum.bm[["IgG"]][, c("UID", "CDR3","CDR3Length","sampleName","MuFreq")]
rs$isotype<-"IgG"
temp<-IG.RecSum.bm[["IgM"]][, c("UID", "CDR3","CDR3Length","sampleName","MuFreq")]
temp$isotype<-"IgM"
rs<-rbind(rs, temp)
rs$tissue<-"Bone Marrow"
temp<-IG.RecSum.sp[["IgM"]][, c("UID", "CDR3","CDR3Length","sampleName","MuFreq")]
temp$isotype<-"IgM"
temp$tissue<-"Spleen"
rs<-rbind(rs, temp)

temp<-IG.RecSum.sp[["IgG"]][, c("UID", "CDR3","CDR3Length","sampleName","MuFreq")]
temp$isotype<-"IgG"
temp$tissue<-"Spleen"
rs<-rbind(rs, temp)
rs$treatment<-conditions[rs$sampleName, "treatment"]

###get rid of porb
rs<-rs[rs$treatment!="OVA+PorB",]
rs<-rs[rs$CDR3Length<100,]


#start plotting
 isotype<-unique(rs$isotype)
tissue<-unique(rs$tissue)
treatment<-unique(rs$treatment)

index<-0
plot(c(min(rs$CDR3Length), 100)#max(rs$CDR3Length))
                            , c(0,0.4) , type="n")
for(i in isotype)
{
    for(j in tissue)
    {
        for(k in treatment)
        {
            rs.temp<-rs[rs$isotype==i&rs$tissue==j& rs$treatment==k,]
            s<-unique(rs.temp$sampleName)
            index<-index+1
            for(m in s)
            {
                lines(density(rs.temp[rs.temp$sampleName==m,"CDR3Length"]), col=index)
            }
        }
    }
}

#get summary of CDR3 length
cdr3.mean<-aggregate(rs$CDR3Length, by=list(rs$tissue, rs$isotype, rs$treatment, rs$sampleName), mean)
#get summary of CDR3 length
cdr3.sd<-aggregate(rs$CDR3Length, by=list(rs$tissue, rs$isotype, rs$treatment, rs$sampleName), sd)
colnames(cdr3.mean)<-c("tissue", "isotype", "treatment", "sampleName","CDR3Length")
colnames(cdr3.sd)<-c("tissue", "isotype", "treatment", "sampleName","CDR3Length")

#boxplot(x~Group.1+Group.2+Group.3, data=cdr3.mean )
cdr3.mean$tissue<-factor(cdr3.mean$tissue, levels=c("Spleen","Bone Marrow"))
cdr3.mean$isotype<-factor(cdr3.mean$isotype, levels=c("IgM","IgG"))
cdr3.mean$treatment<-factor(cdr3.mean$treatment, levels=c("PBS","OVA","OVA+CpG","OVA+Alum"))
cdr3.mean$sampleName<-factor(cdr3.mean$sampleName)
#setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/CDR3IgGSub/")
save(file=here(data.dir,"CDR3_length_data.RData"), cdr3.mean, cdr3.sd)


options(contrasts = c("contr.sum", "contr.poly"))
lmCDR3<-lm(CDR3Length~tissue*isotype*treatment, data=cdr3.mean)
Anova(lmCDR3, type=3)

#doing the  simplified model
lmCDR3.2<-lm(CDR3Length~isotype*treatment, data=cdr3.mean)
Anova(lmCDR3.2, type=3)

#doing the 3 way analysis 
#
#dt.mod<-data.frame(ilr=intraD$ilr, treatment=intraD$treatment, isotype=intraD$isotype, tissue=intraD$tissue)
#mod<-lm(ilr~treatment*isotype*tissue, data=dt.mod)
em_tr<-emmeans(lmCDR3, ~treatment|tissue*isotype)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,0,1), "OVA/CpG - PBS"=c(-1,0,0,1,0), "OVA/PorB - PBS"=c(-1,0,1,0,0), "OVA - PBS"=c(-1, 1,0,0,0)
    "OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), "OVA - PBS"=c(-1, 1,0,0)#, 
      #  "OVA/CpG - OVA "=c(0, -1,1,0), "OVA/Alum- OVA "=c(0, -1,0,1)
  )
emmeans::contrast (em_tr, Set1, adjust ='F')
ec<-emmeans::contrast(em_tr, Set1, adjust='F')

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,0,1), "OVA/CpG - PBS"=c(-1,0,0,1,0), "OVA/PorB - PBS"=c(-1,0,1,0,0), "OVA - PBS"=c(-1, 1,0,0,0)
    "OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), 
    "OVA/Alum - OVA" = c(0, -1,0,1), "OVA/CpG - OVA"=c(0,-1,1,0) 
      #  "OVA/CpG - OVA "=c(0, -1,1,0), "OVA/Alum- OVA "=c(0, -1,0,1)
  )
emmeans::contrast (em_tr, Set1, adjust ='F')
ec.imm<-emmeans::contrast(em_tr, Set1, adjust='F')


Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,0,1), "OVA/CpG - PBS"=c(-1,0,0,1,0), "OVA/PorB - PBS"=c(-1,0,1,0,0), "OVA - PBS"=c(-1, 1,0,0,0)
    "OVA/CpG - OVA" = c(0,-1,1,0), "OVA/CpG - PBS"=c(-1,0,1,0), 
    "OVA/CpG - OVA/Alum" = c(0, 0,1,-1) 
      #  "OVA/CpG - OVA "=c(0, -1,1,0), "OVA/Alum- OVA "=c(0, -1,0,1)
  )
emmeans::contrast (em_tr, Set1, adjust ='F')
ec.CpG<-emmeans::contrast(em_tr, Set1, adjust='F')

#for tissue
em_tr<-emmeans(lmCDR3, ~tissue|treatment*isotype)
#em_tr %>% test(joint=F)

Set1 <- list(
  "Bone Marrow - Spleen" = c(-1,1))#, "OVA/CpG - PBS"=c(-1,0,0,1,0), "OVA/PorB - PBS"=c(-1,0,1,0,0), "OVA - PBS"=c(-1, 1,0,0,0))
emmeans::contrast (em_tr, Set1, adjust ='none')
ec<-emmeans::contrast(em_tr, Set1, adjust='none')


#for isotype
em_tr<-emmeans(lmCDR3, ~isotype|treatment*tissue)
#em_tr %>% test(joint=F)

Set1 <- list(
  "IgG - IgM" = c(-1,1))#, "OVA/CpG - PBS"=c(-1,0,0,1,0), "OVA/PorB - PBS"=c(-1,0,1,0,0), "OVA - PBS"=c(-1, 1,0,0,0))
emmeans::contrast (em_tr, Set1, adjust ='none')
ec<-emmeans::contrast(em_tr, Set1, adjust='none')

png(file=here(output.dir, "Figure4_CDR3_treatment.png"),   #"CDR3_treatment_withZeroMu.png"), 
  width=800, height=550)
cdr3<-ggplot(cdr3.mean, aes(y=CDR3Length, x=isotype, col=treatment))+
         geom_boxplot(aes(shape=treatment), notch=F, notchwidth=0.1)+
        #lims(x=c(0.00, 0.034),y=c(0,0.04))+
        guides(shape="none",colour=guide_legend(title="Immunization"))+
        labs(x="Isotype", y="CDR3 Length")+#, title=name.array[i]) +
       facet_grid(.~tissue)+
       geom_segment(data=data.frame(treatment=0.0,CDR3Length=41, isotype="IgG", tissue="Spleen"),
                    aes(x=1.7, xend=2.1, y=42.75, yend=42.75), colour="black",size=0.9) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
         geom_segment(data=data.frame(treatment=0.0,CDR3Length=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=1.7, xend=2.1, y=37.5, yend=37.5),color="black", size=0.9) +   
            geom_segment(data=data.frame(treatment=0.0,CDR3Length=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=1.7, xend=2.3, y=38, yend=38),color="black", size=0.9)+
#            geom_segment(data=data.frame(treatment=0.0,CDR3Length=0.01, isotype="IgG", tissue="Bone Marrow"),
#                    aes(x=1.7, xend=1.9, y=36.5, yend=36.5),color="black", size=0.9)+
            geom_text(data=data.frame(treatment=c(2.65,3.1,1.7),CDR3Length=c(43, 37.75,38.25), isotype=c(1.9,1.9,2.0)
                                                    , tissue=c("Spleen","Bone Marrow", "Bone Marrow")),
                    color="black", size=4, aes(label=c("p=0.002","p=0.001", "p=0.001")))+
       geom_segment(data=data.frame(treatment=0.0,CDR3Length=41, isotype="IgG", tissue="Spleen"),
                    aes(x=1.9, xend=2.1, y=43.25, yend=43.25), colour="black",size=0.9) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
         geom_segment(data=data.frame(treatment=0.0,CDR3Length=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=1.9, xend=2.1, y=38.75, yend=38.75),color="black", size=0.9) +   
            geom_segment(data=data.frame(treatment=0.0,CDR3Length=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=1.9, xend=2.3, y=39.25, yend=39.25),color="black", size=0.9)+
         geom_segment(data=data.frame(treatment=0.0,CDR3Length=0.01, isotype="IgG", tissue="Spleen"),
                    aes(x=2.3, xend=2.1, y=43.75, yend=43.75),color="black", size=0.9)+
#            geom_segment(data=data.frame(treatment=0.0,CDR3Length=0.01, isotype="IgG", tissue="Bone Marrow"),
#                    aes(x=1.7, xend=1.9, y=36.5, yend=36.5),color="black", size=0.9)+
            geom_text(data=data.frame(treatment=c(2.65,3.1,1.7,1.9),CDR3Length=c(43.5, 39,39.5,44), isotype=c(2.0,1.9,2.0,2.2)
                                                    , tissue=c("Spleen","Bone Marrow", "Bone Marrow","Spleen")),
                    color="black", size=4, aes(label=c("p=0.004","p=0.03", "p=0.02", "p=0.002")))             
cdr3 
dev.off()

save(cdr3, cdr3.mean, file=here(data.dir,"CDR3_treatment_fig_withZeroMu.RData"))
#now doing the isotype
#m.3way<-aggregate(cdr3.mean$CDR3Length, by=list(treatment=cdr3.mean$treatment, isotype=cdr3.mean$isotype, tissue=cdr3.mean$tissue), mean)
#v.3way<-aggregate(cdr3.mean$CDR3Length, by=list(treatment=cdr3.mean$treatment, isotype=cdr3.mean$isotype, tissue=cdr3.mean$tissue), FUN=function(x){sd(x)/sqrt(length(x))})
#mv<-data.frame(mean=m.3way$x, var=v.3way$x, treatment=m.3way$treatment, isotype=m.3way$isotype, tissue=m.3way$tissue) 
png(file=here(output.dir,"CDR3_isotype_withZeroMu.png"), 
  width=800, height=550)
ggplot(data=cdr3.mean, aes(y=CDR3Length, x=treatment, col=isotype))+
         geom_boxplot(aes(shape=treatment), notch=F, notchwidth=0.1)+
        #lims(x=c(0.00, 0.034),y=c(0,0.04))+
        guides(shape="none")+
        labs(x="Treatment", y="CDR3 Length")+#, title=name.array[i]) +
        
       facet_grid(.~tissue)+
       geom_segment(data=data.frame(treatment=0.0,CDR3Length=42.5, isotype="IgG", tissue="Spleen"),
                    aes(x=2.8, xend=3.2, y=42.7, yend=42.7), colour="black",size=0.9) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
#         geom_segment(data=data.frame(treatment=0.0,CDR3Length=0.01, isotype="IgG", tissue="Bone Marrow"),
#                    aes(x=2.8, xend=3.2, y=37.5, yend=37.5),color="black", size=0.9) +   
            geom_segment(data=data.frame(treatment=0.0,CDR3Length=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=3.8, xend=4.2, y=37.5, yend=37.5),color="black", size=0.9)+
            geom_segment(data=data.frame(treatment=0.0,CDR3Length=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=0.8, xend=1.2, y=36.5, yend=36.5),color="black", size=0.9)+
            geom_text(data=data.frame(treatment=c(3,3,4.0),CDR3Length=c(43, 38,38), isotype=c(1.9,1.9,2.0)
                                                    , tissue=c("Spleen", "Bone Marrow", "Bone Marrow")),
                    color="black", size=4, aes(label=c("p=0.001","p=0.06", "p=0.05")))
dev.off()

png(file=here(output.dir,"CDR3_tissue_withZeroMu.png"), 
  width=800, height=550)
ggplot(data=cdr3.mean, aes(y=CDR3Length, x=treatment, col=tissue))+
         geom_boxplot(aes(shape=treatment), notch=F, notchwidth=0.1)+
        #lims(x=c(0.00, 0.034),y=c(0,0.04))+
        guides(shape="none")+
        labs(x="Treatment", y="CDR3 Length")+#, title=name.array[i]) +
        
       facet_grid(.~isotype)+
       geom_segment(data=data.frame(treatment=0.0,CDR3Length=42.5, isotype="IgG", tissue="Spleen"),
                    aes(x=3.8, xend=4.2, y=37.5, yend=37.5), colour="black",size=0.9) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
         
            geom_text(data=data.frame(treatment=c(4.0),CDR3Length=c(38), isotype=c("IgG")
                                                    , tissue=c("Bone Marrow")),
                    color="black", size=4, aes(label=c("p=0.05")))
dev.off()
