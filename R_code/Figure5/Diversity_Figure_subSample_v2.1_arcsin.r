#R code to estimate the repertoire size
#--- 
# Make sure you have read ./ReadMe.txt for instructions
# of running scripts and preparing data.
# 
#----------------------

library(SpadeR)
library(gridBase)
library(grid)
library(here)


library(ggplot2)
library(ggalt)
library(ggpubr)
library(car)
library(emmeans)


data.dir<-"Data/Figure5"
output.dir<-"R_code/Figure5"
#now the hardest part to generate the data that can be used by the packages, SpadeR and iNEXT
#see the user guides for the two packages.
#
#read the saved data, which were saved by DataAnalysis_v1.0_clone_data.R
#read the files 15 mice (5 group x 3 mice each and two tissue)

load(here(data.dir,"clones.Rdata")) # loaded with BM.clones, SP.clones and conditions, (list(IgG and IgM), list (IgG+IgM), data.frame)
                                               #clone.inc data , clone incidence freq data (not count data yet)

                    
load(file=here(data.dir,"diversity_subSample2.RData"))#generated in Diversity_subsample_RCMD_v2.0.r on linux.
                        #two pieces of data are loaded : ReadMe.text and div.tbl (this one is a data format table holding hill numbers 0,4)

#now add conditions/grouping information
conditions$ID<-as.character(conditions$ID)
rownames(conditions)<-conditions$ID
div.tbl$sampleName<-as.character(div.tbl$sampleName)

div.tbl<-cbind(div.tbl, conditions[div.tbl$sampleName,c("snum", "ID","treatment")])
div.tbl$tissue<-as.character(div.tbl$tissue)
div.tbl$treatment<-as.character(div.tbl$treatment)

###now save the diversity data for doing MFA and clustering


      #+++++++++++++++++++++++++++++++++++++++++++     
      ######do unnormalized data plotting first.                  ++
      #---------------------------------------------------------------------
      ####
      
div.t0<-div.tbl[div.tbl$q==0&div.tbl$treatment!="OVA+PorB",]
rownames(div.t0)<-paste0(div.t0$tissue, div.t0$isotype, div.t0$snum)

div.t3<-data.frame()
for(i in seq(0,4,0.25))
{
    temp<-div.tbl[div.tbl$q==i&div.tbl$treatment!="OVA+PorB",]
    rownames(temp)<-paste0(temp$tissue, temp$isotype, temp$snum)

    temp<-temp[rownames(div.t0),]

    temp$proportion<-temp$ChaoJost#/div.t0$ChaoJost
    div.t3<-rbind(div.t3,temp)
}
rownames(div.t3)<-paste0(div.t3$tissue, div.t3$isotype, div.t3$snum, "q",div.t3$q)
div.t3$trans<-log(div.t3$proportion)#-log(1-div.t3$proportion)##(asin(sqrt(div.t3$proportion)))#div.t3$proportion # #
div.t3$treatment=factor(div.t3$treatment, levels=c("PBS","OVA", "OVA+CpG","OVA+Alum"))
div.t3$isotype=factor(div.t3$isotype, levels=c("IgM","IgG")) 
div.t3$tissue=factor(div.t3$tissue, levels=c("Spleen","Bone Marrow")) 
m.div.t3<-aggregate(div.t3$trans, by=list(div.t3$isotype, div.t3$treatment, div.t3$tissue, div.t3$q),
                        FUN=mean)
                        

m<-aggregate(div.t3$trans, by=list(div.t3$q, div.t3$tissue, div.t3$isotype, div.t3$treatment), mean)
v<-aggregate(div.t3$trans, by=list(div.t3$q, div.t3$tissue, div.t3$isotype, div.t3$treatment), FUN=function(x){sd(x)/sqrt(length(x))})
names(m)<-c("q", "tissue", "isotype", "treatment", "mean")
names(v)<-c("q", "tissue", "isotype", "treatment", "std")
m<-merge(m,v)
div.m<-merge(x=div.t3, y=m)
#div.m$treatment<-factor(div.m$treatment, levels=c("PBS", "OVA", "OVA+PorB", "OVA+CpG", "OVA+Alum"))
#div.m$tissue<-factor(div.m$tissue, levels=c("Spleen", "Bone Marrow"))
div.m$ymax<-div.m$mean+div.m$std
div.m$ymin<-div.m$mean-div.m$std
div.m$ChaoJost_log<-div.m$trans
#div.m[div.m$ymin<4,"ymin"]<-4

#pdf("HillNumbers_errorbar_SP.pdf")
#temp<-div.m[div.m$tissue=="Bone Marrow"&div.m$treatment!="OVA+PorB",]
temp<-div.m[div.m$treatment=="OVA+Alum",]


div.0<-temp[temp$q==0,]
      div.0<-rbind(div.0, temp[temp$q==1,])
      div.0<-rbind(div.0, temp[temp$q==2,])
      div.0<-rbind(div.0, temp[temp$q==3,])
      div.0<-rbind(div.0, temp[temp$q==4,])                        
     
       iso.unnorm<-ggplot(data=temp, aes(x=q, y=ChaoJost_log, shape =(treatment),color=(treatment)))+
        geom_line(aes(x = q, y = mean, color=(treatment)), size = 1.2 , 
                            #position= position_dodge2(width=0.2), 
                            linetype=1#"dashed"#col = "red"
                        )+#scale_y_continuous(trans="log10")+
        geom_errorbar(data=div.0,aes(ymin=ymin, ymax=ymax, color=(treatment)),size=0.5, width=0.05
                    #,position= position_dodge2(width=0.2)
                    ) +
                         theme(text = element_text(size=12),legend.title = element_blank())+
      xlab("order (q)") + ylab("Diversity (HILL NUMBERS)") + 
      theme(plot.title = element_text(hjust = 0.5, size=12), legend.position = c(0.2,0.850),  
      strip.text.x = element_text(size = 11, colour = "blue"),strip.text.y = element_text(size = 11, colour = "blue"))+
      facet_grid(tissue~isotype )






##################################***********
##                         start doing normalized 
###########################################
div.t0<-div.tbl[div.tbl$q==0&div.tbl$treatment!="OVA+PorB",]
rownames(div.t0)<-paste0(div.t0$tissue, div.t0$isotype, div.t0$snum)

div.t3<-data.frame()
div.t3.der<-data.frame()
for(i in seq(0,4,0.25))
{
    
    temp<-div.tbl[div.tbl$q==i&div.tbl$treatment!="OVA+PorB",]
    rownames(temp)<-paste0(temp$tissue, temp$isotype, temp$snum)

    temp<-temp[rownames(div.t0),]

    temp$proportion<-temp$ChaoJost/div.t0$ChaoJost
    div.t3<-rbind(div.t3,temp)
    #doing derivative  <--- new code on march 14th. but don't use below
    if(i!=4)
    {
    cat("round i:", i,"\n")
    div.tNext<-div.tbl[div.tbl$q==i+0.25&div.tbl$treatment!="OVA+PorB",]
    temp.der<-div.tbl[div.tbl$q==i&div.tbl$treatment!="OVA+PorB",]
    rownames(temp.der)<-paste0(temp.der$tissue, temp$isotype, temp.der$snum)
    temp.der$derivative<-(temp.der$ChaoJost-div.tNext$ChaoJost)/temp.der$ChaoJost
    div.t3.der<-rbind(div.t3.der, temp.der)
    }
    
}
###==> div.t3
#div.t3<-div.t3.der
rownames(div.t3)<-paste0(div.t3$tissue, div.t3$isotype, div.t3$snum, "q",div.t3$q)
div.t3$trans<-(asin(sqrt(div.t3$proportion)))#log(div.t3$proportion)#-log(1-div.t3$proportion)##div.t3$proportion # #
#div.t3$trans<-log(1/div.t3$proportion)
#div.t3$trans<-log(div.t3$derivative)
m.div.t3<-aggregate(div.t3$trans, by=list(div.t3$isotype, div.t3$treatment, div.t3$tissue, div.t3$q),
                        FUN=mean)

#setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/figure5/subsample/")

save(file=here(data.dir,"Diversity_table_all_subsample_normArcsin.RData"), div.t3)



options(contrasts = c("contr.sum", "contr.poly"))  
div.t3$treatment=factor(div.t3$treatment, levels=c("PBS","OVA", "OVA+CpG","OVA+Alum"))
div.t3$isotype=factor(div.t3$isotype, levels=c("IgM","IgG")) 
div.t3$tissue=factor(div.t3$tissue, levels=c("Spleen","Bone Marrow")) 
x.lm<-lm(data=div.t3[div.t3$q==4,], trans~tissue*isotype*treatment)
Anova(x.lm, type=3)


#dt.mod<-data.frame(q3=div.qs[,"q3"], treatment=div.qs$treatment, isotype=div.qs$isotype, tissue=div.qs$tissue)

#mod<-lm(q3~treatment*isotype*tissue, data=dt.mod)
em_tr<-emmeans(x.lm, ~treatment|tissue*isotype)
#em_tr<-emmeans(mod, ~tissue|isotype*treatment)
#em_tr %>% test(joint=F)

#Set1 <- list(  "OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0),  "OVA - PBS"=c(-1, 1,0,0))
Set1 <- list(  "OVA/Alum - OVA" = c(0,-1,0,1), "OVA/CpG - OVA"=c(0,-1,1,0),  "OVA - PBS"=c(-1, 1,0,0),
                    "OVA/Alum - PBS"=c(-1,0,0,1),"OVA/CpG - PBS"=c(-1,0,1,0),"OVA/Alum - OVA/CpG"=c(0,0,-1,1))
#Set1 <- list("Spleen - Bone Marrow" = c(-1,1))#, "OVA/CpG - PBS"=c(-1,0,1,0),  "OVA - PBS"=c(-1, 1,0,0))
#Set1 <- list("IgG - IgM" = c(-1,1))#, "OVA/CpG - PBS"=c(-1,0,1,0),  "OVA - PBS"=c(-1, 1,0,0))
ec<-emmeans::contrast (em_tr, Set1, adjust ='F')

#em_tr %>% test(joint=F)
em_tr<-emmeans(x.lm, ~isotype|tissue*treatment)
#em_tr<-emmeans(mod, ~tissue|isotype*treatment)
#em_tr %>% test(joint=F)

#Set1 <- list(  "OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0),  "OVA - PBS"=c(-1, 1,0,0))
#Set1 <- list("Spleen - Bone Marrow" = c(-1,1))#, "OVA/CpG - PBS"=c(-1,0,1,0),  "OVA - PBS"=c(-1, 1,0,0))
Set1 <- list("IgG - IgM" = c(-1,1))#, "OVA/CpG - PBS"=c(-1,0,1,0),  "OVA - PBS"=c(-1, 1,0,0))
ec2<-emmeans::contrast (em_tr, Set1, adjust ='F')

em_tr<-emmeans(x.lm, ~tissue|isotype*treatment)
#em_tr<-emmeans(mod, ~tissue|isotype*treatment)
#Set1 <- list(  "OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0),  "OVA - PBS"=c(-1, 1,0,0))
Set1 <- list("Bone Marrow - Spleen" = c(-1,1))#, "OVA/CpG - PBS"=c(-1,0,1,0),  "OVA - PBS"=c(-1, 1,0,0))
#Set1 <- list("IgG - IgM" = c(-1,1))#, "OVA/CpG - PBS"=c(-1,0,1,0),  "OVA - PBS"=c(-1, 1,0,0))
ec3<-emmeans::contrast (em_tr, Set1, adjust ='F')


m<-aggregate(div.t3$trans, by=list(div.t3$q, div.t3$tissue, div.t3$isotype, div.t3$treatment), mean)
v<-aggregate(div.t3$trans, by=list(div.t3$q, div.t3$tissue, div.t3$isotype, div.t3$treatment), FUN=function(x){sd(x)/sqrt(length(x))})
names(m)<-c("q", "tissue", "isotype", "treatment", "mean")
names(v)<-c("q", "tissue", "isotype", "treatment", "std")
m<-merge(m,v)
div.m<-merge(x=div.t3, y=m)
#div.m$treatment<-factor(div.m$treatment, levels=c("PBS", "OVA", "OVA+PorB", "OVA+CpG", "OVA+Alum"))
#div.m$tissue<-factor(div.m$tissue, levels=c("Spleen", "Bone Marrow"))
div.m$ymax<-div.m$mean+div.m$std
div.m$ymin<-div.m$mean-div.m$std
div.m$ChaoJost_log<-div.m$trans
#div.m[div.m$ymin<4,"ymin"]<-4

#pdf("HillNumbers_errorbar_SP.pdf")
#temp<-div.m[div.m$tissue=="Bone Marrow"&div.m$treatment!="OVA+PorB",]
temp<-div.m#[div.m$treatment!="OVA+PorB",]


div.0<-temp[temp$q==0,]
      div.0<-rbind(div.0, temp[temp$q==1,])
      div.0<-rbind(div.0, temp[temp$q==2,])
      div.0<-rbind(div.0, temp[temp$q==3,])
      div.0<-rbind(div.0, temp[temp$q==4,])
tr<-ggplot(data=temp, aes(x=q, y=ChaoJost_log, shape =(treatment),color=(treatment)))+
        geom_line(aes(x = q, y = mean, color=(treatment)), size = 1.2 , 
                            #position= position_dodge2(width=0.2), 
                            linetype=1#"dashed"#col = "red"
                        )+
        geom_errorbar(data=div.0, aes(ymin=ymin, ymax=ymax, color=(treatment)),size=0.5, width=0.05
                    #,position= position_dodge2(width=0.2)
                    ) +
        #            geom_point(size =3.5, position= position_dodge2(width=0.2)
        #                            ) + 
                        #scale_y_continuous(trans="log10")+ 
                                          #,breaks = trans_breaks("log10", 
                                          #         function(x) 10^x), 
                                          #labels = trans_format("log10", 
                                          #          math_format(10^.x))
                         #                 ) + 
                         theme(text = element_text(size=12),legend.title = element_blank())+
 #     annotate("text", label="Spleen", x=0, y=13.5, size=6,  hjust=0, colour="red") + 
      xlab("order (q)") + ylab("Normalized Diversity (HILL NUMBERS)") +
      theme(plot.title = element_text(hjust = 0.5, size=12), legend.position = c(0.4,0.93),  
      strip.text.x = element_text(size = 11, colour = "blue"),strip.text.y = element_text(size = 11, colour = "blue"))+
      facet_grid(tissue~isotype )+
      geom_text(data=data.frame(treatment=rep("PBS",4),
                            ChaoJost_log=rep(0.7,4),
                        isotype=c("IgM", "IgM","IgG","IgG"), 
                        q=rep(4.0,4),
                     tissue=c("Spleen", "Bone Marrow","Spleen","Bone Marrow")),
                    color="black", size=7, aes(label=c("","*","*","***")),hjust = rep(0.5,4))+
        geom_text(data=data.frame(treatment=rep("PBS",4),
                            ChaoJost_log=rep(0.85,4),
                        isotype=c("IgM", "IgM","IgG","IgG"), 
                        q=rep(3.0,4),
                     tissue=c("Spleen", "Bone Marrow","Spleen","Bone Marrow")),
                    color="black", size=7, aes(label=c("","**","","**")),hjust = rep(0.5,4))+
        geom_text(data=data.frame(treatment=rep("PBS",4),
                            ChaoJost_log=rep(0.95,4),
                        isotype=c("IgM", "IgM","IgG","IgG"), 
                        q=rep(2.0,4),
                     tissue=c("Spleen", "Bone Marrow","Spleen","Bone Marrow")),
                    color="black", size=7, aes(label=c("*","***","","*")),hjust = rep(0.5,4))+
       geom_text(data=data.frame(treatment=rep("PBS",4),
                            ChaoJost_log=rep(1.25,4),
                        isotype=c("IgM", "IgM","IgG","IgG"), 
                        q=rep(1.0,4),
                     tissue=c("Spleen", "Bone Marrow","Spleen","Bone Marrow")),
                    color="black", size=7, aes(label=c("","****","","")),hjust = rep(0.5,4))
 tiff(here(output.dir,"Diversity3way_tr_arc_supp.tiff"), 
  width=800, height=500)
tr
dev.off();

#### 
 iso<-ggplot(data=temp, aes(x=q, y=ChaoJost_log, shape =(isotype),color=(isotype)))+
        geom_line(aes(x = q, y = mean, color=(isotype)), size = 1.2 , 
                            #position= position_dodge2(width=0.2), 
                            linetype=1#"dashed"#col = "red"
                        )+#scale_y_continuous(trans="log10")+
        geom_errorbar(data=div.0,aes(ymin=ymin, ymax=ymax, color=(isotype)),size=0.5, width=0.05
                    #,position= position_dodge2(width=0.2)
                    ) +
                         theme(text = element_text(size=12),legend.title = element_blank())+
      xlab("order (q)") + ylab("Normalized Diversity (HILL NUMBERS)") + 
      theme(plot.title = element_text(hjust = 0.5, size=12), legend.position = c(0.2,0.750),  
      strip.text.x = element_text(size = 11, colour = "blue"),strip.text.y = element_text(size = 11, colour = "blue"))+
      facet_grid(tissue~treatment )+
      geom_text(data=data.frame(treatment=c("PBS", "PBS",
                                                                                                    "OVA","OVA","OVA","OVA",
                                                                                                        "OVA+CpG","OVA+CpG",
                                                        "OVA+Alum","OVA+Alum"),
                            ChaoJost_log=c(1.3, 0.7,
                                                                    1.3,0.7,0.6, 0.5,
                                                                    1.3, 0.7,
                                                                    1.3,0.7),
                        isotype="IgG", q=c(1.0,2.0, 
                                                                        1.0,2.0,3.0, 4.0,
                                                                        1, 4,
                                                                        1,4)
                    , tissue="Bone Marrow"),
                    color="black", size=7, aes(label=c("*","*",
                                                                "*","****","****","***",
                                                                    "*","*",
                                                                    "*","*")),hjust = rep(0.5,10))+
     geom_text(data=data.frame(treatment=c("PBS", "OVA","OVA+CpG","OVA+Alum","OVA+Alum"),
                            ChaoJost_log=c(1.3,1.3,1.3,1.3,1.0),
                        isotype="IgG", q=c(1.0,1.0,1.0,1.0,2.0)
                    , tissue="Spleen"),
                    color="black", size=7, aes(label=c("****","****","****","****","**")),hjust = rep(0.5,5))               
        
        tis<-ggplot(data=temp, aes(x=q, y=ChaoJost_log, shape =(tissue),color=(tissue)))+
        geom_line(aes(x = q, y = mean, color=(tissue)), size = 1.2 , 
                            #position= position_dodge2(width=0.2), 
                            linetype=1#"dashed"#col = "red"
                        )+
        geom_errorbar(data=div.0,aes(ymin=ymin, ymax=ymax, color=(tissue)),size=0.5, width=0.05
                    #,position= position_dodge2(width=0.2)
                    ) +
                         theme(text = element_text(size=12),legend.title = element_blank())+
      xlab("order (q)") + ylab("Normalized Diversity (HILL NUMBERS)") + 
      theme(plot.title = element_text(hjust = 0.5, size=12), legend.position = c(0.75,0.7),  
      strip.text.x = element_text(size = 11, colour = "blue"),strip.text.y = element_text(size = 11, colour = "blue"))+
      facet_grid(isotype~treatment )+
      geom_text(data=data.frame(treatment=c("PBS", "OVA","OVA","OVA","OVA"),
                            ChaoJost_log=c(1.2,1.2,0.8,0.6, 0.5),
                        isotype="IgM", q=c(1.0,1.0,2.0,3.0,4)
                    , tissue="Spleen"),
                    color="black", size=7, aes(label=c("****","***","***","**","*")),hjust = rep(0.5,5))   +
       geom_text(data=data.frame(treatment=c("PBS", "OVA","OVA+CpG"
                                                            ,"OVA+Alum","OVA+Alum","OVA+Alum","OVA+Alum"),
                            ChaoJost_log=c(1.0,1.0,1.0,1.0,0.8,0.6,  0.5),
                        isotype="IgG", q=c(1.0,1.0,1.0,1,2,3,4)
                    , tissue="Spleen"),
                    color="black", size=7, aes(label=c("*","**","*","**","**", "**","**")),hjust = rep(0.5,7))              
 tiff(here(output.dir,"Diversity3way_iso_arc_supp.tiff"), 
  width=800, height=500)
iso
dev.off();     
 tiff(here(output.dir,"Diversity3way_tis_arc_supp.tiff"), 
  width=800, height=500)
tis
dev.off();              

#partial with pbs
    temp<-div.m[div.m$treatment=="PBS",]

div.0<-temp[temp$q==0,]
      div.0<-rbind(div.0, temp[temp$q==1,])
      div.0<-rbind(div.0, temp[temp$q==2,])
      div.0<-rbind(div.0, temp[temp$q==3,])
      div.0<-rbind(div.0, temp[temp$q==4,])
    iso.pbs<-ggplot(data=temp, aes(x=q, y=ChaoJost_log, shape =(isotype),color=(isotype)))+
        geom_line(aes(x = q, y = mean, color=(isotype)), size = 1.2 , 
                            #position= position_dodge2(width=0.2), 
                            linetype=1#"dashed"#col = "red"
                        )+#scale_y_continuous(trans="log10")+
        geom_errorbar(data=div.0,aes(ymin=ymin, ymax=ymax, color=(isotype)),size=0.5, width=0.05
                    #,position= position_dodge2(width=0.2)
                    ) +
                         theme(text = element_text(size=12),legend.title = element_blank())+
      xlab("order (q)") + ylab("Normalized Diversity (HILL NUMBERS)") + 
      theme(plot.title = element_text(hjust = 0.5, size=12), legend.position = c(0.65,0.750),  
      strip.text.x = element_text(size = 11, colour = "blue"),strip.text.y = element_text(size = 11, colour = "blue"))+
      facet_grid(tissue~treatment )+
      geom_text(data=data.frame(treatment=c("PBS", "PBS"
                                                                                                    ),
                            ChaoJost_log=c(1.3, 0.7
                                                                    ),
                        isotype="IgG", q=c(1.0,2.0)
                    , tissue="Bone Marrow"),
                    color="black", size=7, aes(label=c("*","*")),hjust = rep(0.5,2))+
     geom_text(data=data.frame(treatment=c("PBS"),
                            ChaoJost_log=c(1.3),
                        isotype="IgG", q=c(1.0)
                    , tissue="Spleen"),
                    color="black", size=7, aes(label=c("****")),hjust = rep(0.5,1))    


        tis.pbs<-ggplot(data=temp, aes(x=q, y=ChaoJost_log, shape =(tissue),color=(tissue)))+
        geom_line(aes(x = q, y = mean, color=(tissue)), size = 1.2 , 
                            #position= position_dodge2(width=0.2), 
                            linetype=1#"dashed"#col = "red"
                        )+
        geom_errorbar(data=div.0,aes(ymin=ymin, ymax=ymax, color=(tissue)),size=0.5, width=0.05
                    #,position= position_dodge2(width=0.2)
                    ) +
                         theme(text = element_text(size=12),legend.title = element_blank())+
      xlab("order (q)") + ylab("Normalized Diversity (HILL NUMBERS)") + 
      theme(plot.title = element_text(hjust = 0.5, size=12), legend.position = c(0.75,0.7),  
      strip.text.x = element_text(size = 11, colour = "blue"),strip.text.y = element_text(size = 11, colour = "blue"))+
      facet_grid(isotype~treatment )+
      geom_text(data=data.frame(treatment=c("PBS"),
                            ChaoJost_log=c(1.2),
                        isotype="IgM", q=c(1.0)
                    , tissue="Spleen"),
                    color="black", size=7, aes(label=c("****")),hjust = rep(0.5,1))   +
       geom_text(data=data.frame(treatment=c("PBS"),
                            ChaoJost_log=c(1.0),
                        isotype="IgG", q=c(1.0)
                    , tissue="Spleen"),
                    color="black", size=7, aes(label=c("*")),hjust = rep(0.5,1)) 
###############3       
    #plot the comparison of treat 
    #####start doing the analysis
      #first get only q=c(0,1,2,3)
      
        dt.mod<-data.frame(q3=div.t3[div.t3$q==4,"trans"], treatment=div.t3[div.t3$q==4,"treatment"], 
                                isotype=div.t3[div.t3$q==4,"isotype"], tissue=div.t3[div.t3$q==4,"tissue"])
#mod<-lm(q3~treatment*isotype*tissue, data=dt.mod)
        #----start plotting the individual effects 
#get the mean , treatment
q3.tr<-ggplot(data=dt.mod, aes(y=q3, x=isotype, color=treatment))+
        #geom_point(aes(shape=treatment))+#lims(x=c(0.00, 0.034),y=c(0,0.04))+theme(legend.position = c(0.79, 0.18))+
        geom_boxplot(aes(shape=isotype), notch=F, notchwidth=0.1, varwidth=T)+
        #lims(x=c(0.00, 0.034),y=c(0,0.04))+
        guides(shape="none")+
        labs( y="Normalized Diversity(q=4)")+#, title=name.array[i]) +
       facet_grid(.~tissue)+ 
        geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=1.7, xend=1.9, y=0.39, yend=0.39), colour="black",size=0.9) +
         geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=1.7, xend=2.1, y=0.41, yend=0.41),color="black", size=0.9) +   
         geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                   aes(x=1.7, xend=2.3, y=0.55, yend=0.55),color="black", size=0.9) +   
          geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                   aes(x=1.9, xend=2.3, y=0.57, yend=0.57),color="black", size=0.9) +
         geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                   aes(x=2.1, xend=2.3, y=0.59, yend=0.59),color="black", size=0.9) +
                   
         geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgM", tissue="Bone Marrow"),
                    aes(x=0.9, xend=1.3, y=0.33, yend=0.33), colour="black",size=0.9) +
         geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgG", tissue="Spleen"),
                    aes(x=1.7, xend=1.9, y=0.42, yend=0.42),color="black", size=0.9) +
         #geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgG", tissue="Spleen"),
         #           aes(x=1.7, xend=1.9, y=3.8, yend=3.8), colour="black",size=0.9) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
         #geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgG", tissue="Spleen"),
         #           aes(x=1.7, xend=2.1, y=4.3, yend=4.3),color="black", size=0.9)+
         #geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgG", tissue="Spleen"),
         #           aes(x=1.7, xend=2.3, y=4.8, yend=4.8),color="black", size=0.9)+
         #geom_point(data=data.frame(treatment=c(2.2, 2.3,2.5,2.6, 1.3 ),ilr=c(0.47,0.47, 0.83,0.83, 0.2), isotype="IgG", tissue="Bone Marrow"),
         #           color="black", size=2, shape=8) +   
         geom_text(data=data.frame(treatment=c(0.9),q3=c(0.395, 0.415, 0.555,0.575,0.595), isotype=c(1.8,1.9,2,2.1,2.2)
                    , tissue="Bone Marrow"),
                    color="black", size=6.5, aes(label=c("*", "*", "***","*","*")),hjust = rep(0.5,5))+
        geom_text(data=data.frame(treatment=c(0.9),q3=c(0.335, 0.425), isotype=c(1.1, 1.8), tissue=c("Bone Marrow","Spleen")),
                    color="black", size=6.5, aes(label=c("*","*")),hjust = rep(0.5,2))+
        #geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgG", tissue="Bone Marrow"),
        #            aes(x=0.9, xend=1.1, y=8, yend=8), colour="black",size=0.9) +
         #geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgG", tissue="Bone Marrow"),
         #           aes(x=0.9, xend=1.3, y=8.5, yend=8.5),color="black", size=0.9) +    
         #geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
         #          aes(x=1.9, xend=2.1, y=6.0, yend=6.0),color="black", size=0.9) +   
          #geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
          #         aes(x=1.9, xend=2.3, y=6.5, yend=6.5),color="black", size=0.9) +   
         #geom_text(data=data.frame(treatment=c(0.9,1.0,2.9,1.0),q3=c(8.2, 8.7, 6.2,6.7), isotype=c(1.0,1.1,2.0,2.1)
          #          , tissue="Bone Marrow"),
          #          color="black", size=5, aes(label=c("p<0.001", "p<0.001", "p<0.05","p<0.05")),hjust = rep(0.5,4))+
        theme(legend.position = c(0.12, 0.8),text= element_text(size=12),axis.title.x=element_blank(),
                    strip.text.x = element_text(size = 11, colour = "blue"))
 
 
 
 
    #plot the pie chart of distribution
#    setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/figure3/")
   load( file=here(data.dir,"cloneSize_distribution13.RData"))
    p1<-p1+theme (plot.title = element_text(size=11))
    p2<-p2+theme (plot.title = element_text(size=11))
    p3<-p3+theme (plot.title = element_text(size=11))
    p4<-p4+theme (plot.title = element_text(size=11))
    
#    setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/figure5/subsample")
tiff(file=here(output.dir,"diversity_proportionNorm_figure5_arcsin_v2.0.tiff"), 
    width=1300, height=850)
ggarrange(
                            ggarrange(iso.unnorm, ggarrange(p2, p1,  p4, p3,  ncol = 2, nrow = 2,
                                                                                        common.legend=T, legend="bottom"),
                                                     iso.pbs,  tis.pbs,ncol=4, nrow=1, labels=c("A","B","C","D")) , 
                            ggarrange(
                                                      tr,  q3.tr   , labels=c("E","F"), 
                            ncol=2, nrow=1, widths=c(2,2))
                            
                            , ncol=1, nrow=2  )
dev.off()
      
#tiff(file="diversity_proportionArcsin.tiff", width=900, height=900)
#ggarrange(tr, tis, iso )
#dev.off()
##########################

#
     dt.mod<-data.frame(q3=div.t3[div.t3$q==1,"trans"], treatment=div.t3[div.t3$q==1,"treatment"], 
                                isotype=div.t3[div.t3$q==1,"isotype"], tissue=div.t3[div.t3$q==1,"tissue"])
#mod<-lm(q3~treatment*isotype*tissue, data=dt.mod)
        #----start plotting the individual effects 
#get the mean , treatment
q1.tr<-ggplot(data=dt.mod, aes(y=q3, x=isotype, color=treatment))+
        #geom_point(aes(shape=treatment))+#lims(x=c(0.00, 0.034),y=c(0,0.04))+theme(legend.position = c(0.79, 0.18))+
        geom_boxplot(aes(shape=isotype), notch=F, notchwidth=0.1, varwidth=T)+
        #lims(x=c(0.00, 0.034),y=c(0,0.04))+
        guides(shape="none")+
        labs( y="Normalized Diversity(q=1)")+#, title=name.array[i]) +
       facet_grid(.~tissue)
       
 tiff(file=here(output.dir,"diversity_proportionNorm_suppl_q1_arcsin.tiff"), 
    width=800, height=650)      
       q1.tr +
 #      dev.off()
       
      # + 
        geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=1.7, xend=1.9, y=0.375, yend=0.375), colour="black",size=0.9) +
         geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=1.7, xend=2.1, y=0.35, yend=0.35),color="black", size=0.9) +   
         geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                   aes(x=1.7, xend=2.3, y=0.275, yend=0.275),color="black", size=0.9) +   
          geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                   aes(x=1.9, xend=2.3, y=0.3, yend=0.3),color="black", size=0.9) +
         geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                   aes(x=2.1, xend=2.3, y=0.325, yend=0.325),color="black", size=0.9) +
         #geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgG", tissue="Spleen"),
         #           aes(x=0.7, xend=0.9, y=7.8, yend=7.8), colour="black",size=0.9) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
         #geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgG", tissue="Spleen"),
         #           aes(x=0.7, xend=1.3, y=8.3, yend=8.3),color="black", size=0.9) +
         #geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgG", tissue="Spleen"),
         #           aes(x=1.7, xend=1.9, y=3.8, yend=3.8), colour="black",size=0.9) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
         #geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgG", tissue="Spleen"),
         #           aes(x=1.7, xend=2.1, y=4.3, yend=4.3),color="black", size=0.9)+
         #geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgG", tissue="Spleen"),
         #           aes(x=1.7, xend=2.3, y=4.8, yend=4.8),color="black", size=0.9)+
         #geom_point(data=data.frame(treatment=c(2.2, 2.3,2.5,2.6, 1.3 ),ilr=c(0.47,0.47, 0.83,0.83, 0.2), isotype="IgG", tissue="Bone Marrow"),
         #           color="black", size=2, shape=8) +   
         geom_text(data=data.frame(treatment=c(0.9),q3=c(0.28, 0.3025, 0.3275,0.3525,0.38), isotype=c(2,2.1,2.2,1.9,1.8)
                    , tissue="Bone Marrow"),
                    color="black", size=6.5, aes(label=c("***", "*", "*","*","*")),hjust = rep(0.5,5))+
        #geom_text(data=data.frame(treatment=c(0.9,1.0,2.8,2.9,1.0),q3=c(8.0, 8.5, 4.0, 4.5,5.0), isotype=c(0.8, 1.0,1.8,
        # 1.9, 2.0), tissue="Spleen"),
        #            color="black", size=5, aes(label=c("p<0.05", "p<0.05", "p<0.05","p<0.05","p<0.05")),hjust = rep(0.5,5))+
        #geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgG", tissue="Bone Marrow"),
        #            aes(x=0.9, xend=1.1, y=8, yend=8), colour="black",size=0.9) +
         #geom_segment(data=data.frame(treatment=0.0,q3=0.01, isotype="IgG", tissue="Bone Marrow"),
         #           aes(x=0.9, xend=1.3, y=8.5, yend=8.5),color="black", size=0.9) +    
         #geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
         #          aes(x=1.9, xend=2.1, y=6.0, yend=6.0),color="black", size=0.9) +   
          #geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
          #         aes(x=1.9, xend=2.3, y=6.5, yend=6.5),color="black", size=0.9) +   
         #geom_text(data=data.frame(treatment=c(0.9,1.0,2.9,1.0),q3=c(8.2, 8.7, 6.2,6.7), isotype=c(1.0,1.1,2.0,2.1)
          #          , tissue="Bone Marrow"),
          #          color="black", size=5, aes(label=c("p<0.001", "p<0.001", "p<0.05","p<0.05")),hjust = rep(0.5,4))+
        theme(legend.position = c(0.1, 0.8),text= element_text(size=12),axis.title.x=element_blank())
   dev.off()
