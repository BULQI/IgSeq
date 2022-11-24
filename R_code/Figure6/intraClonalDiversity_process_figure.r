#R code to process intraclonal diversity results   ------- for figure 6 of the manuscript

#library
library(plyr)
library(ggplot2)
library(frLib)
library(rcPkg)
library(vegan)
library(compositions)
library(grid)
library(gridBase)

 
 load("./clones.df.RData")
 topn<-40#40#30#15 #20
 
 #
 load("intraClonalDiversity_BM.RData")  #idis.BM.IgM, idis.BM.IgG
#generated in ./intraClonal_batch.R, 
 load("intraClonalDiversity_SP.RData")  #idis.SP.IgM, idis.SP.IgG
#generated in ./intraClonal_batch.R  


                    #data.array<-list(mi.BM.IgM.top20, mi.BM.IgG.top20, mi.SP.IgM.top20, mi.SP.IgG.top20)
                    #name.array<-c("Bone Marrow (IgM)", "Bone Marrow (IgG)", "Spleen (IgM)", "Spleen (IgG)");
                    #save(data.array, conditions, name.array, file="intraClonalDiversityDataArray.RData")
#load("intraClonalDiversityDataArray_top75.RData")
#load("intraClonalDiversityDataArray_top50.RData")#       
load("intraClonalDiversityDataArray_top40.RData")            


#first do linear regression, 
library(MASS)

#remove "porB" samples for doing the manuscript
for(i in 1:4)
{
    data.array[[i]]<-data.array[[i]][!is.element(data.array[[i]]$sampleName, c("SP7", "SP8", "SP9", "MB7", "MB8", "MB9")),]
}
conditions<-conditions[conditions$treatment !="OVA+PorB",]
conditions$treatment<-factor(conditions$treatment, levels=c("PBS", "OVA", "OVA+CpG", "OVA+Alum"))
mi.temp<-data.array[[1]]   #<-------------build threshold based on  BM IgM clones.
x<-rlm(MeanMuFreq~idi, data=mi.temp,psi=psi.huber)

#####Note: now the data don't include PorB.
#
png("selection_gate.png", width=610, height=600)
plot(mi.temp$idi, mi.temp$MeanMuFreq, type="p", 
            xlab="Intra-clonal dissimilarity", ylab="Mean mutation frequency", main="Clonal Selection",
            xlim=c(0, 0.038), ylim=c(0,0.038))
cx<-seq(0, 0.05, by=0.001)
cy<-cx*x$coefficients[2]+x$coefficients[1]
lines(cx,cy, col=2, lwd=2, lty=2)
cy2<-cx*x$coefficients[2]+0.01
lines(cx,cy2, col=3, lwd=2, lty=3)
dev.off();

#calculating the correlation
#get the unselected data to do correction estimation
categorizeCloneIntraDiversity_unselected.data<-function(dat, slope, intercept)
{
    c_x<-dat[,1]
    c_y<-dat[,2]
    return(dat[slope*c_x+intercept>=c_y,])
}

x.data<-categorizeCloneIntraDiversity_unselected.data(mi.temp[,c("idi", "MeanMuFreq")], x$coefficients[2],0.01 )
cor(x.data[,1],x.data[,2])
#### [1] 0.9647525
#fit a data using negative clone data.
x.negative<-rlm(MeanMuFreq~idi, data=x.data,psi=psi.huber )
        ####
        #> x.negative
         #           Call:
         #           rlm(formula = MeanMuFreq ~ idi, data = x.data, psi = psi.huber)
         #           Converged in 9 iterations
#
         #            Coefficients:
         #           (Intercept)         idi 
         #           0.003220012 0.582808543 
        #
         #           Degrees of freedom: 234 total; 232 residual
         #           Scale estimate: 0.000812 
         #           > x
        #            Call:
        #            rlm(formula = MeanMuFreq ~ idi, data = mi.temp, psi = psi.huber)
         #           Converged in 10 iterations
        #    
        #            Coefficients:
         #           (Intercept)         idi 
         #           0.003252463 0.584879257 


# with regular 
x.negative.reg<-lm(MeanMuFreq~idi, data=x.data)
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
    legend.label<-as.character(conditions$treatment)[order(as.integer(conditions$treatment))]
    legend(x=0.025, y=0.015, legend=unique(legend.label), col=unique(as.integer(conditions$treatment))
                    , pch=unique(as.integer(conditions$treatment)), cex=0.8)
}
par(op)

library(ggalt)
library(ggpubr)
#now let's do ggplot2 with polygon range around the points
#op<-par(mfrow=c(2,2), mar=c(3,2.5,2,0.2), mgp=c(1.5,0.4,0.0))
#tiff("intraClonalDiversity_threshold.tiff", width=800, height=700)
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
library(ggpubr)
figure<-ggarrange(gp[[1]], gp[[2]], gp[[3]], gp[[4]], ncol=2, nrow=2)
#par(op)
figure
#dev.off()

load("intraClonalDiversity_intraD_top40.RData")
#load("intraClonalDiversity_intraD_top75.RData")
#interaction


tiff("interaction.tiff", width=1100, height=700)
op<-par(mfrow=c(2,2), mar=c(3,2.5,2,1.2), mgp=c(1.5,0.4,0.0), cex=1.25)
interaction.plot(intraD$isotype, intraD$tissue,intraD.ilr[,1], xlab="isotype", 
            ylab="% diversified clone (transformed)" )
interaction.plot( intraD$treatment,intraD$isotype,intraD.ilr[,1], xlab="treatment", 
            ylab="% diversified clone (transformed)" )
interaction.plot(intraD$treatment,intraD$tissue, intraD.ilr[,1], xlab="treatment", 
            ylab="% diversified clone (transformed)" , cex=1.5)
par(op)
dev.off()

intraD<-intraD[intraD$treatment!="OVA+PorB",]
#test with the full model.
options(contrasts = c("contr.sum", "contr.poly"))
lmD<-lm(ilr~tissue*isotype*treatment, data=intraD)
library(car)
Anova(lmD, type=3) #type II SSE
#lmD<-lm(intraD.clr[,1]~intraD$tissue:intraD$isotype*intraD$treatment)

#IgG BM, treatment
#intraD.ilr.bg<-intraD.ilr[intraD$tissue=="Spleen"&intraD$isotype=="IgG"]
#intraD.ilr.cond<-intraD[intraD$tissue=="Spleen"&intraD$isotype=="IgG",]
#lmD.bg<-lm(intraD.ilr.bg~intraD.ilr.cond$treatment)

#intraD.ilr.bm<-intraD.ilr[intraD$tissue=="Spleen"]
#intraD.ilr.cond<-intraD[intraD$tissue=="Spleen",]
#lmD.bm<-lm(intraD.ilr.bm~intraD.ilr.cond$treatment*intraD.ilr.cond$isotype)

#emmeans
library(emmeans)

dt.mod<-intraD[,c("ilr","treatment", "tissue", "isotype")]#data.frame(ilr=intraD$ilr, treatment=intraD$treatment, isotype=intraD$isotype, tissue=intraD$tissue)
mod<-lm(ilr~treatment*isotype*tissue, data=dt.mod)
em_tr<-emmeans(mod, ~treatment|tissue*isotype)
#em_tr %>% test(joint=F)

Set1 <- list(
  #"OVA/Alum - PBS" = c(-1,0,0,0,1), "OVA/CpG - PBS"=c(-1,0,0,1,0), "OVA/PorB - PBS"=c(-1,0,1,0,0), "OVA - PBS"=c(-1, 1,0,0,0)
    "OVA/Alum - PBS" = c(-1,0,0,1), "OVA/CpG - PBS"=c(-1,0,1,0), "OVA - PBS"=c(-1, 1,0,0),
        "OVA/Alum - OVA" = c(0,-1,0,1), "OVA/CpG - OVA"=c(0,-1,1,0)
      #  "OVA/CpG - OVA "=c(0, -1,1,0), "OVA/Alum- OVA "=c(0, -1,0,1)
  )
emmeans::contrast (em_tr, Set1, adjust ='F')
ec<-emmeans::contrast(em_tr, Set1, adjust='F')
png("intra_treatmentEffect.png", width=850, height=1000)
plot(ec)+geom_segment(x=0, xend=0,y=0, yend=3.8, linetype=2, colour="red", size=1.5)+
    theme(text = element_text(size=25))
dev.off()

#for tissue
em_tr<-emmeans(mod, ~tissue|treatment*isotype)
#em_tr %>% test(joint=F)

Set1 <- list(
  "Bone Marrow - Spleen" = c(-1,1))#, "OVA/CpG - PBS"=c(-1,0,0,1,0), "OVA/PorB - PBS"=c(-1,0,1,0,0), "OVA - PBS"=c(-1, 1,0,0,0))
emmeans::contrast (em_tr, Set1, adjust ='none')
ec<-emmeans::contrast(em_tr, Set1, adjust='none')
png("intra_tissueEffect.png", width=850, height=1000)
plot(ec)+geom_segment(x=0, xend=0,y=0, yend=3.8, linetype=2, colour="red", size=1.5)+
    theme(text = element_text(size=20))
dev.off()

#for isotype
em_tr<-emmeans(mod, ~isotype|treatment*tissue)
#em_tr %>% test(joint=F)

Set1 <- list(
  "IgG - IgM" = c(-1,1))#, "OVA/CpG - PBS"=c(-1,0,0,1,0), "OVA/PorB - PBS"=c(-1,0,1,0,0), "OVA - PBS"=c(-1, 1,0,0,0))
emmeans::contrast (em_tr, Set1, adjust ='none')
ec<-emmeans::contrast(em_tr, Set1, adjust='none')
png("intra_isotypeEffect.png", width=850, height=1000)
plot(ec)+geom_segment(x=0, xend=0,y=0, yend=3.8, linetype=2, colour="red", size=1.5)+
    theme(text = element_text(size=20))
dev.off()




#now plot the difference
tiff("treatmentEffect_top40.tiff", width=900, height=800)
te<-ggplot(data=dt.mod, aes(y=ilr, x=isotype, color=treatment))+
        #geom_point(aes(shape=treatment))+#lims(x=c(0.00, 0.034),y=c(0,0.04))+theme(legend.position = c(0.79, 0.18))+
        geom_boxplot(aes(shape=treatment), notch=F, notchwidth=0.1)+
        #lims(x=c(0.00, 0.034),y=c(0,0.04))+
        guides(shape="none",color=guide_legend(title=element_blank()))+
        labs(x=element_blank(), y="Selection (logRatio %)")+#, title=name.array[i]) +
       facet_grid(.~tissue)+ 
       #geom_hline(data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"), aes(yintercept=yintercept),
       #                     linetype="dashed", color="red")
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=1.7, xend=2.1, y=0.4, yend=0.4), colour="black",size=0.9) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
         geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=1.7, xend=2.3, y=0.6, yend=0.6),color="black", size=0.9) +   
            geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=1.9, xend=2.1, y=0.15, yend=0.15),color="black", size=0.9) +   
         #geom_point(data=data.frame(treatment=c(2.2, 2.3,2.5,2.6, 1.3 ),ilr=c(0.47,0.47, 0.83,0.83, 0.2), isotype="IgG", tissue="Bone Marrow"),
         #           color="black", size=2, shape=8) +   
         geom_text(data=data.frame(treatment=c(2.65,3.1,1.7),ilr=c(0.48, 0.68,0.25), isotype=c(1.9,2.0, 2), tissue="Bone Marrow"),
                    color="black", size=4, aes(label=c("p=0.0098", "p=0.05", "p=0.05")))+
        theme(legend.position = c(0.15, 0.88),text= element_text(size=18))
te
dev.off()


            
##################################
#  ----testing the tissue difference-------#
##################################
#mod<-lm(ilr~treatment*isotype*tissue, data=dt.mod)
em_tr<-emmeans(mod, ~tissue|treatment*isotype)
Set1 <- list(
  "Bone Marrow - Spleen" = c(-1,1))#, "OVA/CpG - PBS"=c(-1,0,0,1,0), "OVA/PorB - PBS"=c(-1,0,1,0,0), "OVA - PBS"=c(-1, 1,0,0,0))
ec<-emmeans::contrast (em_tr, Set1, adjust ='none')
tiff("treatmentTissue_top40.tiff", width=900, height=800)
ggplot(data=dt.mod, aes(y=ilr, x=treatment, color=tissue))+
        #geom_point(aes(shape=treatment))+#lims(x=c(0.00, 0.034),y=c(0,0.04))+theme(legend.position = c(0.79, 0.18))+
        geom_boxplot(aes(shape=tissue), notch=F, notchwidth=0.1)+
        #lims(x=c(0.00, 0.034),y=c(0,0.04))+
        guides(shape=FALSE)+
        labs(x="Treatment", y="Selection (logRatio %)")+#, title=name.array[i]) +
       facet_grid(.~isotype)+ 
        geom_segment(data=data.frame(treatment="OVA+Alum",ilr=0.01, isotype="IgG", tissue=2),
                    aes(x=3.75, xend=4.25, y=0.4, yend=0.4), colour="black",size=0.8) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=2.75, xend=3.25, y=0.4, yend=0.4),color="black", size=0.9) +  
        #geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgM", tissue="Bone Marrow"),
        #            aes(x=1.75, xend=2.25, y=-0.5, yend=-0.5), color="black", size=0.9) +                      
        
         geom_text(data=data.frame(treatment=c("OVA+CpG","OVA+Alum"),ilr=c(0.5,0.5), 
                        isotype=c("IgG","IgG"), tissue=rep(1.5,2)),
                    color="black", size=4, aes(label=c("p=0.02", "p=0.001")))+
        theme(legend.position = c(0.89, 0.15), text= element_text(size=18))
dev.off()
####isotype
em_tr<-emmeans(mod, ~isotype|treatment*tissue)
Set1 <- list(
  "IgG - IgM" = c(-1,1))#, "OVA/CpG - PBS"=c(-1,0,0,1,0), "OVA/PorB - PBS"=c(-1,0,1,0,0), "OVA - PBS"=c(-1, 1,0,0,0))
emmeans::contrast (em_tr, Set1, adjust ='none')

tiff("treatmentIsotype_top40.tiff", width=900, height=800)
ggplot(data=dt.mod, aes(y=ilr, x=treatment, color=isotype))+
        #geom_point(aes(shape=treatment))+#lims(x=c(0.00, 0.034),y=c(0,0.04))+theme(legend.position = c(0.79, 0.18))+
        geom_boxplot(aes(shape=isotype), notch=F, notchwidth=0.1)+
        #lims(x=c(0.00, 0.034),y=c(0,0.04))+
        guides(shape=FALSE)+
        labs(x="Treatment", y="Selection (logRatio %)")+#, title=name.array[i]) +
       facet_grid(.~tissue)+ 
        geom_segment(data=data.frame(treatment="OVA+Alum",ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=3.75, xend=4.25, y=0.4, yend=0.4), colour="black",size=0.8) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
        geom_text(data=data.frame(treatment="OVA+Alum",ilr=0.45, isotype="IgG", tissue="Bone Marrow"),
                    color="black", size=4, aes(label=c("p<0.0001")))+
        geom_segment(data=data.frame(treatment="OVA+CpG",ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=2.75, xend=3.25, y=0.4, yend=0.4), colour="black",size=0.8) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
        geom_text(data=data.frame(treatment="OVA+CpG",ilr=0.45, isotype="IgG", tissue="Bone Marrow"),
                    color="black", size=4, aes(label=c("p<0.0001")))+
        geom_segment(data=data.frame(treatment="OVA+CpG",ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=1.75, xend=2.25, y=0.15, yend=0.15), colour="black",size=0.8) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
        geom_text(data=data.frame(treatment="OVA",ilr=0.23, isotype="IgG", tissue="Bone Marrow"),
                    color="black", size=4, aes(label=c("p<0.0001")))+
        geom_segment(data=data.frame(treatment="PBS",ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=0.75, xend=1.25, y=-0.5, yend=-0.5), colour="black",size=0.8) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
        geom_text(data=data.frame(treatment="PBS",ilr=-0.2, isotype="IgG", tissue="Bone Marrow"),
                    color="black", size=4, aes(label=c("p=0.01")))+
        #geom_segment(data=data.frame(treatment="OVA",ilr=0.03, isotype="IgG", tissue="Spleen"),
        #            aes(x=1.75, xend=2.25, y=-0.35, yend=-0.35), colour="black",size=0.8) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
        #geom_text(data=data.frame(treatment="OVA",ilr=-0.28, isotype="IgG", tissue="Spleen"),
        #            color="black", size=4, aes(label=c("p=0.02")))+
       # geom_segment(data=data.frame(treatment="PBS",ilr=0.01, isotype="IgG", tissue="Spleen"),
       #             aes(x=0.75, xend=1.25, y=-0.55, yend=-0.55), colour="black",size=0.8) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
        #geom_text(data=data.frame(treatment="PBS",ilr=-0.45, isotype="IgG", tissue="Spleen"),
         #           color="black", size=4, aes(label=c("p=0.02")))+
        geom_segment(data=data.frame(treatment="OVA+Alum",ilr=0.01, isotype="IgG", tissue="Spleen"),
                    aes(x=2.75, xend=3.25, y=-0.35, yend=-0.35), colour="black",size=0.8) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
        geom_text(data=data.frame(treatment="OVA+CpG",ilr=-0.28, isotype="IgG", tissue="Spleen"),
                    color="black", size=4, aes(label=c("p=0.01")))+
        #geom_point(data=data.frame(treatment=c(2.5, 2.55,3,3.05),ilr=c(0.47,0.47, 0.82,0.82), isotype="IgG", tissue="Bone Marrow"),
        #            color="black", size=2, shape=8) +   
        theme(legend.position = c(0.97, 0.10), text= element_text(size=18))
dev.off()

####start doing the ploting of the figure 6
png("figure6_top40.png", width=800, height=1100)
grid.newpage()
pushViewport(viewport(layout = grid.layout(7, 6)))
#pushViewport(viewport(layout.pos.row = ceiling(i/2), layout.pos.col = c(1:2)+(1-i%%2)*2))
#plotting the base plot 
for (i in 1:4)
{
    pushViewport(viewport(layout.pos.row = ((ceiling(i/2)-1)*2+1):((ceiling(i/2)-1)*2+2), layout.pos.col = c(1:3)+(1-i%%2)*3))
    if(i==1){
            par(omi=gridOMI(), mar=c(2,2,0.4,0.1), mgp=c(1,0.2,0), tcl=-0.25, cex=1.3)
    }
    else
    {
        par(omi=gridOMI(),mar=c(2,2,0.4,0.1), mgp=c(1,0.2,0), tcl=-0.25, new=TRUE, cex=1.3)
    }

    plot(x=c(0.00, 0.034), y=c(0,0.04), type="n", xlab="Intra-Clonal Dissimilarity", 
                                                    ylab="Mean Mutation Frequency",main="")
             title(name.array[i], line=-1)                                       
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
    legend.label<-as.character(conditions$treatment)[order(as.integer(conditions$treatment))]
    legend(x=0.025, y=0.015, legend=unique(legend.label), col=unique(as.integer(conditions$treatment))
                    , pch=unique(as.integer(conditions$treatment)), cex=0.8)
    par(op)
    popViewport()
}
     
#plot first treatment effects 
te2<-ggplot(data=dt.mod[dt.mod$isotype=="IgG"&dt.mod$tissue=="Bone Marrow",], aes(y=ilr, x=treatment, color=treatment, group=treatment))+
        #geom_point(aes(shape=treatment))+#lims(x=c(0.00, 0.034),y=c(0,0.04))+theme(legend.position = c(0.79, 0.18))+
        geom_boxplot(aes(shape=treatment), notch=F, notchwidth=0.1)+
        #lims(x=c(0.00, 0.034),y=c(0,0.04))+
        guides(shape=FALSE)+
        labs(x="", y="Selection (logRation %)")+#, title=name.array[i]) +
       #facet_grid(isotype~tissue)+ 
       #geom_hline(data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"), aes(yintercept=yintercept),
       #                     linetype="dashed", color="red")
        geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=1, xend=3, y=0.4, yend=0.4), colour="black",size=0.9) +#data=data.frame(yintercept=-0.02, isotype="IgG", tissue="Bone Marrow"),
         geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=1, xend=4, y=0.75, yend=0.75),color="black", size=0.9) +   
            geom_segment(data=data.frame(treatment=0.0,ilr=0.01, isotype="IgG", tissue="Bone Marrow"),
                    aes(x=1, xend=2, y=0.08, yend=0.08),color="black", size=0.9) +   
         geom_point(data=data.frame(treatment=c(2.0, 2.1,2.5,2.6, 1.3 ),ilr=c(0.47,0.47, 0.83,0.83, 0.2), isotype="IgG", tissue="Bone Marrow"),
                    color="black", size=2, shape=8) +   
         geom_text(data=data.frame(treatment=c(2.65,3.1,1.7),ilr=c(0.48, 0.83,0.2), isotype="IgG", tissue="Bone Marrow"),
                    color="black", size=4, aes(label=c("p=0.011", "p=0.009", "p=0.03")))+
         annotate("text", x = 3, y = -1.5, label = "Bone Marrow IgG", size=5)+
        theme(legend.position ="none",text= element_text(size=15),axis.text.x = element_text(angle = 30), axis.title.x = element_blank())
vp <- pushViewport(viewport(layout.pos.row =5:7, layout.pos.col =1:6))

####we are plot the one in above, not this current one.
print(te, vp=vp, newpage=FALSE) #<-plot the non-productive reads
popViewport()
dev.off()



