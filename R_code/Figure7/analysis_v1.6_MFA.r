 # R code to read data and then do MFA (PCA)
 #########
 #----- 1/26/2022
 #          add partial point analysis to the grouping variables (supplementary)
 #
 ##    12/9/21---
 #    now decide to run pc 3 (3 pcs of gene usage )
 #          and also save everything to clustering/pc3 
 #          remove mutation sd and cdr3 sd (keep only cdr3 length). 
 ###---------------------
 #  --------------7/18/2021
 #           doing first 2 pcs of gene usage for MFA. 
 #                          check the file analysis_v1.5_MFA.r for first 4 pcs of gene usage 
 ########
 ###---------7/16/2021------
 #   doing analysis for manuscript without the porinB group
 #          copied from file:///home/feng/Windows/windowsD/feng/LAB/MSI/MousePorB_lee_IgSeq/WL03/WL03R2/newPipeline/clustering/MFA/analysis_v1.5_MFA.r
 ####
 #=================
 #  7/10/2021
 # start doing  all variables without BM09
 #          copied from ../analysis_v1.0_MFA.r and ../analysis_v1.5_cluster.r
 #  Also we need to limited the gene usage affects.
 #========================
 #-----Feng 4/13/2021
 
 #read the data (saved in analysis_v1.0_data.R)
 # using FactoMineR
 
 library(FactoMineR)
 library(missMDA)
 library(factoextra)
 library(ggpubr)
 library(ggrepel)
 library(ggalt)
 library(here)
 library(ggplot2)
 
 data.dir<-"Data/Figure7"
 output.dir<-"R_code/Figure7"
 #read data.
#setwd("/home/feng/Windows/windowsD/feng/LAB/MSI/MousePorB_lee_IgSeq/WL03/WL03R2/newPipeline/clustering")
#setwd("E://feng/LAB/MSI/MousePorB_lee_IgSeq/WL03/WL03R2/newPipeline/clustering")
load(here(data.dir,"all_data.RData")) #read in  data containing everthing.  d.all, 
                                                        #saved in "dataPrep.R"
                                                        
#need to do PCA first on gene usage to reduce dimension of usage, otherwise it has too much impact on the output
#
d.all<-d.all[d.all$treatment!="OVA+PorB",]
#conds<-conds[rownames(conds)!="MB9_IgG"&rownames(conds)!="MB9_IgM",]
#d.all.conds<-cbind(d.all, conds[,c("tissue","treatment", "isotype")])

d.all$isotype<-factor(d.all$isotype)
#conds<-d.all[,c(94:99)]
conds<-d.all[,c(86:91)] #<-factors, qualitive factors
#geneUsage<-d.all[,c(1:93)]
geneUsage<-d.all[,c(1:85)]#<-gene usage.
#d.all<-d.all[,-c(94:99, 109:148)]
d.all<-d.all[,-c(86:91, 103:142)] #<----get rid of factor (to move them to the end) and also by clone information 103:142)
d.all.conds<-cbind(d.all, conds[,c("tissue","treatment", "isotype")])

#####################3
#               PCA on gene usage 
####################
pca.VGene.clrd<-prcomp(geneUsage, scale=F, center=T)#previously were center=T and scale=T
contr<-pca.VGene.clrd$sdev^2/sum(pca.VGene.clrd$sdev^2)
gu_pcs<-3#2

#put them back
#d.all<-d.all[,-c(1:93)]
d.all<-d.all[,-c(1:85)]
d.all<-d.all[,-c(2,14)]#get rid of sd of mutation at 2, and cdr3sd at 12.
d.all<-cbind(pca.VGene.clrd$x[,c(1:gu_pcs)], d.all)
d.all.conds<-cbind(d.all, conds[, c("tissue","treatment","isotype")])

#impute first   
#group =c(93  ,   #"usage"
#               2, #"mutation
#               7, #"diversity
#               10, #"CDR length"   
#               10, #"idi"
#               10, #"meanMuFreq
#               10, #"S_index
#               1, #"selection ilr
#)
#
group =c(gu_pcs  ,   #"usage"
               1, #"mutation
               9, #"diversity
#               10, #"CDR length"   
#               10, #"idi"
#               10, #"meanMuFreq
#               10, #"S_index
               1, #"selection ilr
               1, #"CDR3"
               1,#tissue, 
               1,#treatment
               1#isotype
)

n_pcs<-12
d.out <- MFA(d.all.conds, group=group, type=c(rep("s", length(group)-3),"n","n","n"),ncp=n_pcs, 
                     name.group = c("GeneUsage", "Mutation", "Diversity", #"CDR3Length", "Dissimiliarity","meanMuFreq", "S_index", 
                     "Selection","CDR3Length","tissue", "treatment","isotype" ), num.group.sup=c(6,7,8),
                     graph=F)

save(file=here(data.dir,"d.out.pc3.RData"), d.out, d.all.conds)                     

     png(file=here(output.dir,"pcs_cont.png"), width=500, height=400)                     
fviz_screeplot(d.out, addlabels=T)
dev.off()

png(file=here(output.dir,"pcs_group.png"), width=1000, height=900)
g1<-fviz_mfa_var(X=d.out, axes=c(1,2),choice="group", geom=c("point", "text"), repel=T, labelsize=5)+scale_size(range = c(17,18))+
            theme(text = element_text(size=15))
g2<-fviz_mfa_var(X=d.out, axes=c(3,4),choice="group", geom=c("point", "text"), repel=T, labelsize=5)+scale_size(range = c(17,18))+
        theme(text = element_text(size=15))
g3<-fviz_mfa_var(X=d.out, axes=c(5,6),choice="group", geom=c("point", "text"), repel=T, labelsize=5)+scale_size(range = c(17,18))+
            theme(text = element_text(size=15))
g4<-fviz_mfa_var(X=d.out, axes=c(7,8),choice="group", geom=c("point", "text"), repel=T, labelsize=5)+scale_size(range = c(17,18))+
        theme(text = element_text(size=15))   
 ggarrange(g1,g2,g3,g4, ncol=2, nrow=2)       
dev.off()
save(file=here(data.dir,"PCA_scree_plot.RData"), g1,g2,g3,g4 )#,d.out is not saved since it was done above. 
png(file=here(output.dir,"pcs_group_quantiVar.png"), width=750, height=500)
h1<-fviz_mfa_var(d.out, palette = "jco",  axes=c(1,2),
             col.var.sup = "violet", repel = TRUE)+
        theme(text = element_text(size=15), legend.text=element_text(size=18))
h2<-fviz_mfa_var(d.out, palette = "jco",  axes=c(3,4),
             col.var.sup = "violet", repel = TRUE) +
     theme(text = element_text(size=15), legend.text=element_text(size=18))
     
h3<-fviz_mfa_var(d.out, palette = "jco",  axes=c(5,6),
             col.var.sup = "violet", repel = TRUE) +
     theme(text = element_text(size=15), legend.text=element_text(size=18))
     
h4<-fviz_mfa_var(d.out, palette = "jco",  axes=c(7,8),
             col.var.sup = "violet", repel = TRUE) +
     theme(text = element_text(size=15), legend.text=element_text(size=18))
     ggarrange(h1,h2, h3, h4,ncol=2, nrow=2)     
dev.off()         
    
# Contributions to dimension 1
png(here(output.dir,"pcs_group_quantiVar_bar.png"), width=900, height=1200)


g1<-fviz_contrib(d.out, choice = "quanti.var", axes = 1, top = 20,
             palette = "jco")   +theme(text = element_text(size=15), legend.text=element_text(size=18))
g2<-fviz_contrib(d.out, choice = "quanti.var", axes = 2, top = 20,
             palette = "jco")      +theme(text = element_text(size=15), legend.text=element_text(size=18))
g3<-fviz_contrib(d.out, choice = "quanti.var", axes =3, top = 20,
             palette = "jco")      +theme(text = element_text(size=15), legend.text=element_text(size=18))          
g4<-fviz_contrib(d.out, choice = "quanti.var", axes =4, top = 20,
             palette = "jco")      +theme(text = element_text(size=15), legend.text=element_text(size=18)) 
g5<-fviz_contrib(d.out, choice = "quanti.var", axes =5, top = 20,
             palette = "jco")      +theme(text = element_text(size=15), legend.text=element_text(size=18))  
g6<-fviz_contrib(d.out, choice = "quanti.var", axes =6, top = 20,
             palette = "jco")      +theme(text = element_text(size=15), legend.text=element_text(size=18))
g7<-fviz_contrib(d.out, choice = "quanti.var", axes =7, top = 20,
             palette = "jco")      +theme(text = element_text(size=15), legend.text=element_text(size=18))  
g8<-fviz_contrib(d.out, choice = "quanti.var", axes =8, top = 20,
             palette = "jco")      +theme(text = element_text(size=15), legend.text=element_text(size=18))             
  ggarrange(g1, g2, g3,g4,g5,g6,g7,g8, ncol=2, nrow=4, common.legend=T)           
dev.off();       

dd<-as.data.frame(d.out$ind$coord)  
dd<-cbind(dd, conds)           
dd$is_ti<-paste0(dd$isotype,"_", dd$tissue)
dd$is_ti_tr<-paste0(dd$isotype,"_", dd$tissue,"_", dd$treatment)
dd$treatment<-factor(dd$treatment, levels=c("PBS", "OVA","OVA+PorB","OVA+CpG","OVA+Alum"))
dd$tissue<-factor(dd$tissue)
dd$isotype<-factor(dd$isotype)
dd$is_ti<-factor(dd$is_ti)
dd$is_ti_tr<-factor(dd$is_ti_tr)

ind1<-fviz_mfa_ind(d.out, col.ind = "contrib", axes=c(1,2), geom = c("point"),
                habillage=as.factor((conds$isotype)),
            addEllipses =T, ellipse.type = "t", ellipse.level=0.90,#level.conf=0.68,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             col.quali.var.sup="white",pointsize=3,
             repel = TRUE) + theme(text = element_text(size=15),legend.position=c(0.09,0.75))+
                ggtitle("")+
                    guides(color = guide_legend(override.aes = list(size = 3) ) )
 library("plotly")
 #set browser to firefox or other 
 options(browser="firefox")
 plot_ly(x=~Dim.1, y=~Dim.2, z=~Dim.3, type="scatter3d", 
                        data=as.data.frame(d.out$ind$coord),
                    mode="markers", color=(conds$isotype), colors=c("red","blue")#, 
                    #name=rownames(d.quali), symbol=~group,symbols=c('circle','diamond','square')
                    )

ind2<-fviz_mfa_ind(d.out, axes=c(1,3),geom = c("point"),#"text"), 
             habillage = as.factor(paste0(conds$tissue,"+",conds$isotype)), # color by groups 
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "green", "purple"), 
             col.quali.var.sup="white",
             addEllipses = TRUE, ellipse.type = "t", ellipse.level=0.68,#level.conf=0.68,
             repel = TRUE ,# Avoid text overlapping
              shape.ind=19, pointsize=3
             ) +  ggtitle("")+
             theme(text = element_text(size=15), #legend.text=element_text(size=18),
                                        legend.position=c(0.09,0.75)
                                                )+
              guides(color = guide_legend(override.aes = list(size = 3) ) )
 plot_ly(x=~Dim.1, y=~Dim.2, z=~Dim.3, type="scatter3d", 
                        data=as.data.frame(d.out$ind$coord),
                    mode="markers", color=paste0(conds$tissue,"+",conds$isotype), colors=c("#00AFBB", "#E7B800", "#FC4E07", "green", "purple")#, 
                    #name=rownames(d.quali), symbol=~group,symbols=c('circle','diamond','square')
                    )
#ggplot(data=dd, aes(x=Dim.1, y=Dim.2, color=is_ti, group=is_ti))+
#        geom_point(aes(shape=is_ti),size=6)+#lims(x=c(0.00, 0.034),y=c(0,0.04))+theme(legend.position = c(0.79, 0.18))+
#        geom_encircle(expand=0)+#stat_ellipse(level=0.5) +  ###add the circles 
#        theme()

#fviz_mfa_ind(d.out, col.ind = "coord", axes=c(1,2), habillage=as.factor(paste0(conds$tissue,"_", conds$isotype)),
#            addEllipses =T,
#            #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#             repel = TRUE)  

ind2.2<-  ggplot(data=dd, aes(x=Dim.1, y=Dim.3, color=tissue, group=is_ti))+
        geom_point(aes(shape=isotype),size=4)+#lims(x=c(0.00, 0.034),y=c(0,0.04))+theme(legend.position = c(0.79, 0.18))+
        geom_encircle(data=dd,expand=0.0, aes())+#stat_ellipse(level=0.5) +  ###add the circles 
        geom_text_repel(label=rownames(dd), size=5)+
        theme_bw()
plot_ly(x=~Dim.1, y=~Dim.2, z=~Dim.3, type="scatter3d", 
                        data=as.data.frame(d.out$ind$coord),
                    mode="markers", color=paste0(conds$tissue,"+",conds$isotype), colors=c("#00AFBB", "#E7B800", "#FC4E07", "green", "purple"), 
                    #name=rownames(d.quali), 
                    symbol=conds$treatment,symbols=c('circle','diamond','square','o','x')
                    )
                    
ind3<-  ggplot(data=dd, aes(x=Dim.1, y=Dim.2, color=treatment, group=is_ti_tr))+
        geom_point(aes(shape=treatment),size=4)+#lims(x=c(0.00, 0.034),y=c(0,0.04))+theme(legend.position = c(0.79, 0.18))+
        geom_encircle(data=dd,expand=0.0, aes())+#stat_ellipse(level=0.5) +  ###add the circles 
        #geom_polygon(alpha=0.0, aes())+
        #geom_text_repel(label=as.character(dd$ID))+
        theme_bw()+theme(text = element_text(size=15), legend.position=c(0.9,0.75))+
        facet_grid(tissue~isotype)
        
        
 png(here(output.dir,"pcs_group_individual.png"), width=1200, height=1000)            
ggarrange(ggarrange(ind1, ind2, nrow=1, ncol=2, labels=c("A", "B")), 
                        ggarrange(ggplot()+theme_bw()+theme(panel.border = element_blank()),ind3, nrow=1, ncol=3, widths=c(1,4,1),labels=c("", "C")), 
            nrow=2, ncol=1,heights=c(2,3))
 dev.off()            



#fviz_mfa_ind(d.out, geom = c("point","text"), 
#             habillage = as.factor(paste0(conds$isotype, conds$tissue, conds$treatment)), # color by groups 
#             #palette = c("#00AFBB", "#E7B800", "#FC4E07", "green", "purple"),
#             addEllipses = TRUE, ellipse.type = "convex", ellipse.level=0.95,#level.conf=0.68,
#             repel = TRUE #,# Avoid text overlapping
#              #shape.ind=19, pointsize=3,select.ind=list(name=paste0("MB", c(1:15),"_IgM")),
#             ) +theme(text = element_text(size=15), legend.text=element_text(size=18))+
#              guides(color = guide_legend(override.aes = list(size = 3) ) )              


ind1.suppl<-fviz_mfa_ind(d.out, col.ind = "contrib", axes=c(1,2), habillage="none",#as.factor((conds.noBM9$isotype)),
            addEllipses =F,# ellipse.type = "t", ellipse.level=0.95,#level.conf=0.68,
            select.ind=list(name=c("not existing")),#select.var=list(name=c("isotype")),
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)+xlim(c(-0.5, 0.5))
d.quali<-as.data.frame(d.out$quali.var.sup$coord)
             d.quali$group<-"treatment"
             d.quali[c(1:2), "group"]<-"tissue"
             d.quali[c(7:8), "group"]<-"isotype"
             d.quali$group<-factor(d.quali$group)
          d.quali<-as.data.frame(d.out$quali.var.sup$coord)
             d.quali$group<-"treatment"
             d.quali[c(1:2), "group"]<-"tissue"
             d.quali[c(7:8), "group"]<-"isotype"
             d.quali$group<-factor(d.quali$group)
             #draw it with clolor "manually"

#library(umap)
#umap_object<- umap(d.quali[,c(1:7)], scale=F, init= "spectral", min_dist=0.01, n_components=2)             
ind1.suppl.m<-ggplot(d.quali, aes(x=Dim.1, y=Dim.2, colour=group, group=group))+
                    geom_point(aes(shape=group), size=6)+ #xlim(c(-0.5, 0.5))+
                    geom_text_repel(
                                label=rownames(d.quali), 
                                #nudge_x = 0.01, nudge_y = 0.1, 
                                 size=5#,repel=T
                              )
ind1.suppl.m3<-ggplot(d.quali, aes(x=Dim.1, y=Dim.3, colour=group, group=group))+
                    geom_point(aes(shape=group), size=6)+ #xlim(c(-0.5, 0.5))+
                    geom_text_repel(
                                label=rownames(d.quali), 
                                #nudge_x = 0.01, nudge_y = 0.1, 
                                 size=5#,repel=T
                              )
 ind1.suppl.m23<-ggplot(d.quali, aes(x=Dim.2, y=Dim.3, colour=group, group=group))+
                    geom_point(aes(shape=group), size=6)+ #xlim(c(-0.5, 0.5))+
                    geom_text_repel(
                                label=rownames(d.quali), 
                                #nudge_x = 0.01, nudge_y = 0.1, 
                                 size=5#,repel=T
                              )
ggarrange(ind1.suppl.m, ind1.suppl.m3, ind1.suppl.m23,nrow=3, ncol=1)

 png(here(output.dir,"pcs_group_individual_quali_noCloneNoBM9.png"), width=1200, height=500)            
ind1.suppl.m
 dev.off()

### plot the partial points.
partial.isotype<-plot.MFA(d.out, axes=c(1,2),choix="ind", partial=c("IgG","IgM"), invisible='ind')
partial.tissue<-plot.MFA(d.out, axes=c(1,2),choix="ind", partial=c("Spleen","Bone Marrow"), invisible='ind')
partial.immune<-plot.MFA(d.out, choix="ind", partial=c("PBS","OVA","OVA+Alum","OVA+CpG"), invisible='ind') 
ggarrange(partial.isotype, partial.tissue, partial.immune, nrow=2, ncol=2, common.legend=T)

#now drawing the partial points manually, since the build one is not doing good.
d.quali<-as.data.frame(d.out$quali.var.sup$coord)  #this is the barycenter of the partial points
d.quali.partial<-as.data.frame(d.out$quali.var.sup$coord.partiel)
#put them together, get the order first
ord<-rownames(d.quali.partial)
ord<-sub(pattern="\\.[A-Za-z0-9]*","",ord, fix=F )
d.quali.all<-d.quali[ord, ]
rownames(d.quali.all)<-rownames(d.quali.partial)
d.quali.all.two<-cbind(d.quali.all[,c(1,2)], d.quali.partial[,c(1,2)])
#create grouping variables
d.quali.all.two$type<-ord 
ord<-rownames(d.quali.partial)
ord<-sub(pattern="[A-Za-z0-9+ ]*\\.","",ord, fix=F )
d.quali.all.two$variable<-ord
#colnames(d.quali.all.two)<-c("Dim.1","Dim.2","Dim.1","Dim.2", "type","variable")
d.quali.all.two$group<-"Immunization"
d.quali.all.two[d.quali.all.two$type=="Bone Marrow"|d.quali.all.two$type=="Spleen","group"]<-"Tissue"
d.quali.all.two[d.quali.all.two$type=="IgG"|d.quali.all.two$type=="IgM","group"]<-"Isotype"
d.quali.all.two$group<-factor(d.quali.all.two$group, levels=c("Isotype", "Tissue", "Immunization"))
#rearrange so that we can plot them into lines
d.quali.all.two.line<-rbind(d.quali.all.two[,c(1,2,5,6,7)], d.quali.all.two[,c(3,4,5,6,7)])
d.quali.all.two.line$groupVar<-paste0(d.quali.all.two.line$type, d.quali.all.two.line$variable)
d.quali.lab<-d.quali
d.quali.lab$group="Immunization"
d.quali.lab$group[c(1,2)]="Tissue"
d.quali.lab$group[c(7,8)]="Isotype"
d.quali$factors<-rownames(d.quali)
d.quali.all.two.line$variable<-factor(d.quali.all.two.line$variable)
d.quali$factors<-factor(d.quali$factors, levels=c("Bone Marrow","Spleen", "PBS","OVA","OVA+CpG","OVA+Alum", "IgM","IgG"))
tiff(file=here(output.dir,"figure6_MFA_partial_suppl.tiff")
        , width=1200, height=1800)
ggplot(data=d.quali.all.two.line)+
            geom_line(aes(x=Dim.1, y=Dim.2, colour=variable, group=groupVar, linetype=variable), size=1.5)+ 
            guides( shape=guide_legend(title="Group"))+
            geom_point(data=d.quali,aes(x=Dim.1, y=Dim.2, shape=factors),  size=5.5)+
            scale_shape_manual(values=c(0,1,18,15,16,17,9,10))+
            scale_linetype_manual(values=c("solid","dashed","twodash", "F1","longdash"))+
            #scale_color_manual(values=c("red","red","red","red","red"))+
            geom_point(aes(x=Dim.1, y=Dim.2, colour=variable, group=groupVar), size=2.5, shape=0)+
            #geom_text(data=d.quali.lab,
            #        color="black", size=8, aes(label=rownames(d.quali.lab),x=Dim.1, y=Dim.2))+
            facet_wrap(.~group, ncol=2)+
            theme_bw(base_size=25)+theme(legend.position=c(0.75,0.3))
dev.off()
save(file=here(data.dir,"fig6_PCA_partial.RData"),
    d.quali.all.two.line, d.quali )
#now plotting.
#plotting 

temp<-d.quali.all.two.line[d.quali.all.two.line$group=="Immunization",]
ggplot(data=temp)+
            geom_line(aes(x=Dim.1, y=Dim.2, colour=variable, group=groupVar, linetype=variable), size=1.5)+ 
            guides( shape=guide_legend(title="Group"))+
            geom_point(data=d.quali,aes(x=Dim.1, y=Dim.2, shape=factors),  size=5.5)+
            scale_shape_manual(values=c(0,1,18,15,16,17,9,10))+
            scale_linetype_manual(values=c("solid","dashed","twodash", "F1","longdash"))+
            #scale_color_manual(values=c("red","red","red","red","red"))+
            #geom_point(aes(x=Dim.1, y=Dim.2, colour=variable, group=groupVar), size=2.5, shape=0)+
            #geom_text(data=d.quali.lab,
            #        color="black", size=8, aes(label=rownames(d.quali.lab),x=Dim.1, y=Dim.2))+
            #facet_wrap(.~group, ncol=2)+
            theme_bw(base_size=25)+theme(legend.position="none")
        #dev.off()
    library("plotly")
    library(plot3D)
    library(scatterplot3d)
    anx<-50
      png(here(output.dir,"pcs_group_individual_quali_3d.png"),
       width=1200,height=500)
  #scatterplot is easier to manipulate the behavior , do3d and rotate3d calling/wrapper of scatterplot3d
  scatterplot3d(angle=anx, xlab="PC1", ylab="PC2",zlab="PC3",
        x=d.quali[,1],y=d.quali[,2], z=d.quali[,3],
                    color=as.integer(as.factor(d.quali$factors))+1, 
                    pch=c(19,19, 8,9,10,12,15,15),#pch=as.integer(as.factor(rownames(d.quali))),
                    cex.symbols=2
                    )
  legend(x=-2.8,y=-0.1, legend=rownames(d.quali),pch=c(19,19, 8,9,10,12,15,15),
            col=as.integer(as.factor(d.quali$factor))+1
  )
  dev.off()
  
  save(file=here(data.dir,"data_for_3d_gifmovie.RData"),
    dd, d.quali)
    
  ##########################################
  ###doing rotation movie gif 
  ##########################################
  # need to go to the temp directory to manually copy over the movie gif.
  

#######################NOTE: be really careful, the following chunk are optimized for pc2 data. don't change it yet.
#############doing hiarchecal 
 library(cluster)
 library(dendextend)
 x.dist<-daisy(d.out$ind$coord, metric="e", stand=F)
hh<-hclust(x.dist, method="ward.D")#"single")
dx<-reorder(as.dendrogram(hh),1:dim(as.matrix(x.dist))[1])
plot(dx)           


d.all.conds$treatment<-factor(d.all.conds$treatment, levels=c("PBS","OVA", "OVA+CpG", "OVA+Alum"))
 #dx<-reorder(as.dendrogram(hh),1:dim(as.matrix(x.dist))[1])
 #####note: about reordering...........#############
 ## it was in order of sp IgG 1:11, sp IgM 12:22, BM IgG 23:34, BM IgM 35:46
 #dx<-reorder(as.dendrogram(hh),c(c(1:5,1,7:8,50,10:11),c(c(12:13)+550,14:22)+200,(23:34)+20,c(38,37,36,35, 39:46)+500), agglo.FUN=mean)
 #          the above is the original order when doing the input using /PorB_2nd_try data
 #dx<-reorder(as.dendrogram(hh),c(c(1:6,7:8,9,10:11)+200,c(c(12:13),14:17,18:22)+800,
#                                                                c(c(23:28)+250,c(29:31), 32+1.5,33:34)+30,c(38,37,36,35, 39:46)+500), agglo.FUN=mean)
dx<-reorder(as.dendrogram(hh),c(c(1:3,4+400,5:6,7:8,9,10:11)+500,c(c(12:13)+80,14:17,18:21,22)+100,
                                                                c(c(23:28)+250,c(29:31), 32+1.5,33, 34+20)+800,c(35:37,38-5, 39:40, 41+100,c(42:43)+50,c(44:45)+100, 46+50)), agglo.FUN=mean)
                                                                
 #dx<-reorder(as.dendrogram(hh),c(c(1:5,1,7:8,50,10:11)+200,c(c(12:13)+550,14:17,18:22),c(23:28,c(29:31)+10, 32:34)+30,c(38,37,36,35, 39:46)+800), agglo.FUN=mean)
 order.x<-order.dendrogram(dx)
 
 #order.x<-order.x[c(c(1:12), c(16,15),c(13,14),c(17:44),c(46,45))]
 #dx<-reorder(as.dendrogram(hh), order.x)
 #order.x<-order.dendrogram(dx)
 bars<-d.all.conds[order.x,c("isotype","tissue","treatment")]

 #bars<-as.data.frame()
 
 cols_2<-colorspace::rainbow_hcl(2, c=70,l=50)
 cols_5<-colorspace::heat_hcl(5, c=70,l=50)
 cols_4<-colorspace::rainbow_hcl(4, c=70,l=50)
 cols_4<-c("blue","black","pink","red")
 dx<-color_branches(dx, col=cols_4[as.factor(paste0(bars$isotype,bars$tissue))])
 dx<-color_labels(dx, col=cols_4[as.factor(paste0(bars$isotype,bars$tissue))])
 
 bars$isotype<-cols_2[bars$isotype]
bars$tissue<-cols_2[bars$tissue]
 bars$treatment<-cols_4[bars$treatment]
 colnames( bars)[3]<-"immunization"
 label.sample<-rownames(d.all.conds)
 label.sample<-substr(label.sample, 0, nchar(label.sample)-4)
 label.sample<-label.sample[order.x]
 #setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/clustering/pc3/")
 png(file=here(output.dir,"hiearchicalCluster.png"), 
    width=1000, height=400)
 dx %>% set("labels", label.sample) %>%set("branches_lwd", 2.5) %>%plot
 #plot(dx)
 colored_bars(colors = bars[,c(2:3)], dend = dx, sort_by_labels_order=F, y_shift=-4.5)
 lines(x=c(23.5,23.5),y=c(-10,30), col="red", lty=2, lwd=4)
text(x=35, y=18, label=c("IgG"), col="blue",cex=2.0 )
text(x=12, y=18, label=c("IgM"), col="red",cex=2.0)
legend(x=40, y=25, legend = levels(d.all.conds$treatment), fill = cols_4, title="Immunization")
legend(x=40,y=15, legend = levels(d.all.conds$tissue), fill = cols_2, title="Tissue")
 
dev.off()
#save(dx, bars,d.all.conds, cols_4, dd, d.out, file=("MFA.hclust.plots.v1.6.rData"))
save(#dx, bars,
                            d.all.conds, #cols_4, 
                            dd, d.out, 
                            file=here(data.dir,"MFA.hclust.plots.v1.6_pc3.rData"))


