#R code to load the previously saved ggplot2 objects for plotting 
#    figures.
# Make sure you have read ./ReadMe.txt for instructions
# of running scripts and preparing data.


#library

 library(FactoMineR)
 library(missMDA)
 library(factoextra)
 library(ggpubr)
 library(ggrepel)
 library(ggalt)
 
  library(cluster)
 library(dendextend)
 library(here)
 
 
data.dir<-"Data/Figure7"
output.dir<-"R_code/Figure7"
# load data, clustering  
load(here(data.dir,"g3.grey.rdata")) # saved in analysis_v1.0_cluster.r
                                                  #load g3.grey and conds



#load data, MFA and hclust
load(here(data.dir,"MFA.hclust.plots.v1.6_pc3.rData"))# saved in analysis_v1.5_MFA.r
                                #load dx and d.all.conds and dd, bars, cols_4

#now plotting them together.

ind1<-fviz_mfa_ind(d.out, col.ind = "contrib", axes=c(1,2), geom = c("point"),
                habillage=as.factor((conds$isotype)),
            addEllipses =T, ellipse.type = "t", ellipse.level=0.90,#level.conf=0.68,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             col.quali.var.sup="white",pointsize=3,
             repel = TRUE) + theme(text = element_text(size=10),legend.position=c(0.79,0.3))+
                ggtitle("")+
                    guides(color = guide_legend(override.aes = list(size = 3) ) )+
                    #labs(tag="A")+
        theme(legend.title = element_blank()#plot.tag=element_text(margin=margin(t=40,b=-20),size=18)
        )

ind2<-fviz_mfa_ind(d.out, axes=c(1,2),geom = c("point"),#"text"), 
             habillage = as.factor(paste0(conds$tissue,"+",conds$isotype)), # color by groups 
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "green", "purple"), 
             col.quali.var.sup="white",
             addEllipses = TRUE, ellipse.type = "t", ellipse.level=0.68,#level.conf=0.68,
             repel = TRUE ,# Avoid text overlapping
              shape.ind=19, pointsize=3
             ) +  ggtitle("")+
             theme(text = element_text(size=10), legend.title = element_blank(),#legend.text=element_text(size=18),
                                        legend.position=c(0.79,0.3)
                                                )+
              guides(color = guide_legend(override.aes = list(size = 3) ) )+#labs(tag="B")+
        theme(#plot.tag=element_text(margin=margin(t=40,b=-20),size=18)
        )

ind3<-  ggplot(data=dd, aes(x=Dim.1, y=Dim.2, color=treatment, group=is_ti_tr))+
        geom_point(aes(shape=treatment),size=4)+#lims(x=c(0.00, 0.034),y=c(0,0.04))+theme(legend.position = c(0.79, 0.18))+
        geom_encircle(data=dd,expand=0.0, aes())+#stat_ellipse(level=0.5) +  ###add the circles 
        #geom_polygon(alpha=0.0, aes())+
        #geom_text_repel(label=as.character(dd$ID))+
        theme_bw()+theme(text = element_text(size=15), legend.position=c(0.85,0.25), legend.title=element_blank())+
        facet_grid(tissue~isotype)+#labs(tag="C")+
        theme(#plot.tag=element_text(margin=margin(t=40,b=-20),size=18)
        )
       

g1<-ggplot(aes(x = X, y = Y), data = tsne_df) +#labs(tag="F")+
        geom_point(aes(shape =conds$tissue, color=conds$isotype),size=6 )+#,size=(x.complications.sum+1)*3)+
        #ggtitle("Isotype+Tissue")+ 
        xlab("umap X")+ylab("umap Y")+
        guides(col=guide_legend(title="Isotype (color)",override.aes = list(size=6)),
         shape=guide_legend(title="Tissue (shape)",override.aes = list(size=6)))+
         theme(legend.title=element_text(size=10), 
    legend.text=element_text(size=10), legend.position="bottom",
      plot.title=element_text(margin=margin(t=40,b=-30), hjust=1)#,
      #plot.tag=element_text(margin=margin(t=40,b=-30),size=18)
      )#+
    #geom_text_repel(label=paste(rownames(x.clean),as.factor(tsne_df$cluster)))
 

  
    g3.grey<-ggplot(aes(x = X, y = Y), data = tsne_df) +labs(tag="G")+
            geom_point(aes(shape=factor(conds$treatment)), colour=cols, size=6)+#, size=(x.complications.sum+1)*2)+
                ggtitle("Clustering")+guides(col=guide_legend(title="Cluster",override.aes = list(size=6)),
                 shape=guide_legend(title="Immunization",override.aes = list(size=6)))+
                 scale_shape_manual(values=c(19,17,15,7))+
                 theme(legend.title=element_text(size=10), 
                                                    legend.text=element_text(size=10), legend.position="bottom",
      plot.title=element_text(margin=margin(t=40,b=-30), hjust=1),
      plot.tag=element_text(margin=margin(t=40,b=-30),size=18),
      plot.margin=margin(t = 0, r = 10, b = 0, l = 0, unit = "pt")
      )
               # geom_text_repel(label=paste(rownames(x.clean),as.factor(v_small)))  
        
        ###drawing cluster by tissue and isotype compartment
        load(file=here(data.dir,
            "figure_cluster_data_by_compartment.RData")) #load g3.ti, g4.ti and data set for plotting clustering                                                                   
                                                                        #the 
                                                                        #saved in analysis_v1.6_cluster.r
                                                                        
         #tissue<-c("Spleen","Bone Marrow")
         #isotype<c("IgM","IgG")
         index<-0
         for(i in tissue)
         {
            for (j in isotype){
                index<-index+1
                temp<-xdat.umap.ti[[index]]
                g4.ti[[index]]<-ggplot(aes(x =X, y = Y), data =temp) +
                        geom_point(aes(shape=treatment,color=cluster),size=5)+#, size=(x.complications.sum+1)*2)+
                        #ggtitle(paste0(i,"+",j))+
                        scale_shape_manual(values=c(19,17,15,8,7))+guides(color="none")+
                        theme_bw(base_size=13)+theme(legend.title = element_blank(),#plot.title = element_text(vjust = - 10)
                        )
                        #geom_text()    
        #
                }
         }
         #===------------------need to take care of this manually
         g4.ti[[1]]<-g4.ti[[1]]+geom_text(data=data.frame(X=c(-1.8),Y=c(0.6)),
                    color="black", size=5, aes(label=c("BM+IgM")))
         g4.ti[[2]]<-g4.ti[[2]]+geom_text(data=data.frame(X=c(1),Y=c(1.3)),
                    color="black", size=5, aes(label=c("BM+IgG")))
        g4.ti[[3]]<-g4.ti[[3]]+geom_text(data=data.frame(X=c(-2.2),Y=c(-1.4)),
                    color="black", size=5, aes(label=c("Spleen+IgM"))) 
        g4.ti[[4]]<-g4.ti[[4]]+geom_text(data=data.frame(X=c(1.7),Y=c(-0.2)),
                    color="black", size=5, aes(label=c("Spleen+IgG")))            
                    
           tiff(file=here(output.dir,
                "suppfigure6.cluster.tiff"), 
            width=600, height=1000)       
 
            ggarrange(g1, ggarrange(g4.ti[[1]],g4.ti[[2]],g4.ti[[3]],g4.ti[[4]], nrow=2, ncol=2, common.legend = TRUE, legend="bottom")
                ,ncol=1, nrow=2, labels=c("A","B"), widths=c(1.5,2))
                              
dev.off()         
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
                    geom_point(aes(shape=group), size=6)+ #xlim(c(-1.5, 1.5))+
                    geom_hline(yintercept=0, linetype="dashed", 
                color = "black", size=0.5)+
                geom_vline(xintercept=0, linetype="dashed", 
                color = "black", size=0.5)+
                    geom_text_repel(
                                label=rownames(d.quali), 
                                #nudge_x = 0.01, nudge_y = 0.1, 
                                 size=5#,repel=T
                              )+theme(legend.title = element_blank())
 
### plot the sub fig E
load(file=here(data.dir,
    "PCA_scree_plot.RData"))#, g1,g2,g3,g4 )#,d.out is loaded in above. the file is saved in analysis_v1.6_MFA.r    

#plot the sub fig F for doing 
load(file=here(data.dir,"fig6_PCA_partial.RData"))#,d.quali.all.two.line, d.quali ) 
temp<-d.quali.all.two.line[d.quali.all.two.line$group=="Isotype"|d.quali.all.two.line$group=="Tissue",]
k<-ggplot(data=temp)+
            geom_line(aes(x=Dim.1, y=Dim.2, colour=variable, group=groupVar, linetype=variable), size=1.5)+ 
            guides( shape=guide_legend(title="Group"))+
            geom_point(data=d.quali,aes(x=Dim.1, y=Dim.2, shape=factors),  size=5.5)+
            scale_shape_manual(values=c(0,1,18,15,16,17,9,10))+
            scale_linetype_manual(values=c("solid","dashed","twodash", "F1","longdash"))+
            #scale_color_manual(values=c("red","red","red","red","red"))+
            #geom_point(aes(x=Dim.1, y=Dim.2, colour=variable, group=groupVar), size=2.5, shape=0)+
            #geom_text(data=d.quali.lab,
            #        color="black", size=8, aes(label=rownames(d.quali.lab),x=Dim.1, y=Dim.2))+
            facet_wrap(.~group, ncol=2)+
            theme_bw(base_size=13)+theme(legend.position="right") 
 tiff(file=here(output.dir,"figure6.plot.v2.0.tiff"), 
        width=1000, height=1200)       
 ggarrange(  #first level
            #drawing pc1 and 2 for doing isotye and tissue 
            ggarrange(ggarrange(ind1, ind2, ncol=1, nrow=2, labels=c("A","B")), ind3, ncol=2, nrow =1, widths=c(1.5,2)
                    , labels=c("","C")),
            #draw by isotype and tissue 
            
            ind1.suppl.m,
            ggarrange(g1, #ggarrange(g4.ti[[1]],g4.ti[[2]],g4.ti[[3]],g4.ti[[4]], nrow=2, ncol=2, common.legend = TRUE, legend="bottom")
                                    k
                ,ncol=2, nrow=1, labels=c("E","F"), widths=c(0.5,1.1)),
                        ncol=1, nrow=3, heights=c(2,1.5,2), labels=c("","D","")
 )       
dev.off()



 

 tiff(file=here(output.dir,"suppfigure6.partial_immune.tiff"), 
        width=700, height=500)       
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
            theme_bw(base_size=25)+theme(legend.position="right")
            
                              
dev.off()

        #####################DONE ********************
        

