#figure 2 plotting module. in here we load the data from previous work and draw the figure.
#R
## Make sure you have read ./ReadMe.txt for instructions
# of running scripts and preparing data.
library(plyr)
library(compositions)
library(ggplot2)
library(MASS)
library(ggpubr)
library(here)
library(ggfortify)

#   setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/figure4/noZeroMu/")
   #setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/figure4")
   #load for pca loadings
   data.dir<-"Data/Figure2"
   load(file=here(data.dir,"pca_loading.RData"))# saved in DataAnalysis_v2.0_geneUsage_compistion_figure4.r
                                                        #loaded with : g, g2, h, J, pca.VGene.clrd, cn, dt1

   #ggbiplot(pca.VGene.clrd)
h<-autoplot(pca.VGene.clrd,x=1,y=2,data=dt.clr.all, colour="treatment",shape="isotype",
				label=F ,frame=F, size=4,#frame.type="t", size=5
			loadings=F, loadings.label=F, loadings.label.size=4, loadings.colour="orange", loadings.label.colour="orange"
		) + scale_colour_discrete(name="Immunization")+
		#scale_colour_manual(name="treatment", values= c("black","pink","forest green", "red3" ))   +# c("forest green", "red3", "dark blue"))+
  ggtitle("PCA VGene Usage")+labs(y="PC2", x="PC1")+
  theme_bw(base_size=13) +
  theme(legend.position = c(0.01,0.45),legend.justification=c(0,1),plot.title = element_text(margin = margin(b = -20))) 


    #load for PC1 drawing
    load(file=here(data.dir,"figure4_pc1_draw.RData"))#saved in DataAnalysis_v2.0_geneUsage_compistion_figure4
                                                          #loaded with : dt.mod, pc1.tr)
    #load for PC1 lines
    load(file=here(data.dir,"figure4_pc1_lines.RData"));#saved in the same module as in above
                                #loaded with :, xt2 , pc1.line , trendline);   
    
    
    #load for PC2 drawing 
    load(file=here(data.dir,"figure4_pc2_draw.RData"));# saved in the same module as in above,
                        #loaded with: dt.mod.pc2, tr.pc2)
    #load for PC2 lines 
    load(file=here(data.dir,"figure4_pc2_lines.RData")); # saved in the same module as in above,
                            #loaded with :pc2.line, xt.pc2, trendline.pc2)
    

    load(file=here(data.dir,"figure4_cronbachAlphaValues.RData"))    #data were calculated and saved PC_consistency_CronbachA.R
                                                                                                                ### loaded with alpha 
    #add alpha to the dataframe
    dt1$alpha<-alpha[1:dim(dt1)[1]]
    cn<- ggplot(data=dt1, aes(x=ind, y=con.ac))+
            geom_bar(data=dt1, aes(x=ind,y=con), stat="identity", fill="turquoise2")+
            scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), limits=c(0,1))+
            geom_point()+ylab("Contribution (%)")+#ylim(c(0,1))+
            xlab("PCs")+ggtitle("PC contributions")+theme_bw(base_size=14)+
            geom_line(linetype=2)+
            geom_point(aes(x=ind, y=alpha), colour="grey")+
            geom_line(aes(x=ind, y=alpha), colour="grey",linetype=2)+
            geom_segment(x=6,xend=6, y=0,yend=1.1, colour="red", linetype=3 )+
            geom_segment(x=0,xend=50, y=0.05,yend=0.05, colour="red", linetype=3 )+
            annotate(geom="text", label=c("Accumulative", "Individual","Cronbach Alpha"), x=c(25,25, 23),y=c(0.9, 0.02,0.5), colour="red",size=7 )+
            theme(plot.title= element_text(margin = margin(b = -20), hjust=0.01))
	#   lines(c(0,64),c(0.015,0.015), lty=2, lwd=2,col=3)
    


tiff(here("Output","figure4_3.tiff"), width=1200, height=1400)
ggarrange(
#level 1 for 
     ggarrange(  cn, g,  h,
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
    ncol = 1, nrow = 3, heights=c(2,2.8,2.8)#, labels=c("A","D","F","H","J"),
    )
dev.off()
