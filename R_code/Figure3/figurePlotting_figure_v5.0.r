#R code to plot the figure 3 in manuscript, 
#   Please run the below script before this
#    DataAnalysis_v1.0_MutationFreq_figure.R file.
#
## Make sure you have read ./ReadMe.txt for instructions
# of running scripts and preparing data.
#
library(ggplot2)
library(ggpubr)
library(here)

#
output.dir<-"R_code/Figure3"
data.dir<-"Data"
data.figure3.dir<-"Data/Figure3"
#load tissue difference mut freq
load(file=here(data.figure3.dir,"mutFreq_tissueOnly_fig.RData"))# mv, m1, mu.means.np 
                                                                                             # saved in DataAnalysis_v1.0.MutationFreq_fiugre.R
                                                                            # saved in the DataAnalysis_v1.0.MutationFreq_fiugre.R
load(file=here(data.figure3.dir,"m3wayIso_split.RData")) # load m3wayIso.pbs, m3wayIso.imm, mv.split

 load(file=here(data.figure3.dir,"m3way_treatment.RData")) # load mv, m3way, mu.means.np, 
                                                                            # saved in the DataAnalysis_v1.0.MutationFreq_fiugre.R
 
 load(file=here(data.figure3.dir,"mutationDifference_fig.RData" ))# md.line, x and w. 
                                                                                            #saved in DataAnalysis_v1.0.MutationFreq_fiugre.R
    md.bar<-ggplot(data=x, aes(y=the.emmean, x=treatment, fill=isotype, group=isotype  ))+
                                        ylab("Mu Freq Difference\n (Bone Marrow - Spleen)")+#geom_point(aes(), size=4)+
                                       geom_col(size=1.5,position=position_dodge(.5), width=0.5)+
                                       scale_x_discrete(expand=expansion(mult=c(0.0,0.0), add=c(.95,0.4)))+
                                       ylim(-0.0005, 0.007)+
                                       theme (text=element_text(size=16),axis.title.x = element_blank())+
                                      # geom_errorbar(aes(ymin=min_y, ymax=max_y), width=.2)+
                                       geom_segment(aes(x=0.2, xend=4.4, y=0.0, yend=0),size=0.8, colour="grey", linetype=2)+
                                       geom_text(#data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0.0005, treatment="OVA", isotype="IgG"),
                                        aes(x=0.05, y=0.0002,
                                            label=c("Bone Marrow")), hjust=0,size=5,colour="red")+
                                     geom_text(#data=data.frame(lower.CL=0, upper.CL=0, the.emmean=0.0005, treatment="OVA", isotype="IgG"),
                                        aes(x=0.05, y=-0.0002,
                                            label=c("Spleen")), hjust=0,size=5, colour="red")+
                                    geom_segment(aes(x=0.7, xend=0.7, y=0.0001, yend=0.0005),
                                        arrow = arrow(length = unit(0.2, "cm"),type="closed"), size=1, colour="red", linetype=1)+
                                    geom_segment(aes(x=0.7, xend=0.7, y=-0.0001, yend=-0.0005),
                                        arrow = arrow(length = unit(0.2, "cm"),type="closed"), size=1, colour="red", linetype=1)
#load(file="mutationDifferenceWOPBS_fig.RData" )# md2, x.2 and w. 
                                                                                            #saved in DataAnalysis_v1.0.MutationFreq_fiugre.R
                                                                                            #now in this one we have the difference plot excluding PBS
load(file=here(data.figure3.dir,"m3way_tisIso.RData")) #m3wayTis2, m3wayIso, also need the data frame mv that is save in the above treatment fig and data frame.
                                                                    #also  m3wayTis.pbs, mv.split.tis,

 png(file=here(output.dir,"figure2_v5.2.png"), width=1050, height=750)
 ggarrange(
                        ggarrange(ggarrange(m3wayIso.pbs, m3wayTis.pbs,
                        nrow=1, ncol=2,widths=c(1,1),labels=c("A","B")), m3way,  
                        nrow=1, ncol=2, widths=c(1.5,1.5), labels=c("","C")), 
                        #cdr3, 
                        ggarrange(NULL, md.line, NULL,nrow=1, ncol=3, labels=c("", "D",""),widths=c(1,4,1) ),  
                       nrow=2, ncol=1
                        )
dev.off()