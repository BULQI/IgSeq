#R code to plot the supplement figure to show the clone size distribution
#reading data from previously saved R data
#
#  prepare the data for figure5 plotting.
#
library(ggplot2)
library(here)
 library(ggpubr)

#read data.
data.dir<-"Data/Figure5"
output.dir<-"R_code/Figure5"
#setwd("/home/feng/Feng/LAB/Wetzler_PorinB/wl03R2/newPipeline/")
load(here(data.dir,"clones.df.RData"))#containing "BM.df.clones" "conditions"   "SP.df.clones" (lists)
mouse<-13
#get SP1 
SP.IgG<-SP.df.clones[["IgG"]]
SP.IgG<-SP.IgG[SP.IgG$sampleName==paste0("SP",mouse),]
SP.IgM<-SP.df.clones[["IgM"]]
SP.IgM<-SP.IgM[SP.IgM$sampleName==paste0("SP",mouse),]
BM.IgG<-BM.df.clones[["IgG"]]
BM.IgG<-BM.IgG[BM.IgG$sampleName==paste0("MB",mouse),]
BM.IgM<-BM.df.clones[["IgM"]]
BM.IgM<-BM.IgM[BM.IgM$sampleName==paste0("MB",mouse),]

min_seqNum<-min(sum(SP.IgG$X.Members), sum(SP.IgM$X.Members),
                        sum(BM.IgG$X.Members),sum(BM.IgG$X.Members))
#Scale the clone size.
SP.IgG$X.Members<-SP.IgG$X.Members/sum(SP.IgG$X.Member)*min_seqNum
SP.IgM$X.Members<-SP.IgM$X.Members/sum(SP.IgM$X.Member)*min_seqNum
BM.IgG$X.Members<-BM.IgG$X.Members/sum(BM.IgG$X.Member)*min_seqNum
BM.IgM$X.Members<-BM.IgM$X.Members/sum(BM.IgM$X.Member)*min_seqNum


#brks<-c(0,1,2,5,10,20,40)  #for smaller sequences 600+ (SP13) and 
#brks<-c(0,1,2,5,10,20,400) #800 + (SP2)
#brks<-c(0,3,6,15,30,60,600)#for SP1
brks<-c(0,2,4,10,20,40,400) #800 + (SP3)
size_lab<-c(paste0("#<=", brks[c(-1,-length(brks))]),paste0("# >",brks[length(brks)-1]))
G<-hist(SP.IgG$X.Members, breaks=brks)
GS<-data.frame(counts=G$counts, breaks=G$breaks[-1], 
        size=size_lab)
GS$size<-factor(GS$size, levels=size_lab)
M<-hist(SP.IgM$X.Members, breaks=brks)
MS<-data.frame(counts=M$counts, breaks=M$breaks[-1], 
        size=size_lab)
MS$size<-factor(MS$size, levels=size_lab)


p1<-ggplot(GS, aes(x = "", y =counts, fill =size)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  #geom_text(aes(y = lab.ypos, label = labels), color = "white")+
  scale_fill_manual(values = rainbow(length(G$breaks))) +labs(title="Spleen IgG")+
  theme_void()+theme(plot.title = element_text(hjust = 0.5),legend.text=element_text(size=12))+
  guides(fill=guide_legend(title="Clone size"))
  
p2<-ggplot(MS, aes(x = "", y =counts, fill =size)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  #geom_text(aes(y = lab.ypos, label = labels), color = "white")+
  scale_fill_manual(values = rainbow(length(M$breaks))) +labs(title="Spleen IgM")+
  theme_void()  +theme(plot.title = element_text(hjust = 0.5),legend.text=element_text(size=12))+
  guides(fill=guide_legend(title="Clone size"))

BG<-hist(BM.IgG$X.Members, breaks=brks)
BGS<-data.frame(counts=BG$counts, breaks=BG$breaks[-1], 
        size=size_lab)
BGS$size<-factor(BGS$size, levels=size_lab)
BM<-hist(BM.IgM$X.Members, breaks=brks)
BMS<-data.frame(counts=BM$counts, breaks=BM$breaks[-1], 
        size=size_lab)
BMS$size<-factor(BMS$size, levels=size_lab)


p3<-ggplot(BGS, aes(x = "", y =counts, fill =size)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  #geom_text(aes(y = lab.ypos, label = labels), color = "white")+
  scale_fill_manual(values = rainbow(length(BG$breaks))) +labs(title="Bone Marrow IgG")+
  theme_void()+theme(plot.title = element_text(hjust = 0.5),legend.text=element_text(size=12))+
  guides(fill=guide_legend(title="Clone size"))
  
p4<-ggplot(BMS, aes(x = "", y =counts, fill =size)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  #geom_text(aes(y = lab.ypos, label = labels), color = "white")+
  scale_fill_manual(values = rainbow(length(BM$breaks))) +labs(title="Bone Marrow IgM")+
  theme_void()  +theme(plot.title = element_text(hjust = 0.5),legend.text=element_text(size=12))+
  guides(fill=guide_legend(title="Clone size"))
   
   
 

 #setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/figure3/")
 save(p1,p2,p3,p4, GS, MS, BGS, BMS, 
    file=here(data.dir,"cloneSize_distribution13.RData"))

tiff(file=here(output.dir,paste0("cloneSizeSpleen_mouse", mouse, ".tiff")), 
  width=650, height=630)
ggarrange(p1, p2,  p3,p4,
                          ncol = 2, nrow = 2,
                          #heights = c(2, 2)),  
          #ncol = 1, nrow = 2, 
           common.legend=T, legend="bottom"#, labels=c("OVA+Alum","","","")
          #heights = c(2, 3)
          )
dev.off()
