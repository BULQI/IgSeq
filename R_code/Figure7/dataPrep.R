###R code to put together the data in order to do clustering and MFA
### ---- updated 1/20/22
#   add code to read latest data
#           diversity : subsampled normalized arcsin transformed, q=1:4
#            selection: 40 largest clones
####################
#-----updated 11/24/21
#   switch data source to the ones saving for plotting the manuscript figures . 
#           in hg/IgSeq_MS/manuscript/figure*
#           check individual section for the source of each piece of data.
#-----------------------
#  need to add the data to a big table.
#
#started 6/14/2021

#libraries
library(here)

#table layout
#    rows : by each sample,   BM1-IgM.... BM15-IgM, BM1-IgG,.... BM15-IgG, SP1-IgM, .... SP15-IgM, SP1-IgG, ... SP14-IgG
#     cols: mutation, usage, expansion and selection by groups
# 

#reading the data section
            ### IgHV gene usage. this is the old one ( BEFORE 11/24/21)
            #     first gene usage with clr transformed, read the mean mutation to use.
            #       the data were saved in DataAnalysis_V2.0_geneUsage_composition_figure4_lee.r   in /analysis_PorinB_2ndtry/
#             setwd("/home/feng/Windows/windowsD/feng/LAB/MSI/MousePorB_lee_IgSeq/WL03/WL03R2/newPipeline/analysis_PorinB_2ndtry/")
#               load(file="geneUsage_all_clr.RData") #load dt.clr with all samples including porin B group with BM 9 sample.
                                            #the data were clr transformed and contains condition for the groups
                                                #DataAnalysis_V2.0_geneUsage_composition_figure4_lee.r   
                                                
            #ending of old reading.
            ###########################
            #first gene usage with clr transformed
#            setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/figure4/")

#reading figure2
            load(file=here("Data/Figure2","geneUsage_all_clr.RData") )#load dt.clr with all samples including porin B group with BM 9 sample.
                                            #the data were clr transformed and contains condition for the groups
                                            #saved in DataAnalysis_v3.0_geneUsage_composition_figure4.r
                                            
                                            
geneUsage<-dt.clr                                
       colnames(geneUsage)<-sub(".", "-",colnames(geneUsage),fix=T)
                             
###now reading the mutation freq. mean mutation freq.
#    the old we get from separate files/runs, now we got the data from the data for manuscript figure.
    # setwd("/home/feng/Windows/windowsD/feng/LAB/MSI/MousePorB_lee_IgSeq/WL03/WL03R2/newPipeline/analysis_PorinB_2ndtry/")
            #the data were save in .DataAnalysis_V2.0_MutationFreq_figure_lee.r
    
#setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/figure4b/discnAllSize_W_zero/")    

load( file=here("Data/Figure3","muFreq_table_all.RData"))#muFreq.table, contains all the data include porin B BM 9 sample and also contain zero mutation seq    
                #the data were save in. DataAnalysis_V1.0_MutationFreq_figure.r, in /discnAllSize_W_zero/
                    # were all data, need to convert to mean and std
                    
     mu.means<-aggregate(MuFreq.table$MuFreq, by=list(tissue=MuFreq.table$tissue, 
                                            treatment=MuFreq.table$treatment, isotype=MuFreq.table$isotype
                                                     , sampleName=MuFreq.table$sampleName), mean)
mu.means$treatment<-factor(mu.means$treatment, levels=c("PBS","OVA","OVA+PorB", "OVA+CpG", "OVA+Alum"))
mu.means$isotype<-factor(mu.means$isotype, levels=c("IgG", "IgM"))
mu.means$tissue<-factor(mu.means$tissue, levels=c("Bone Marrow", "Spleen"))  
rownames(mu.means)<-paste0(mu.means$sampleName, "_",mu.means$isotype)
colnames(mu.means)[5]<-"mean"
mu.sd<-aggregate(MuFreq.table$MuFreq, by=list(tissue=MuFreq.table$tissue, 
                                            treatment=MuFreq.table$treatment, isotype=MuFreq.table$isotype
                                                     , sampleName=MuFreq.table$sampleName), sd)
  rownames(mu.sd)<-paste0(mu.sd$sampleName, "_",mu.sd$isotype)    

mu.means$sd<-mu.sd[rownames(mu.means),"x"]  

#now doing CDR3 lengths 
#setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/CDR3IgGSub/")
load(file=here("Data/Figure4","CDR3_length_data.RData"))#saved in CDR3IgGSub/recomb_v1.0_withZeroMu.r
                                                                                #two tables with  cdr3.mean, cdr3.sd)
    cdr3.mean$sd<-cdr3.sd$CDR3Length
    colnames(cdr3.mean)[6]<-"CDR3Sd"
rownames(cdr3.mean)<-paste0(cdr3.mean$sampleName, "_",cdr3.mean$isotype)
###now doing the diversity
###
##############++++++++++++++++++++++++++++++++++++
#      Originally we load the data from the below folder.    +
#                       newly updated to read from the manuscript figure folder.
#    ++++++++++++++++++++++++++++++++++++++++++++++
#setwd("/home/feng/Windows/windowsD/feng/LAB/MSI/MousePorB_lee_IgSeq/WL03/WL03R2/newPipeline/analysis_PorinB_2ndtry/")
#setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/figure5/")


    ########################33
    ##  this section is not necessary. we later switched to use subsampled arcsin transformed
    #################################
#load( file="Diversity_table_all.RData")#div.tbl, contains all the data include porin B BM 9 sample  
#    # from file:///home/feng/Feng/LAB/Wetzler_PorinB/wl03R2/newPipeline/Diversity/Diversity_Figure.r


#            div.m<-div.tbl
#            div.m$treatment<-factor(div.m$treatment, levels=c("PBS", "OVA", "OVA+PorB", "OVA+CpG", "OVA+Alum"))
#            div.m$tissue<-factor(div.m$tissue, levels=c("Spleen", "Bone Marrow"))
#            div.m$ChaoJost_t<-log(div.m$ChaoJost)
#            div.0<-div.m[div.m$q==0,c("ChaoJost_t","tissue", "isotype", "treatment", "ID")]
#            div.0.5<-div.m[div.m$q==0.5,c("ChaoJost_t","tissue", "isotype", "treatment", "ID")]
#                div.1<-div.m[div.m$q==1,c("ChaoJost_t","tissue", "isotype", "treatment", "ID")]
#                div.1.5<-div.m[div.m$q==1.5,c("ChaoJost_t","tissue", "isotype", "treatment", "ID")]
#                div.2<-div.m[div.m$q==2,c("ChaoJost_t","tissue", "isotype", "treatment", "ID")]
#                div.2.5<-div.m[div.m$q==2.5,c("ChaoJost_t","tissue", "isotype", "treatment", "ID")]
#                div.3<-div.m[div.m$q==3,c("ChaoJost_t","tissue", "isotype", "treatment", "ID")]
#                
#                #merge
#                div.qs<-merge(x=div.0, y=div.1, by=c("tissue", "isotype", "treatment", "ID"))
#                    names(div.qs)[c(5,6)]<-c("q0", "q1")
#                div.qs<-merge(x=div.qs, y=div.2, by=c("tissue", "isotype", "treatment", "ID"))
#                    names(div.qs)[c(7)]<-c("q2")
#                div.qs<-merge(x=div.qs, y=div.3, by=c("tissue", "isotype", "treatment", "ID"))
#                    names(div.qs)[c(8)]<-c("q3")    
#                div.qs<-merge(x=div.qs, y=div.0.5, by=c("tissue", "isotype", "treatment", "ID"))
#                    names(div.qs)[c(9)]<-c("q0.5")
#                div.qs<-merge(x=div.qs, y=div.1.5, by=c("tissue", "isotype", "treatment", "ID"))
#                    names(div.qs)[c(10)]<-c("q1.5")
#                    div.qs<-merge(x=div.qs, y=div.2.5, by=c("tissue", "isotype", "treatment", "ID"))
#                    names(div.qs)[c(11)]<-c("q2.5")
#                    rownames(div.qs)<-paste0(div.qs$ID, "_", div.qs$isotype)
#                    
        ##################################################
        #           END OF not necessary sectoin
        ###################################################

 ###we do another loading with subsampled and normalized and arcsin transformed
#        setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/figure5/subsample")
        load( file=here("Data/Figure5","Diversity_table_all_subsample_normArcsin.RData"))#, div.t3, saved in Diversity_Figure_subSample_v2.0_arcsin.r
            div.m<-div.t3
        div.m$treatment<-factor(div.m$treatment, levels=c("PBS", "OVA", "OVA+PorB", "OVA+CpG", "OVA+Alum"))
        div.m$tissue<-factor(div.m$tissue, levels=c("Spleen", "Bone Marrow"))
        div.m$ChaoJost_t<-div.m$trans
         div.0<-div.m[div.m$q==0,c("ChaoJost_t","tissue", "isotype", "treatment", "ID")]
          div.0.5<-div.m[div.m$q==0.5,c("ChaoJost_t","tissue", "isotype", "treatment", "ID")]
              div.1<-div.m[div.m$q==1,c("ChaoJost_t","tissue", "isotype", "treatment", "ID")]
              div.1.5<-div.m[div.m$q==1.5,c("ChaoJost_t","tissue", "isotype", "treatment", "ID")]
      div.2<-div.m[div.m$q==2,c("ChaoJost_t","tissue", "isotype", "treatment", "ID")]
      div.2.5<-div.m[div.m$q==2.5,c("ChaoJost_t","tissue", "isotype", "treatment", "ID")]
      div.3<-div.m[div.m$q==3,c("ChaoJost_t","tissue", "isotype", "treatment", "ID")]
      div.3.5<-div.m[div.m$q==3.5,c("ChaoJost_t","tissue", "isotype", "treatment", "ID")]
      div.4<-div.m[div.m$q==4,c("ChaoJost_t","tissue", "isotype", "treatment", "ID")]
      
      #merge
      div.qs.sna<-merge(x=div.0, y=div.1, by=c("tissue", "isotype", "treatment", "ID"))
        names(div.qs.sna)[c(5,6)]<-c("q0", "q1")
     div.qs.sna<-merge(x=div.qs.sna, y=div.2, by=c("tissue", "isotype", "treatment", "ID"))
        names(div.qs.sna)[c(7)]<-c("q2")
    div.qs.sna<-merge(x=div.qs.sna, y=div.3, by=c("tissue", "isotype", "treatment", "ID"))
        names(div.qs.sna)[c(8)]<-c("q3")    
    div.qs.sna<-merge(x=div.qs.sna, y=div.0.5, by=c("tissue", "isotype", "treatment", "ID"))
        names(div.qs.sna)[c(9)]<-c("q0.5")
     div.qs.sna<-merge(x=div.qs.sna, y=div.1.5, by=c("tissue", "isotype", "treatment", "ID"))
        names(div.qs.sna)[c(10)]<-c("q1.5")
        div.qs.sna<-merge(x=div.qs.sna, y=div.2.5, by=c("tissue", "isotype", "treatment", "ID"))
        names(div.qs.sna)[c(11)]<-c("q2.5")
    div.qs.sna<-merge(x=div.qs.sna, y=div.3.5, by=c("tissue", "isotype", "treatment", "ID"))
        names(div.qs.sna)[c(12)]<-c("q3.5")
    div.qs.sna<-merge(x=div.qs.sna, y=div.4, by=c("tissue", "isotype", "treatment", "ID"))
        names(div.qs.sna)[c(13)]<-c("q4")
        
        rownames(div.qs.sna)<-paste0(div.qs.sna$ID, "_", div.qs.sna$isotype)
####doing selection
###########
#==============================================
###originally the data were loaded from                            +
#++++++++++++++++++++++++++++++++++++++++++++
#setwd("/home/feng/Windows/windowsD/feng/LAB/MSI/MousePorB_lee_IgSeq/WL03/WL03R2/newPipeline/analysis_PorinB/")
                                    ###############<----note this is from  /analysis_PorinB folder !!!!
                                    #file:///home/feng/Windows/windowsD/feng/LAB/MSI/MousePorB_lee_IgSeq/WL03/WL03R2/newPipeline/analysis_PorinB/intraClonalDiversity_process_lee.R
                                    
                                    
                                    #now they were loaded from
#setwd("/home/feng/Feng/LAB/Wetzler_PorinB/wl03R2/newPipeline/IntraClonalDiversity")

#load("intraClonalDiversity_intraD_top15.RData")       #intraD and intraD.ilr
load(file=here("Data/Figure6","intraClonalDiversity_intraD_top40.RData"))       #intraD and intraD.ilr 
                            #data were saved by file:///home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/figure6/intraClonalDiversity_process_figure.r                             

#load(file="intraClonalDiversity_intraD_top20.RData") # intraD and intraD.ilr  <---old version
            #intraD.ilr contains % of selected clones
####
#setwd("/home/feng/Windows/windowsD/feng/LAB/MSI/MousePorB_lee_IgSeq/WL03/WL03R2/newPipeline/analysis_PorinB/")
#load( file="intraClonalDiversityDataArray_top20.RData")# load data.array and its name array
#setwd("/home/feng/Feng/LAB/Wetzler_PorinB/wl03R2/newPipeline/IntraClonalDiversity")
#load("intraClonalDiversityDataArray_top15.RData")#loading data with data.array , name.array and conditions
load(file=here("Data/Figure6","intraClonalDiversityDataArray_top40.RData"))#loading data with data.array , name.array and conditions
      #data.array i=1, BM IgM, i=2, BM IgG, i=3, SP IgM, i=4 SP IgG
                #saved in file:///home/feng/Feng/LAB/Wetzler_PorinB/wl03R2/newPipeline/IntraClonalDiversity/intraClonalDiversity_process.R
      
        ###############<----note this is from  /analysis_PorinB folder !!!!
        #this data.array is holding the original data of idi and mean mutation of clones (top 20 clones)
        #going through to get top10 clones for comparison
       sel.dt<-data.frame()
       size.topClone<-10
       isotype<-substr(name.array, nchar(name.array)-3, nchar(name.array)-1)
       for(i in 1:length(data.array))
        {
            temp<-data.array[[i]]
            sampleArray<-unique(temp$sampleName)
               for(j in sampleArray)
               {
                    tp.sample<-temp[temp$sampleName==j, ]
                    tp.sample<-tp.sample[order(tp.sample$xMembers, decreasing=T)[1:size.topClone],]
                    tp.sample$clone.order<-1:size.topClone
                    tp.sample$isotype<-isotype[i]
                    sel.dt<-rbind(sel.dt, tp.sample)
               }
               #sel.dt$isotype<-isotype[i]
            
        }     
         #remove NA
         sel.dt<-sel.dt[!is.na(sel.dt$sampleName),]

        #now turn it into sample vs clone 1 2.... 10 we will do clone idi, clone s.index, clone mean mutation and clone CDR3 length
        #
            sel.tbl<-data.frame()
            for (i in 1:size.topClone)
            {   
                temp<-sel.dt[sel.dt$clone.order==i,
                        c("sampleName", "cloneID", "idi","CDR3Length",  "S.Index", "MeanMuFreq","isotype")]
                 rownames(temp)<-paste0(temp$sampleName,temp$isotype)
                 colnames(temp)[c(3,4,5,6)]<-paste0(colnames(temp)[c(3,4,5,6)],".c",i)
                 if(i==1)
                 {
                    sel.tbl<-temp
                 }
                 else {
                    sel.tbl<-cbind(sel.tbl, temp[rownames(sel.tbl),c(3,4,5,6)])
                 }
            }
    #sort the columns to put things together
    sel.tbl<-sel.tbl[,order(colnames(sel.tbl))]
    sel.tbl<-cbind(sel.tbl,intraD.ilr[rownames(sel.tbl),])
    colnames(sel.tbl)[dim(sel.tbl)[2]]<-"selectionIndex"
    rownames(sel.tbl)<-paste0(sel.tbl$sampleName, "_", sel.tbl$isotype)
    
####now put them together in order to do MFA
d.all<-geneUsage
d.all<-cbind(d.all,mu.means[rownames(d.all),c("mean", "sd")])
#d.all<-cbind(d.all,div.qs[rownames(d.all),-c(1:4)])  #<-original diversity 
d.all<-cbind(d.all,div.qs.sna[rownames(d.all),-c(1:4)]) #<-subsampled, normalized and arcsin transformed.
d.all<-cbind(d.all, sel.tbl[rownames(d.all),-c(11,22,43)])
d.all<-cbind(d.all, cdr3.mean[rownames(d.all),-c(1:4)])

###save the data
#setwd("/home/feng/Windows/windowsD/feng/LAB/MSI/MousePorB_lee_IgSeq/WL03/WL03R2/newPipeline/clustering")
save(d.all,file=here("Data/Figure7","all_data.RData"))
