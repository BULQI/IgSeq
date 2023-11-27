#in this module we will write code to batch-calculate intra-clonal diversities of samples
##
# Make sure you have read ./ReadMe.txt for instructions
# of running scripts and preparing data.
# 
##library
library(plyr)
library(frLib)
library(rcPkg)
library(here)

data.dir<-"Data/Figure6"
#setwd("/projectnb/keplab/feng/Sequencing/130.64.74.72/200121-0384H_Feng_Feng/analysis_2/merge")
load(here(data.dir,"mutations_df.RData"));  
	#now we should have two objects/lists, df.mutations.bm and df.mutations.sp
	#in DataAnalysis_v1.0_clone_data.R
load(here(data.dir,"clones.df.RData"))
	#now we should have two lists, BM.clones and SP.clones
	#in DataAnalysis_v1.0_clone_data.R
load(here(data.dir,"cloneAssignment_df.RData"))
	#now we should have two data frames, df.cloneAssignment.sp and df.cloneAssignment.bm
	#in DataAnalysis_v1.0_clone_data.R
load(here(data.dir,"Ig.good.RecSum.df.RData"))# recombination data frame 
    #IG.RecSum.bm, IG.RecSum.sp
	#DataAnalysis_v1.0_geneUsage_data.R
load(here(data.dir,"df.vlength.RData"));
	#df.vlength.sp, df.vlength.bm 
	#in processAignmentLenOfV.R
	
# BM IgG
df.clones<-BM.df.clones[["IgG"]]
df.cloneAssign<-BM.df.cloneAssign[["IgG"]]
IG.RecSum<-IG.RecSum.bm[["IgG"]]
df.mutations<-BM.df.mutations[["IgG"]]
df.vlength<-BM.df.vlength[["IgG"]]

df.clone.test<-df.clones [,c("sampleName", "CloneID", "X.Members")]
 df.cloneAssign.test<-df.cloneAssign[, c("sampleName", "CloneID", "ReadID")]
df.mutations.test<-df.mutations[, c("ReadID", "Type", "Position")]
 IG.RecSum.test<-IG.RecSum[,c("UID", "X.VBases")]
 df.vlength.test<-df.vlength[, c("ReadID", "totalVBase", "discPosStart", "discPosEnd")]; 
 idis<-intraClonalDiversities_w(clones=df.clone.test, #clone summary, 
                       clone.assigns=df.cloneAssign.test, #clone assignments for sequences
                      mutations=df.mutations.test, #seq muations VRG
                       Ig.RecSums=IG.RecSum.test, #Ig Recsum for VBases
                       vlengths=df.vlength.test,  #total lengths of v 
                       mCloneSize=2, #have all sizes 
                       indel.penalty=1 #for penalty of indel.
                      ) #
                      
 idis.BM.IgG<-idis
 
 
 # BM IgM
df.clones<-BM.df.clones[["IgM"]]
df.cloneAssign<-BM.df.cloneAssign[["IgM"]]
IG.RecSum<-IG.RecSum.bm[["IgM"]]
df.mutations<-BM.df.mutations[["IgM"]]
df.vlength<-BM.df.vlength[["IgM"]]

df.clone.test<-df.clones [,c("sampleName", "CloneID", "X.Members")]
 df.cloneAssign.test<-df.cloneAssign[, c("sampleName", "CloneID", "ReadID")]
df.mutations.test<-df.mutations[, c("ReadID", "Type", "Position")]
 IG.RecSum.test<-IG.RecSum[,c("UID", "X.VBases")]
 df.vlength.test<-df.vlength[, c("ReadID", "totalVBase", "discPosStart", "discPosEnd")]; 
 idis<-intraClonalDiversities_w(clones=df.clone.test, #clone summary, 
                       clone.assigns=df.cloneAssign.test, #clone assignments for sequences
                      mutations=df.mutations.test, #seq muations VRG
                       Ig.RecSums=IG.RecSum.test, #Ig Recsum for VBases
                       vlengths=df.vlength.test,  #total lengths of v 
                       mCloneSize=2, #have all sizes 
                       indel.penalty=1 #for penalty of indel.
                      ) #
                      
 idis.BM.IgM<-idis
 #setwd("/projectnb/keplab/feng/Sequencing/130.64.74.72/200121-0384H_Feng_Feng/analysis_2/merge/IntraClonalDiversity")

 save( idis.BM.IgM, idis.BM.IgG, 
  file=here(data.dir,"intraClonalDiversity_BM.RData"))

# SP IgG
df.clones<-SP.df.clones[["IgG"]]
df.cloneAssign<-SP.df.cloneAssign[["IgG"]]
IG.RecSum<-IG.RecSum.sp[["IgG"]]
df.mutations<-SP.df.mutations[["IgG"]]
df.vlength<-SP.df.vlength[["IgG"]]

df.clone.test<-df.clones [,c("sampleName", "CloneID", "X.Members")]
 df.cloneAssign.test<-df.cloneAssign[, c("sampleName", "CloneID", "ReadID")]
df.mutations.test<-df.mutations[, c("ReadID", "Type", "Position")]
 IG.RecSum.test<-IG.RecSum[,c("UID", "X.VBases")]
 df.vlength.test<-df.vlength[, c("ReadID", "totalVBase", "discPosStart", "discPosEnd")]; 
 idis<-intraClonalDiversities_w(clones=df.clone.test, #clone summary, 
                       clone.assigns=df.cloneAssign.test, #clone assignments for sequences
                      mutations=df.mutations.test, #seq muations VRG
                       Ig.RecSums=IG.RecSum.test, #Ig Recsum for VBases
                       vlengths=df.vlength.test,  #total lengths of v 
                       mCloneSize=2, #have all sizes 
                       indel.penalty=1 #for penalty of indel.
                      ) #
                      
 idis.SP.IgG<-idis
 
 
 # SP IgM
df.clones<-SP.df.clones[["IgM"]]
df.cloneAssign<-SP.df.cloneAssign[["IgM"]]
IG.RecSum<-IG.RecSum.sp[["IgM"]]
df.mutations<-SP.df.mutations[["IgM"]]
df.vlength<-SP.df.vlength[["IgM"]]

df.clone.test<-df.clones [,c("sampleName", "CloneID", "X.Members")]
 df.cloneAssign.test<-df.cloneAssign[, c("sampleName", "CloneID", "ReadID")]
df.mutations.test<-df.mutations[, c("ReadID", "Type", "Position")]
 IG.RecSum.test<-IG.RecSum[,c("UID", "X.VBases")]
 df.vlength.test<-df.vlength[, c("ReadID", "totalVBase", "discPosStart", "discPosEnd")]; 
 idis<-intraClonalDiversities_w(clones=df.clone.test, #clone summary, 
                       clone.assigns=df.cloneAssign.test, #clone assignments for sequences
                      mutations=df.mutations.test, #seq muations VRG
                       Ig.RecSums=IG.RecSum.test, #Ig Recsum for VBases
                       vlengths=df.vlength.test,  #total lengths of v 
                       mCloneSize=2, #have all sizes 
                       indel.penalty=1 #for penalty of indel.
                      ) #
                      
 idis.SP.IgM<-idis
 
 #setwd("/projectnb/keplab/feng/Sequencing/130.64.74.72/200121-0384H_Feng_Feng/analysis_2/merge/IntraClonalDiversity")
save(idis.SP.IgM, idis.SP.IgG, 
  file=here(data.dir,"intraClonalDiversity_SP.RData"))

cat("Done saving the output.\n")
