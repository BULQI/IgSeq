#R code to determine the diversity HILL estimation.
#  #this is based on Diversity_sampleDepth_simulation.R
#  this one we will do subsampling using the smallest number of 
#   sequences. so we can compare the diversity profile.
# this can be run on linux.
#===================old notes down++++++++++++++++++++=
# 1/6/2022, we do subsampling to show the effect of sample depth on the diversity HILL numbers.
#    we also try to show the D(infinity) which is also called Berger Parker Index
#

#library(iNEXT)
library(SpadeR)
library(here)
#library(frLib)

###define functions used for doing the job
##
#turn abundance vector into accumulative sum vector 
accumulativeSum<-function(clones)
{
    #temp<-clones
    for(i in 2:length(clones))
    {
            clones[i]<-clones[i]+clones[i-1]
    }
    return(clones)
}


#determine the clone abunance based on subsampled sequences
#'@param clones this is the clone abundance vector as c(1,2,10,1,1,20,). each item is for one clone and the number is 
#'          the number of sequences in the specific clone
#'@param num.samples, the number of sequence sampled. we assume the total number of sequences is the sum(clones)
#'@values return a vector of number which is the abunance of clones in a format as the input clones
#          
determineClones<-function(clones, num.samples)
{
    #set.seed(1) #for doing samplings
    if(num.samples>=sum(clones))
    {
        warning("the total number of subsampled is more then the total number of clones specified!\n return the original sequence abundance")
        return(clones)
    }
    #doing subsampling first using sample
    f_totalnum<-sum(clones)

    f_subsamples<-sample(1:f_totalnum, num.samples)
    #cat("subsamples:",f_subsamples,"\n")
    
    #now determine for the abundence for 
    clones_acc<-accumulativeSum(clones)
    #cat("abunance:",clones_acc,"\n")
    
    #now calling histgram to do summary of each "clone"
    f_abundence<-hist(f_subsamples, breaks=c(0,clones_acc),plot=F)
    return (f_abundence$counts)
}

#'@description we are doing subsampling of sequences (no colons) to estimate the diversity with repeats
#'@param abundance the abundance vector each item is for one clone and number is the abundance
#'@param sub.sample.num number of sequences to do subsampling.
#'@param num.rep number of repeats to do subsampling.
#'@param q the order vector that we want to estimate diversity /hill numbers.
#'@values return a matrix holding the results each column is for one run.
diversity_subSample<-function(abundance, sub.sample.num, num.rep, q=seq(0,4,0.25))
{
    div.sub<-matrix(0, nrow=length(q), ncol=num.rep)
    rownames(div.sub)<-q
    sub.abundance.matrix<-matrix(0, nrow=length(abundance),ncol=num.rep)
    i<-1
    while(i <=num.rep)
    {
        cat("doing round i:", i, "........\n")
        x<-determineClones(abundance, sub.sample.num)
        sub.abundance.matrix[,i]<-x
        #do diversity
        temp<-NULL
        skip_to_next<-FALSE
        tryCatch(temp<-Diversity(x,datatype="abundance",q=q),
                        error=function(e){print("inside error block,set flag"); 
                            skip_to_next <<- TRUE;
                            print(e);})
        if(is.na(sum(temp$Hill_numbers[,2])))
        {
            cat("\t NaN error !!!.....skip \n")
            next
        }
        if(skip_to_next)
            {
                cat("\t found error try skip\n")
                next
            }
            cat("\tno skipping, a good one. updating.........")
        div.sub[,i]=temp$Hill_numbers[,2]
        i<-i+1
    }
    return (div.sub)
}


#read the data first

#read the saved data, which were saved by DataAnalysis_v1.0_clone_data.R
#read the files 15 mice (5 group x 3 mice each and two tissue)
#setwd("/home/feng/Feng/LAB/Wetzler_PorinB/wl03R2/newPipeline/")
#setwd(".")
load(here("Data/Figure5","clones.Rdata")) # loaded with BM.clones, SP.clones and conditions, (list(IgG and IgM), list (IgG+IgM), data.frame)
                                               #clone.inc data , clone incidence freq data (not count data yet)

isotypes<-c("IgG","IgM")
tissue<-c("Bone Marrow", "Spleen")
#first determine the smallest number of sequences in the samples
subnum<-1E9
sub.sampleName<-"none"
clones.list<-list("Bone Marrow"=BM.clones, Spleen=SP.clones )

#determine number of sequence to subsample 
cat("detere the smallest number of sequences to do subsampling..........\n")
for(k in 1:length(clones.list))
{
    clones.tissue<-clones.list[[k]]
    for(i in 1: length(clones.tissue)){
        clones.isotype<-clones.tissue[[isotypes[i] ]]
        for(j in 1:length(clones.isotype))
        {
            clones<-clones.isotype[[j ]]
            cat("doing.....sample=",unique(clones$sampleName),"; tissue:", names(clones.list)[k], 
                                            ";isotype=", isotypes[i],"\n")
            temp.seqNum<-sum(clones$X.Members)
            if(j==15&(i==2|i==1)&k==1)
            {
                cat("skip.....sample=",unique(clones$sampleName),"; tissue:", names(clones.list)[k], 
                                            ";isotype=", isotypes[i],"\n")
                next;
            }
            
            if(temp.seqNum<subnum)
            {
                subnum<-temp.seqNum
                sub.sampleName<-unique(clones$sampleName)
                sub.sampleName<-paste0(names(clones.list)[k], isotypes[i],sub.sampleName)
            }
            cat("\tsubnum:", subnum,"\n")
        }
    }
}
cat("done ............\n")
cat("\t smallest one:", subnum , "in sample ", sub.sampleName,"\n")

# we do abunance data to subsample  and deteremine the diversity 

cat("do subsampling to deteremine diversity ..........\n")
set.seed(1)
#now let run subsampling, with 1/100 number of sequences
num.rep<-100; subnum<-1000 #<-----------------------------
q<-seq(0,4,0.25)
div.tbl<-data.frame()
for(k in 1:length(clones.list))
{
    cat("tissue:", names(clones.list)[k],"\n")
    clones.tissue<-clones.list[[k]]
    for(i in 1: length(clones.tissue)){
        cat("\tisotype:", isotypes[i],"\n")
        clones.isotype<-clones.tissue[[isotypes[i]]]
        for(j in 1:length(clones.isotype))
        {
            
            clones<-clones.isotype[[j]]
            cat("\t\tsample:", unique(clones$sampleName),"\n" )
            if(unique(clones$sampleName)=="SP7"|unique(clones$sampleName)=="SP8"|unique(clones$sampleName)=="SP9"|
                                unique(clones$sampleName)=="MB7"|unique(clones$sampleName)=="MB8"|unique(clones$sampleName)=="MB9")
            {
                cat("skip.....BM=",unique(clones$sampleName),"\n")
                next;
            }
            #temp.seqNum<-sum(clones$X.Members)
            temp.abundance<-clones$X.Members
            run.subnum<-subnum
            if(sum(temp.abundance)<subnum)
            {
                run.subnum<-sum(temp.abundance)
                cat("***calling on run subnum:", run.subnum,"\n")
            }
            temp.div<-diversity_subSample(abundance=temp.abundance, sub.sample.num=run.subnum,num.rep=num.rep, q=q)
            #now we need to do averaging.
            temp.mean<-apply(temp.div, 1,mean)
            #
            temp.tbl<-data.frame(q=q,ChaoJost=temp.mean, sampleName=clones$sampleName[1], tissue=names(clones.list)[k], 
                                                            isotype=isotypes[i])
            div.tbl<-rbind(div.tbl, temp.tbl)
            ReadMe.text<-paste0("done to do subsampling with 1000 sequences (with exceptions), q=seq(0,4,0.25), repeat 100x!!  sample :",
                            unique(clones$sampleName),"\n")
            save(file=here("Data/Figure5","diversity_subSample2.RData"), div.tbl, ReadMe.text)

        }
    }
}
#ReadMe.text<-"done to do subsampling with 1000 sequences (with exceptions), q=seq(0,4,0.25), repeat 100x!! all done!!!"
#save(file="diversity_subSample.RData", div.tbl, ReadMe.text)




