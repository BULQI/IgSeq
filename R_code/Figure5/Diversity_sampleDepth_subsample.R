#R code to determine the sample depth on diversity estimation.
# Do subsampling to show the effect of sample depth on the diversity HILL numbers.
#    we also try to show the D(infinity) which is also called Berger Parker Index
#
#   Make sure you have read ./ReadMe.txt for 
#   preparing data input and running this module.
#
#library(iNEXT)
library(SpadeR)
library(frLib)
library(here)
data.dir<-"Data/Figure5"
output.dir<-"R_code/Figure5"

#read the data first

load(here(data.dir,"clones.Rdata")) # loaded with BM.clones, SP.clones and conditions, (list(IgG and IgM), list (IgG+IgM), data.frame)
                                               #clone.inc data , clone incidence freq data (not count data yet)
#now we need to turn the data into the ones that can be used by SpadeR and iNext.
# we do abunance data first. for spleen IgM
clones<-SP.clones[["IgM"]][[3]]

abundance<-clones$X.Members
div.true<-c()
 system.time({
        temp<-Diversity(abundance,datatype="abundance",q=seq(0,3, 0.25))
        div.true<-temp$Hill_numbers
})      
       
#now start calling to get Diversity


#now start subsampling, first what is the total number of samples
totalnum<-sum(abundance)

subnum<-totalnum/100
subsamples<-sample(1:totalnum, subnum)

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


set.seed(1)
#now let run subsampling, with 1/100 number of sequences
num.rep<-100
div.sub<-matrix(0, nrow=dim(div.true)[1], ncol=num.rep)
sub.abundance.matrix<-matrix(0, nrow=length(abundance),ncol=num.rep)
for(i in 1:num.rep)
{
    cat("doing round i:", i, "........\n")
    x<-determineClones(abundance, floor(sum(abundance)/100))
    sub.abundance.matrix[,i]<-x
    #do diversity
    temp<-Diversity(x,datatype="abundance",q=seq(0,3, 0.25))
    div.sub[,i]=temp$Hill_numbers[,2]
}

set.seed(100)
#subsampling 1/10x sequences
div.sub.10<-matrix(0, nrow=dim(div.true)[1], ncol=num.rep)
sub.abundance.matrix<-matrix(0, nrow=length(abundance),ncol=num.rep)
for(i in 1:num.rep)
{
    cat("doing round i:", i, "........\n")
    x<-determineClones(abundance, floor(sum(abundance)/10))
    sub.abundance.matrix[,i]<-x
    #do diversity
    temp<-Diversity(x,datatype="abundance",q=seq(0,3, 0.25))
    div.sub.10[,i]=temp$Hill_numbers[,2]
}
#setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/figure5/")
png(file=here(output.dir,"SampleSizeDiversity.png"), 
        width=650, height=650)
#plot the results
plot(div.true[,1], div.true[,2], log="y", main="Diversity with Hill Numbers", 
                    xlab="order (q)", ylab="Hill Numbers (Log)")
for(i in 1:num.rep)
{
    lines(div.true[,1],div.sub[,i], col="black",lty=2)
}
lines(div.true[,1], div.true[,2], col="red", lwd=2)

#plot the results
#plot(div.true[,1], div.true[,2], log="y")
for(i in 1:num.rep)
{
    lines(div.true[,1],div.sub.10[,i], col="light blue")
}
points(div.true[,1], div.true[,2], col="black", lwd=2)
lines(div.true[,1], div.true[,2], col="red", lwd=2)
legend(x=1.6, y=5e4, legend=c("True", "1/10 SubSampling", "1/100 SubSampling"), 
                                    pch=c(1,-1,-1), lty=c(1,1,2),col=c("red","light blue", "black"))
                                    
                                    
dev.off()
