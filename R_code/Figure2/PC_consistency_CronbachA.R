#PC_consistency_CronbachA.R
#R code to load the data and then calculate the cronbach alapha for each PC.

# Make sure to run the below script before this one
#           DataAnalysis_V2.0_geneUsage_compsition_figure4.r

#We did three things: Cronbach Alpha, split-half reliability on IGVH and samples.
#       split-half did not work, since it tests whether the variation depends on IGVH and samples. 
#                   of course it depends on them.!!!

library(here)
data.dir<-"Data/Figure2"
output.dir<-"R_code/Figure2"
load(here(data.dir,"pcaData_cronbachAlapha.RData"))    #data were saved DataAnalysis_V2.0_geneUsage_compsition_figure4_noZeroMu.r 
                    #we should have pca.VGene.clrd <-- the pca object by prcomp
                       # dt.clr.data <-original data scaled and center 
###
#start caculating Cronbach's alpha 
#loop to get alpha for each PC
num_PC<-dim(pca.VGene.clrd$rotation)[2]
rtn<-pca.VGene.clrd$rotation
num_sample<-dim(dt.clr.data)[1]
mdt<-as.matrix(dt.clr.data)
mdts<-scale(mdt, center=T, scale=T)
var.Sx<-rep(0, num_PC);
alpha<-rep(0, num_PC)
png(file=here(output.dir,"sample_withinVar_perPC.png"))
plot(x=c(0,46),y=c(-0.10,0.15),type="n", log="", xlab="sample ID", ylab="value(mean and var)")
sampleNames<-rownames(dt.clr.data)
#title(xlab=sampleNames)
for(i in 1:46)#1:num_PC) #loop through PCs
{
    #var.Sx<-0
    item.var<-rep(0,num_sample)
    item.mean<-rep(0,num_sample)
    for(j in 1:num_sample )
    {
        item<-mdts[j,]*rtn[,i]
        item.var[j]<-var(item)  
        item.mean[j]<-mean(item)
    }
    names(item.var)<-sampleNames
    lines(item.var, col=3, lty=2)
    lines(item.mean, type="l", col=2, lty=3 )
    points(sqrt(item.var)/item.mean, col=4)

    var.Sx[i]<-sum(item.var)
    alpha[i]<-num_sample/(num_sample-1)*(1-var.Sx[i]/var(pca.VGene.clrd$x[,i]))
}
lines(x=c(11.5,11.5),y=c(-0.1,0.15), lty=1, lwd=3, col=4)
lines(x=c(22.5,22.5),y=c(-0.1,0.15), lty=1, lwd=3, col=4)
lines(x=c(34.5,34.5),y=c(-0.1,0.15), lty=1, lwd=3, col=5)
lines(x=c(46,46),y=c(-0.1,0.15), lty=1, lwd=3, col=5)
legend(x=35, y=-0.072, legend=c("mean","variance"), col=c(2,3), lty=c(3,2),lwd=c(2,2))
text(x=18,y=0.05, labels="spleen IgM",col=4,cex=2)
text(x=39,y=0.05, labels="BM IgM",col=5,cex=2)

dev.off()
save(alpha, file=here(data.dir,"figure4_cronbachAlphaValues.RData"))
pdf(file=here(output.dir,"varOfEachItemScore.pdf"))
plot(var.Sx)
dev.off()
pdf(file=here(output.dir,"varOfEachPC.pdf"))
plot(pca.VGene.clrd$sdev)
dev.off()
pdf(file=here(output.dir,"cronbachAlpha.pdf"))
plot(y=c(0,1), x=c(1,45), type="n")
points(alpha[1:45])
dev.off()


#now we do split-half validation
#  note: it did not work since we have to few samples/too parameters.
#we run many times for two half of the prcomp scores and then calculate the correlation 
library(MASS)

#sample
runs<-1000
ncols <-dim(mdt)[2]
cs<-NULL
for ( i in 1:runs)
{
    cat("i: ",i,"........\n")
    #sample half 
    
    rows<-sample(ncols)
    temp<-mdt[,rows]
    #split into half
    t1<-temp[,c(1:(ncols/2))]
    t2<-temp[,c((ncols/2+1):ncols)]
    
    pr.t1<-prcomp(t1, scale=T)
    pr.t2<-prcomp(t2, scale=T)
    
    #corrlation
    cc<-cor(pr.t1$x, pr.t2$x)
    ccs<-rep(0,ncols/2)
    for(j in 1:(ncols/2)){
            ccs[j]<-cc[j,j]
    }
    cs<-rbind(cs,ccs)
}


plot(c(0,42),c(min(cs),max(cs)),type="n")

for(i in 1:(runs))
{
    points(1:(ncols/2), abs(cs[i,]))
    
}
pdf(file=here(output.dir,"splithalf_IGHV.pdf"))
plot(density(cs[,1]))
dev.off()
apply((cs), 2, median)


####now we do split-half on samples, so reliable on all instead of some subjects.
nsample<-c(3,2,3,3,   #IgG, spleen
                            3,2,3,3,  #IgM, 
                            3,3,3,3, #IgG BM
                            3,3,3,3    #IgM
                            )
cs<-NULL
for ( i in 1:runs)
{
    cat("i: ",i,"........\n")
    #sample half 
    t1<-NULL
    t2<-NULL
    for(k in 1:length(nsample))
    {
        rows<-sample(nsample[k])
        sm<-0
        if(k>1)
        {
            sm=sum(nsample[1:(k-1)])
        }
        temp<-mdt[,rows+sm]
        #split into half
        t1.temp<-temp[,1]
        t2.temp<-temp[,2]
        
        t1<-rbind(t1, t1.temp)
        t2<-rbind(t2, t2.temp)
    }
    pr.t1<-prcomp(t1, scale=T)
    pr.t2<-prcomp(t2, scale=T)
    
    #corrlation
    cc<-cor(pr.t1$x, pr.t2$x)
    ccs<-rep(0,length(nsample))
    for(j in 1:length(nsample)){
            ccs[j]<-cc[j,j]
    }
    cs<-rbind(cs,ccs)
}

pdf(file=here(output.dir,"splithalf_sample.pdf"))
plot(density(cs[,1]))
dev.off()
apply((cs), 2, median)
