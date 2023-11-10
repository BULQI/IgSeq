#R code to do the clustering analysis.
#####-----7/18/2021
#           doing pc of gene usage of 2 maximum ones.
###                 see analysis_v1.0_cluster.r for first 4 pcs of gene usage
####
#--------------------------------------------------
##      updated   7/17/2021-----
##       add code to  cluster without PorinB group
###           copy the code from file:///home/feng/Windows/windowsD/feng/LAB/MSI/MousePorB_lee_IgSeq/WL03/WL03R2/newPipeline/clustering/clusterWithLimitGU/GU_PC4/analysis_v1.0_cluster.r
####
#--------------------------------------------------------------
###
##-------------------7/10/2021
#           adding code to limit the effects of gene usage.
##                  run pca first and then select first 4 to do MFA.
#--------------------------------------
####
#      --------7/9/2021----
#cluster using all mfa n_components=55, explaining all variations
# run cluster with on all pcs

#--------
#   update 6/24/2021
##    add umap. and PCA results for doing cluster.
################################
#reading data from pca processed data. (see analysis_v1.0_MFA.R)
#
#we first try gower with missing value and then PAM
 ########the file copied from PJmaglione fold of doing cluster.
 #          ----6/22/2021
 #########################
 library(FactoMineR)
 library(missMDA)
 library(factoextra)
 library(cluster)
 
 library(pheatmap)  #only if you want to draw heatmap
library(dplyr)
#library(reshape2)  #<-no need??
library(igraph)   # for plotting groups
#library(linkcomm)  
library(kernlab)   #<- for specc and spirals data

library(rARPACK)     #<-for getting the ***generalized*** eigen value and eigen vector 
library(ggpubr)
library(here)
 #read data.
#setwd("/home/feng/Windows/windowsD/feng/LAB/MSI/MousePorB_lee_IgSeq/WL03/WL03R2/newPipeline/clustering")

#
data.dir<-"Data/Figure7"
output.dir<-"R_code/Figure7"
load(here(data.dir,"all_data.RData")) #read in  data containing everthing.  d.all


#need to do PCA first on gene usage to reduce dimension of usage, otherwise it has too much impact on the output
#
d.all<-d.all[d.all$treatment!="OVA+PorB",]

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
gu_pcs<-3

#put them back
#d.all<-d.all[,-c(1:93)]
d.all<-d.all[,-c(1:85)]
d.all<-d.all[,-c(2,14)]#get rid of sd of mutation at 2, and cdr3sd at 12.
d.all<-cbind(pca.VGene.clrd$x[,c(1:gu_pcs)], d.all)
d.all.conds<-cbind(d.all, conds[, c("tissue","treatment","isotype")])
#d.all.conds<-d.all.conds[,-c(5,9:12)]
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

n_pcs<-15
d.out <- MFA(d.all.conds, group=group, type=c(rep("s", length(group)-3),"n","n","n"),ncp=n_pcs, 
                     name.group = c("GeneUsage", "Mutation", "Diversity", #"CDR3Length", "Dissimiliarity","meanMuFreq", "S_index", 
                     "Selection","CDR3Length","tissue", "treatment","isotype" ), num.group.sup=c(6,7,8),
                     graph=F)


#setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/clustering/")
#setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/clustering/pc3")                   
#load(file="d.out.pc2.RData") 
#load(file="d.out.pc3.RData") 
                        #loaded  d.out and d.all.conds
                        #d.out is the MFA output with 2 gene usage 
                        # d.all.conds contains original data 
                        #the data were saved by analysis_v1.6_MFA.r at ../MFA
n_pcs<-14
######clean up the data 
#x.clean<-d.all   #<-----
x.clean<-as.data.frame(d.out$ind$coord[,c(1:n_pcs)])

#gower distance

#x.dist<-daisy(x.clean, metric="euclidean", stand=F)
d_method<-"minkowski"#"manhattan"#"gower"#"manhattan"
pmink<-5
#x.dist<-daisy(x.clean, metric=d_method, stand=F)
x.dist<-dist(x.clean, method=d_method, p=pmink)

silhouette <- c()
#silhouette = c(silhouette, NA)
for(i in 2:30){
  pam_clusters = pam(as.matrix(x.dist),
                 diss = TRUE,
                 k = i)
  silhouette = c(silhouette ,pam_clusters$silinfo$avg.width)
}
#setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/clustering/")
#setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/clustering/pc3")                   
png(file=here(output.dir,"PAM_silhouette.png"), 
    width=1050, height=750)
plot(2:30, silhouette,
     xlab = "Clusters",
     ylab = "Silhouette Width")
lines(2:30, silhouette)
dev.off()

pam_german  = pam(as.matrix(x.dist), diss =TRUE, k =12)

####
x.is_ti<-paste0(conds$isotype,"_", conds$tissue)#apply(conds.complications,2,as.integer)
x.is_ti<-factor(x.is_ti)
#for(i in 1:dim(x.complications)[2])
#{
#    x.complications[is.na(x.complications[,i]),i]<-0
#    }

#x.complications.sum<-apply(x.complications, 1, sum)

#visualization
library(Rtsne)
library(ggplot2)
library(ggrepel)
library(umap)
#set.seed(3456)#11234)#12345) <--- this is the one using the first 4 dimensions only for umap
#gower distance
#x.dist<-daisy(x.clean, metric=d_method, stand=F)
                          #x.dist<-dist(x.clean, method="e") <----------
#set.seed(1234)#<---this is the one for the pca components of 35

######
#set.seed(5679) 
 #set.seed (4321) 
 #set.seed(2345)#
 #set.seed(12345)#
#
#set.seed(2)#<-good
#set.seed(100)
#set.seed(1)
#set.seed(2)
set.seed(4)
   #tsne_object <- Rtsne(x.dist, is_distance = TRUE, perplexity=5)
   umap_object<- umap(x.clean, scale=F#, init= "spectral"
        #, n_neigbors=10,min_dist=0.01, n_components=2 
      )
  # plot(umap_object$layout[,1], umap_object$layout[,2])
#save(tsne_object, file="tsne.RData")
#load(file="tsne.RData")
#tsne_df <-tsne_object$Y %>%
tsne_df <- umap_object$layout[,c(1,2)] %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_german$clustering))

  tsne_df<-cbind(tsne_df, conds)
  png(here(output.dir,"umap_withPAM_cluster.png"), 
    width=1100, height=700)
  
  g1<-ggplot(aes(x = X, y = Y), data = tsne_df) +
  geom_point(aes(color =conds$tissue, shape=conds$isotype),size=6 )+#,size=(x.complications.sum+1)*3)+
        ggtitle("Isotype+Tissue")+guides(col=guide_legend(title="Tissue",override.aes = list(size=6)),
         shape=guide_legend(title="Isotype",override.aes = list(size=6)))+
         theme(legend.title=element_text(size=6), 
    legend.text=element_text(size=6))+
    geom_text_repel(label=paste(rownames(x.clean),as.factor(tsne_df$cluster)))
  g1
  dev.off()
  
  #by tissue 
  g1.it<-ggplot(aes(x = X, y = Y), data = tsne_df) +
  geom_point(aes(color =conds$isotype, shape=conds$tissue),size=6 )+#,size=(x.complications.sum+1)*3)+
        ggtitle("Tissue+Isotype")+
        guides(col=guide_legend(title="Isotype (colour)",override.aes = list(size=6)),
         shape=guide_legend(title="Tissue(shape)",override.aes = list(size=6)))+
         theme(legend.title=element_text(size=10), legend.position="bottom",
                            legend.text=element_text(size=10))#+
    #geom_text_repel(label=paste(rownames(x.clean),as.factor(tsne_df$cluster)))
  g1.it


###using umap object as te coord
#x.dist<-daisy(as.data.frame(umap_object$layout), metric="euclidean", stand=F)
#now doing ap clustering, affinity propogation
library(apcluster)  #<-show ap cluster affinity propogation
#x.clean<-as.data.frame(umap_object$layout)
x.dist<-dist(x.clean[, 1:4], method=d_method,p=pmink)
#x.dist<-daisy(x.clean, metric=d_method, stand=F)
x.sim<-exp(-1*(as.matrix(x.dist)))
#x.sim<-negDistMat(r=2,x.clean)
#c<-apcluster(x.sim, details=T)#, K=5)
c<-apcluster(negDistMat(r=1), x.clean, details=T)#, K=5)
#c<-apclusterK(x.sim, details=T, K=11)
clusters<-rep(1,dim(x.sim)[1])
#clusters[c@clusters[[1]]]<-1

for(i in 1:length(c@clusters))
{
    clusters[c@clusters[[i]]]<-i
}

g2<-ggplot(aes(x = X, y = Y, label=paste0(rownames(tsne_df),"+",clusters)), data = tsne_df) +
  geom_point(aes(shape =conds$tissue, color=factor(clusters)),size=6 )+#, size=(x.complications.sum+1)*2)+
            ggtitle("APC Clustering")+guides(col="none",
         shape=guide_legend(title="Tissue",override.aes = list(size=6)))+
         theme(legend.title=element_text(size=10), legend.position="bottom",
    legend.text=element_text(size=10))+
            geom_text_repel()
png(here(output.dir,"umap_withAPC_cluster.png"), 
  width=1100, height=700)
g2
dev.off() 

#--------------------start doing spectral clustering RData
############to generate an affinity matrix/adjacency matrix based on S similarity matrix.
#'@description A function, make.affinity to compute a restricted (or filtered) “affinity” between vertices using k-nearest neighbors.
#'          the function basically go through the rows and pick/filter out the nearest neighbors
#'@param S similarity matrix
#'@param n.neightboors number of nearest neighbors
make.affinity <- function(S, n.neighboors=2) {
  N <- length(S[,1])

  if (n.neighboors >= N) {  # fully connected
    A <- S
  } else {
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighboors]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j] <- S[i,j]
        A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
      }
    }
  }
  A
}




#Calculate the affinity matrix, A, from S by applying k-nearest neighbors algorithm
            A <- make.affinity(x.sim,2)  # use 3 neighbors (includes self)   
                #NOTE: neighbor # is very important, #=4, it will fail!!! and also in this spiral case 
                #           epsilon neighborhood method won't work!!! it is hard to find a good neighborhood to make spiral connected.




        #++++++++++++++++++++++++++++++
        eigen_m <- eigs(A,1,which="SR")
        toadd = ceiling(abs(min(eigen_m$values)))
        A.norm<-A
          diag(A.norm) = diag(A.norm) + toadd
        
        normalizeKernel = function(K) {
          # from Shawe-Taylor & Cristianini's "Kernel Methods for Pattern Analysis", p113
          # original kernel matrix stored in variable K
          # output uses the same variable K
          # D is a diagonal matrix storing the inverse of the norms
          # Based on MATLAB script from: www.kernel-methods.net
          rnames = rownames(K)
          cnames = colnames(K)
          D = diag(1/sqrt(diag(K)))
          K = D %*% K %*% D
          rownames(K) = rnames
          colnames(K) = cnames
          return(K)
         }

          A.norm = normalizeKernel(A.norm)
          
          mytriangle <- function(coords, v=NULL, params) {
              vertex.color <- params("vertex", "color")
              if (length(vertex.color) != 1 && !is.null(v)) {
                vertex.color <- vertex.color[v]
              }
              vertex.size <- 1/200 * params("vertex", "size")
              if (length(vertex.size) != 1 && !is.null(v)) {
                vertex.size <- vertex.size[v]
              }

              symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
                      stars=cbind(vertex.size, vertex.size, vertex.size),
                      add=TRUE, inches=FALSE)
            }
            # clips as a circle
            add_shape("triangle", clip=shapes("circle")$clip,
                             plot=mytriangle)
          
           library(igraph)
           png(here(output.dir,"umap_withNetwork_cluster.png"), 
              width=1100, height=700)
            g = graph_from_adjacency_matrix(A.norm, mode = "undirected", weighted = TRUE, diag = FALSE)
             plot.igraph(g, edge.width = E(g)$weight/(max(E(g)$weight)/5), vertex.label =rownames(x.sim), vertex.label.color="red",
                                    vertex.size = 5, edge.curved = T)
             dev.off()
             
            #cluster size
            csize=9
            mc2_specc = specc(as.kernelMatrix(A.norm), centers = csize)
            v_small = (mc2_specc@.Data)
            
            v_small.list<-vector("list", length=csize)
            names(v_small)<-rownames(x.sim)

            for(i in 1:csize){
                v_small.list[[i]]<-which(v_small==i)
            }
#            v_small.list[[1]]<-which(v_small==1)
#            v_small.list[[2]]<-which(v_small==2)
#            v_small.list[[3]]<-which(v_small==3)
#            v_small.list[[4]]<-which(v_small==4)
#            v_small.list[[5]]<-which(v_small==5)
            
            #reassign the clusters
             
             #adding triangle
             # triangle vertex shape
            
            
            g = graph_from_adjacency_matrix(A.norm, mode = "undirected", weighted = TRUE, diag = FALSE)
            shapes<-c("circle","triangle", "square", "sphere", "sphere")
            png(here(output.dir,"umap_withSpectral_cluster_1.png"), 
                width=1100, height=700)   
            h<- plot.igraph(g, edge.width = E(g)$weight/(max(E(g)$weight)/5), edge.color="black",
                    vertex.label =v_small, 
                    vertex.label.color="red",vertex.label.dist=1.5, vertex.label.font=4,vertex.label.cex=1.2,
                    vertex.size = as.integer(conds$isotype)*5, 
                    edge.curved = T, vertex.shape=shapes[as.integer(conds$treatment)],#[as.integer(factor(paste0(conds$tissue,"+",conds$isotype)))],
                    vertex.color=as.integer(conds$tissue)
                )
                lines(x=c(-1.0,2.0),y=c(-0.3,-0.3),lty=2,lwd=3) #<-to separate isotypes
            dev.off()
            
            
            g3<-ggplot(aes(x = X, y = Y), data = tsne_df) +
  geom_point(aes(shape=factor(paste0(conds$tissue,"+",conds$isotype)), color=factor(v_small)),size=10)+#, size=(x.complications.sum+1)*2)+
                ggtitle("spectral clustering")+guides(col=guide_legend(title="Cluster",override.aes = list(size=6)),
         shape=guide_legend(title="Tissue+Isotype",override.aes = list(size=6)))+
         theme(legend.title=element_text(size=6), 
                                            legend.text=element_text(size=6))+
                geom_text_repel(label=paste(rownames(x.clean),as.factor(v_small)))
                
                
png(here(output.dir,"umap_withSpectral_cluster2.png"), 
    width=1100, height=700)   
   #ggarrange(h, g3,ncol=1, nrow=2) 
    g3
   dev.off()
   
   ###drawing the bone marrow treatment without grey dots as background
     #v_small<-clusters
     
     cols_n<-colorspace::rainbow_hcl(length(unique(v_small[conds$tissue=="Bone Marrow"])), c=105,l=65)
  #        cols_n<-colorspace::heat_hcl(length(unique(v_small)))#, c.=50,l=50, h=20)
     cols<-rep("grey", length(v_small))
     cols[conds$tissue=="Bone Marrow"]<-cols_n[factor(v_small[conds$tissue=="Bone Marrow"])]
   
 
     g3.grey<-ggplot(aes(x = X, y = Y), data = tsne_df) +
            geom_point(aes(shape=factor(conds$treatment)), colour=cols, size=6)+#, size=(x.complications.sum+1)*2)+
                ggtitle("Clustering")+guides(col=guide_legend(title="Cluster",override.aes = list(size=6)),
                 shape=guide_legend(title="treatment",override.aes = list(size=6)))+
                 scale_shape_manual(values=c(19,17,15,7))+
                 theme(legend.title=element_text(size=6), 
                                                    legend.text=element_text(size=6))#+
               # geom_text_repel(label=paste(rownames(x.clean),as.factor(v_small)))
   save(g3.grey, g1.it, g1, conds,tsne_df, 
      cols,file=here(data.dir,"g3.grey.rdata") )  
   
   
   
                               #####################################
                               #     HDBScan  --- DO NOT RUN. for code demo
                               ####################################3
                               library(dbscan)
                               #xdat<-as.data.frame(umap_object$layout)
                               xdat<-as.data.frame(x.clean)
                               cl<-hdbscan(xdat,minPts=3)  # minPts the minimal # of points to form a cluster.
                               
                               #
                               kNNdistplot(xdat,k=4) #used to check for the eps. eps is used to determine the cut of hiearchical tree, or the level of the density for a real cluster.
                                                        # check dbscan manual for explanation about how to use kNNdistplot to determine eps
                               cld<-dbscan(x.dist, minPts=3, eps=1)
                             
                                            g4<-ggplot(aes(x = X, y = Y), data = tsne_df) +
                                                geom_point(aes(shape=factor(paste0(conds$tissue,"+",conds$isotype)), color=factor(cld$cluster+1L)),size=10)+#, size=(x.complications.sum+1)*2)+
                                                ggtitle("spectral clustering")+guides(col=guide_legend(title="Cluster",override.aes = list(size=5)),
                                     shape=guide_legend(title="Tissue+Isotype",override.aes = list(size=5)))+
                                            scale_color_manual(labels = c("noise", "1","2","3","4","5"), values=c(1:6))+ 
                                            geom_text_repel(label=paste(rownames(x.clean),as.factor(cld$cluster)))
                                            
                            png("umap_withdbscan_cluster2.png", width=1100, height=700)   
                               #ggarrange(h, g3,ncol=1, nrow=2) 
                                g4
                               dev.off()
                               png("umap_withhdbscan_cluster1.png", width=1100, height=700)   
                               #ggarrange(h, g3,ncol=1, nrow=2) 
                                plot(cl,show_flat=T)
                               dev.off()
                               
                               #======================================================
                               #################doing the cluster by tissue and Isotype
                               #+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                                    #==== dbscan======
                               #xdat<-as.data.frame(umap_object$layout)
                               xdat<-as.data.frame(x.clean)
                               xdat<-cbind(xdat, conds[rownames(xdat),c("isotype", "tissue","treatment")])
                               xdat.ti<-xdat[xdat$isotype=="IgG"& xdat$tissue=="Bone Marrow",]
                               cl<-dbscan(xdat.ti[,c(1:n_pcs)],minPts=3,eps=0.4)  # minPts the minimal # of points to form a cluster.
                               colnames(xdat.ti)[1:2]=c("X","Y")
                               
                               gt<-ggplot(aes(x = X, y = Y), data = xdat.ti) +
                                                geom_point(aes(shape=factor(paste0(tissue,"+",isotype)), color=factor(cl$cluster+1L)),size=10)+#, size=(x.complications.sum+1)*2)+
                                                ggtitle("dbscan clustering")+guides(col=guide_legend(title="Cluster",override.aes = list(size=5)),
                                     shape=guide_legend(title="Tissue+Isotype",override.aes = list(size=5)))+
                                            scale_color_manual(labels = c("noise", "1","2","3","4","5"), values=c(1:6))+ 
                                            geom_text_repel(label=paste(rownames(xdat.ti),as.factor(cl$cluster)))
                                            
           ###########======================
           #            END of HDBScan section ######
           #################################
                
                ######========spectral clustering.
    #get data first
  #  xdat<-as.data.frame(umap_object$layout)
xdat<-as.data.frame(x.clean)#d.out$ind$coord)
colnames(xdat)<-paste0("V",c(1:dim(x.clean)[2]))#c("V1","V2","V3","V4"))
   xdat<-cbind(xdat, conds[rownames(xdat),c("isotype", "tissue","treatment")])
   tissue<-c( "Bone Marrow","Spleen")
   isotype<-c("IgM","IgG" )
   h.ti<-vector("list", length=4)
      g3.ti<-vector("list", length=4)  #remembering by compartment clustering by first 2pcs
      g4.ti<-vector("list", length=4)# rember by compartment clusters by umap 2 d
      xdat.pc.ti<-vector("list", length=4)#remembering data to do plotting in another unit
      xdat.umap.ti<-vector("list", length=4)
      v_small<-vector("list", length=4)
      pv_small<-vector("list", length=4)
      index<-0
   for(i in tissue)
   {
        for(j in isotype)
        {
            index<-index+1
            xdat.ti<-xdat[as.character(xdat$isotype)==j& as.character(xdat$tissue)==i,]
            
            #x.dist<-daisy(as.data.frame(xdat.ti[,c(1:(dim(xdat.ti)[2]-3))]), metric=d_method, stand=F)
            x.dist<-dist(as.data.frame(xdat.ti[,c(1:(dim(xdat.ti)[2]-3-0))]), method=d_method, p=pmink)
 
            x.sim<-exp(-1*(as.matrix(x.dist)))
            A <- make.affinity(x.sim,2)  # use 3 neighbors (includes self)   
              
            eigen_m <- eigs(A,1,which="SR")
            toadd = ceiling(abs(min(eigen_m$values)))
            A.norm<-A
              diag(A.norm) = diag(A.norm) + toadd
              A.norm = normalizeKernel(A.norm)
           
            csize<-4
           #mc2_specc<-apclusterK(x.sim, details=T, K=csize, maxit=2000)
            #mc2_specc<-apcluster(x.sim, details=T)#, K=csize)
            #v_small.list<-rep(1,dim(x.sim)[1])
            #for(k in 1:csize){
             #   v_small.list[mc2_specc@clusters[[k]]]<-k
            #}
            #v_small[[index]]<-v_small.list
            pamc = pam(as.matrix(x.dist), diss =T, k =4)
            pv_small[[index]]<-pamc$clustering
            mc2_specc = specc(as.kernelMatrix(A.norm), centers = 4)
           v_small [[index]]= (mc2_specc@.Data)
           #v_small[[index]]<-pv_small[[index]]
            g = graph_from_adjacency_matrix(A.norm, mode = "undirected", weighted = TRUE, diag = FALSE)
            
            h.ti[[index]]<- plot.igraph(g, edge.width = E(g)$weight/(max(E(g)$weight)/5), vertex.label =rownames(x.sim), 
                    vertex.label.color="red",vertex.label.dist=1.5, vertex.label.font=4,
                    vertex.size = as.integer(xdat.ti$isotype)*5, 
                    edge.curved = T, vertex.shape=shapes[as.integer(factor(paste0(xdat.ti$tissue,"+",xdat.ti$isotype)))],
                    vertex.color=v_small[[index]])
   
        xdat.ti$cluster<-as.factor(v_small[[index]])        
              g3.ti[[index]]<-ggplot(aes(x =V1, y = V2), data =xdat.ti) +
  geom_point(aes(shape=treatment,color=cluster),size=10)+#, size=(x.complications.sum+1)*2)+
                ggtitle("spectral clustering")+scale_shape_manual(values=c(19,17,15,8,7))+
                geom_text_repel(label=paste0(rownames(xdat.ti),xdat.ti$cluster))
                
                temp<-tsne_df[tsne_df$isotype==j&tsne_df$tissue==i,]
                temp$treatment<-factor(temp$treatment, levels=c("PBS","OVA", "OVA+PorB","OVA+CpG","OVA+Alum"))
                temp$cluster<-as.factor(v_small[[index]])  
                
              g4.ti[[index]]<-ggplot(aes(x =X, y = Y), data =temp) +
                geom_point(aes(shape=treatment,color=cluster),size=5)+#, size=(x.complications.sum+1)*2)+
                ggtitle(paste0(i,"+",j))+scale_shape_manual(values=c(19,17,15,8,7))+
                theme_bw(base_size=13)#+theme(legend.position="none")
                #geom_text_repel(label=paste0(rownames(temp),temp$cluster))    
                xdat.pc.ti[[index]]<-xdat.ti
                xdat.umap.ti[[index]]<-temp
        }
   }
   
    png(file=here(output.dir,"cluster_byIsotypeTissue_PCA1_2.png"), 
        width=1200, height=1000)
    ggarrange(g3.ti[[1]],g3.ti[[2]],g3.ti[[3]],g3.ti[[4]], nrow=2, ncol=2, common.legend = TRUE)
    dev.off()
    
    png(file=here(output.dir,"cluster_byIsotypeTissue_UMap1_2.png"), 
        width=1200, height=1000)
     ggarrange(g4.ti[[1]],g4.ti[[2]],g4.ti[[3]],g4.ti[[4]], nrow=2, ncol=2, common.legend = TRUE, legend="bottom")
     dev.off()
        
        #save the data and 
        #setwd("/home/feng/Windows/windowsD/feng/LAB/hg/IgSeq_MS/manuscript/clustering/pc3")
        
        save(g3.ti, g4.ti, xdat.pc.ti, 
            xdat.umap.ti, isotype, tissue, 
            file=here(data.dir,"figure_cluster_data_by_compartment.RData")
            )




