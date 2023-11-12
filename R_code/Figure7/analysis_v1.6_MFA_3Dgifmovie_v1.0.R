#analysis to draw 3D gif movie to 
# first 3 PCs of data arrangement

# # Make sure you have read ./ReadMe.txt for instructions
# of running scripts and preparing data.


 library(ggpubr)
 library(ggrepel)
# library(ggalt)
 library(here)
 library(ggplot2)

 library(rgl)
 library(knitr)
 library(rglwidget)
 library(magick)
 library(webshot2)
 
 #set up to let rgl to use the web browser to
 # show output.
 options(rgl.printRglwidget = TRUE)
 options(browser = "firefox") #<- make sure to set this correctly
 

 data.dir<-"Data/Figure7"
 output.dir<-"R_code/Figure7"

 #read data.
load(here(data.dir,"data_for_3d_gifmovie.RData")) #read in  data containing everthing.  d.all, 
   
   #start doing it. 
   device_id<-open3d()# start a new fresh drawing
   plot3d( xlab="PC1", ylab="PC2",zlab="PC3",
        x=d.quali[,1],y=d.quali[,2], z=d.quali[,3],
                    col=1,#as.integer(as.factor(d.quali$factor))+1, 
                    pch=-1, #c(1:8),#pch=as.integer(as.factor(rownames(d.quali))),
                    size=3
                    )
            
    #play3d(spin3d( rpm=4), duration=10)
    #set up pch for plotting.
    d.quali$pch<-15 #square 
    d.quali[2,"pch"]<-16 #circle 
    d.quali[8,"pch"]<-16# circle 
    d.quali[4,"pch"]<-16 #circle 
    d.quali[5,"pch"]<-17#triangle
    d.quali[6,"pch"]<-18#diamond 
    #legend, OVA +Alum: circle; OVA : square;  diamond : PBS ; triangle: OVA+CpG
        ##pch3d keeps calling on a new interface (
        # browser page). we can not stop it. 
        # so the effects are that we open a new face but 
        # redrawing the pchs/points to the preexisting
        #   plots
        #
        #but we can call explicitly to open a new 
        # interface (without see the previous plottings)
        pch3d( xlab="PC1", ylab="PC2",zlab="PC3",
                x=d.quali[,1],y=d.quali[,2], z=d.quali[,3],
                    col=as.integer(as.factor(d.quali$factor))+3, 
                    pch=as.integer(d.quali$pch),#pch=as.integer(as.factor(rownames(d.quali))),
                    cex=1
                    )
    # method 1:
    #the following line to set up how to rotate
    # and duration
    play3d(spin3d( axis=c(0,0,1),rpm=4), duration=10)
    
    #this next line to generate the movie(pngs) and record
    # by calling underlying facility to generate gif
    # the movie is in "tempdir()" by default
    movie3d(spin3d( rpm=3), duration=5,
        dir=tempdir())
    
    #
    cat("the generated movie is ",
        paste0(tempdir(),"/movie.gif"),"\n"
        ) 
    #just need the belew image as a legend to the movie   
   png(file=here(output.dir,"MFA_movie3d_quali_legend.png"), width=600, height=600)
    plot(x=c(0,1),y=c(0,1), type="n", yaxt="n", xaxt="n", xlab="", ylab="")
    legend (x=0.75, y=0.4, legend=c("IgG","IgM","PBS", "OVA","OVA+CpG", "OVA+Alum", "Spleen","Bone Marrow")  
                                            ,col=c(2,2,4,4,4,4,3,3), pch=c(15,16, 18,15,17,16, 16,15), cex=1.3)
    dev.off()
    #save(file=here(data.dir,"3d_dd.RData"), dd)
    
    #------------------------------
    ###here we do another movie for rotating the points in 3d (NOT supplimentary variables)
    # open3d()  <- open a new window.
    #pch3d(x1, y1, z1, color="red", bg="red", pch=point.styles)
    #pch3d(x2, y2, z2, color="blue", bg="blue", pch=point.styles)
    #pch3d(x3, y3, z3, color="green", bg="green", pch=point.styles)

    #To do movie we need plot3d and pch3d (have different shape of chars)
    #drawing all individual points with different chars 
    #plot 3d to set scale and limits
    device_id<-open3d()  #<- open a new window.
     plot3d( xlab="PC1", ylab="PC2",zlab="PC3",
        x=dd[,1],y=dd[,2], z=dd[,3],
                        col=1, 
                        pch=-1,#as.integer(dd.iso$treatment),#pch=as.integer(as.factor(rownames(d.quali))),
                        size=0
                        )
                        

    dd.iso<-dd#[dd$isotype=="IgG"&dd$tissue=="Spleen",]
    dd.iso$pch<-18
    dd.iso[dd.iso$treatment=="OVA","pch"]<-15
    dd.iso[dd.iso$treatment=="OVA+CpG","pch"]<-16
    dd.iso[dd.iso$treatment=="OVA+Alum","pch"]<-17

     pch3d( xlab="PC1", ylab="PC2",zlab="PC3",
            x=dd.iso[,1],y=dd.iso[,2], z=dd.iso[,3],
            col=as.integer(as.factor(paste0(dd.iso$isotype,"+",dd.iso$tissue))), 
            pch=as.integer(dd.iso$pch),#pch=as.integer(as.factor(rownames(d.quali))),
            cex=0.5
        )
              
    #axes3d()  <-not necessary since we calling plot3d to set up the canvas

    play3d(spin3d( rpm=4), duration=10)
    
    movie3d(spin3d( rpm=3), duration=5,
        movie="movieTissueIso", dir=tempdir())
    cat("the generated movie is ",
        paste0(tempdir(),"/movieTissueIso.gif"),"\n"
        )
    #draw legend as a separate image
    png(file=here(output.dir,"MFA_movie3d_ind_legend.png")
        , width=600, height=600)
        plot(x=c(0,1),y=c(0,1), type="n", yaxt="n", xaxt="n", xlab="", ylab="")
        legend (x=0.65, y=0.8, legend=c("IgG+Spleen", "IgM+Spleen", "IgG+Bone Marrow", "IgM+Bone Marrow")  
                                                    ,col=c(2,4,1,3), pch=15, cex=1)
        legend (x=0.65, y=0.4, legend=c("PBS", "OVA", "OVA+CpG", "OVA+Alum")  
                                                    ,col=1, pch=c(18,15, 16,17), cex=1)
     dev.off()


    ############
    # method 2: more recent code and have a better control
    #    1) set up the size of the dimension of output
    #    2) control to loop indefinitely the gif
    ##############

    #first do by tissue
    device_id<-open3d() 
    plot3d( xlab="PC1", ylab="PC2",zlab="PC3",
        x=d.quali[,1],y=d.quali[,2], z=d.quali[,3],
                    col=1,#as.integer(as.factor(d.quali$factor))+1, 
                    pch=-1, #c(1:8),#pch=as.integer(as.factor(rownames(d.quali))),
                    size=1
                    )
            
    #play3d(spin3d( rpm=4), duration=10)
    #set up pch for plotting.
    d.quali$pch<-15 #square 
    d.quali[2,"pch"]<-16 #circle 
    d.quali[8,"pch"]<-16# circle 
    d.quali[4,"pch"]<-16 #circle 
    d.quali[5,"pch"]<-17#triangle
    d.quali[6,"pch"]<-18#diamond 
    #legend, OVA +Alum: circle; OVA : square;  diamond : PBS ; triangle: OVA+CpG
    ##pch3d keeps calling on a new interface (
    # browser page). we can not stop it. 
    # so the effects are that we open a new face but 
    # redrawing the pchs/points to the preexisting
    #   plots
    #
    #but we can call explicitly to open a new 
    # interface (without see the previous plottings)
    pch3d( xlab="PC1", ylab="PC2",zlab="PC3",
            x=d.quali[,1],y=d.quali[,2], z=d.quali[,3],
                col=as.integer(as.factor(d.quali$factor))+3, 
                pch=as.integer(d.quali$pch),#pch=as.integer(as.factor(rownames(d.quali))),
                cex=0.5
                )
    rgl.bringtotop()
    
    olddir <- setwd(tempdir())
    angle.rotate=5                  
    for (i in 1:(180/angle.rotate)) {
        cat("i......",i, " and current device :",cur3d(),"\n")
        
        set3d(device_id)
       view3d(i*angle.rotate, 10)  #this is set the view point, here we do rotating around vertical axis. check the manual                              
       filename <- paste("pic", formatC(i, digits = 1, flag = "0"), ".png", sep = "")             
       snapshot3d(filename, width=750, height=550)                      
     }                  
     cat("\n")

     system("convert -delay 30 pic*.png -loop 0 pic.gif")
     system("rm -f pic*.png")
     setwd(olddir)                           
     cat("the generated movie is ",
        paste0(tempdir(),"/pic.gif"),"\n"
        )
    # draw a different movie by isotype and tissue
    device_id<-open3d()  #<- open a new window.
     plot3d( xlab="PC1", ylab="PC2",zlab="PC3",
        x=dd[,1],y=dd[,2], z=dd[,3],
                        col=1, 
                        pch=-1,#as.integer(dd.iso$treatment),#pch=as.integer(as.factor(rownames(d.quali))),
                        size=0
                        )
                        

    dd.iso<-dd#[dd$isotype=="IgG"&dd$tissue=="Spleen",]
    dd.iso$pch<-18
    dd.iso[dd.iso$treatment=="OVA","pch"]<-15
    dd.iso[dd.iso$treatment=="OVA+CpG","pch"]<-16
    dd.iso[dd.iso$treatment=="OVA+Alum","pch"]<-17

     pch3d( xlab="PC1", ylab="PC2",zlab="PC3",
            x=dd.iso[,1],y=dd.iso[,2], z=dd.iso[,3],
            col=as.integer(as.factor(paste0(dd.iso$isotype,"+",dd.iso$tissue))), 
            pch=as.integer(dd.iso$pch),#pch=as.integer(as.factor(rownames(d.quali))),
            cex=0.5
        )

    rgl.bringtotop()
    
    olddir <- setwd(tempdir())
    angle.rotate=5                  
    for (i in 1:(180/angle.rotate)) {
        cat("i......",i, " and current device :",cur3d(),"\n")
        
        set3d(device_id)
       view3d(i*angle.rotate, 0)  #this is set the view point, here we do rotating around vertical axis. check the manual                              
       filename <- paste("pic", formatC(i, digits = 1, flag = "0"), ".png", sep = "")             
       snapshot3d(filename, width=750, height=550)                      
     }                  
     cat("\n")

     system("convert -delay 30 pic*.png -loop 0 pic_tissueIso.gif")
     system("rm -f pic*.png")
     setwd(olddir) 

     cat("the generated movie is ",
        paste0(tempdir(),"/pic_tissueIso.gif"),"\n"
        )