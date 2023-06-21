#Relict dace - contract report analyses and figures

#Load packages
library(ggplot2)

#Import sample metadata
rd <- read.csv("https://raw.githubusercontent.com/meflacava/RelictDace/main/metadata/AllProjects_RD_Metadata.csv")


#### Range-wide ####

## Sample metadata
#Import bamlist
rd_meta <- read.table("https://raw.githubusercontent.com/meflacava/RelictDace/main/bamlists/bam.rwpopgen.txt",header=F)

#add metadata to rd_meta
names(rd_meta) <- "SampleID"
rd_meta$SampleID <- sub("_sorted_proper_rmdup.bam.*","",rd_meta$SampleID)
rd_meta$SampleID <- sub(".*/home/mlacava/align/RAD_alignments/","",rd_meta$SampleID)
for (i in 1:nrow(rd_meta)){
  rd_meta$Valley[i] <- rd$Valley[rd$SampleID==rd_meta$SampleID[i]]
  rd_meta$Site[i] <- rd$Site[rd$SampleID==rd_meta$SampleID[i]]
  rd_meta$SiteWaypoint[i] <- rd$SiteWaypoint[rd$SampleID==rd_meta$SampleID[i]]
}
rd_meta$Valley <- factor(rd_meta$Valley)
rd_meta$Site <- factor(rd_meta$Site,
                       levels=c("RBFL","RBFR","RBRL","BTOC","BTQP","BTST","GSTS","GSJS",
                                "STCO","STPR","STMR","STRT","STWA","SPSH","SPKR"))
rd_meta$SiteWaypoint <- factor(rd_meta$SiteWaypoint,
                               levels=c("RBFL","RBFR","RBRL_P226","RBRL_P246","BTOC","BTQP","BTST",
                                        "BTST_NSR","GSTS","GSJS_BS","GSJS_BSN","GSJS_CS02","GSJS_CS03",
                                        "GSJS_CS04","GSJS_CS11","GSJS_CS16","GSJS_CS22A","GSJS_D01",
                                        "GSJS_D09","GSJS_NS03","GSJS_NS05","GSJS_SS19","STCO",
                                        "STPR_BT2","STPR_FT1","STMR","STRT","STWA","SPSH","SPKR"))



###### PC1 x PC2 (color by site) ####

#Import covariance matrix
rd_mat <- as.matrix(read.table("https://raw.githubusercontent.com/meflacava/RelictDace/main/results/rwpopgen_PCA.covMat"))

#prep matrix for plotting
eig <- eigen(rd_mat, symm=T)
eig$val <- eig$val/sum(eig$val)
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))

#adding metadata parameters
PC$Valley <- rd_meta$Valley
PC$Site <- rd_meta$Site
PC$SiteWaypoint <- rd_meta$SiteWaypoint

#plot PCA
pdf("rw_pca1x2.pdf",
    width=5,height=4.5) #default=7
ggplot2::ggplot(PC, ggplot2::aes(x = PC1,y = PC2, color = Site)) + 
  ggplot2::geom_point(cex = 4, alpha = 0.75) + 
  ggplot2::theme_classic() + 
  ggplot2::xlab(paste0("PC1, ", round(eig$val[1] * 100, 2), "% variance explained")) + 
  ggplot2::ylab(paste0("PC2, ", round(eig$val[2] * 100, 2), "% variance explained")) + 
  ggplot2::guides(color = ggplot2::guide_legend(title = "Site"))
dev.off()



##### Admix K=2 ####

i<-2 #set k value you are plotting
pdf(paste0("rw_k",i,".pdf"),
    width=9,height=3) #default=7
q <- as.matrix(read.table(paste0("https://raw.githubusercontent.com/meflacava/RelictDace/main/results/admix/rw_allsamples_2017snps_k",i,"_run1.qopt")))
group <- rd_meta$SiteWaypoint
grp_name <- levels(factor(group))
idx <- list() #separate samples into groups
for (n in grp_name) {
  x <- which(group == n)
  x <- x[order(q[x, 1L], decreasing = TRUE)] #order q matrix by decreasing ancestry proportions
  idx <- c(idx, list(x))
}
grp_len <- lengths(idx, use.names = FALSE)
xl <- c(0, cumsum(grp_len))
x <- xl[-1L] - 0.5 * grp_len
q <- q[unlist(idx),] #not sure what this does, but without it, site dividing lines don't work
#barplot(t(q),col=hcl.colors(i,palette="Spectral"),space=0,border=NA,yaxt="n",
#        xlab="",ylab="Ancestry proportion") #as many colors in spectrum as there are clusters
barplot(t(q),col=hcl.colors(3,palette="Spectral")[c(1,3)],space=0,border=NA,yaxt="n", #spectral palette messed up for n=2, so subset from n=3
        xlab="",ylab="Ancestry proportion") #as many colors in spectrum as there are clusters
abline(v=xl,col="black")
axis(1,x,grp_name,las=3,cex.axis=0.75,tick=F,hadj=0.8)
axis(2,at=seq(0,1,0.25),las=1,cex.axis=0.74) 
dev.off()





#### Western valleys ####

## Sample metadata
#Import bamlist
rd_meta <- read.table("https://raw.githubusercontent.com/meflacava/RelictDace/main/bamlists/bam.west.txt",header=F)

#add metadata to rd_meta
names(rd_meta) <- "SampleID"
rd_meta$SampleID <- sub("_sorted_proper_rmdup.bam.*","",rd_meta$SampleID)
rd_meta$SampleID <- sub(".*/home/mlacava/align/RAD_alignments/","",rd_meta$SampleID)
for (i in 1:nrow(rd_meta)){
  rd_meta$Valley[i] <- rd$Valley[rd$SampleID==rd_meta$SampleID[i]]
  rd_meta$Site[i] <- rd$Site[rd$SampleID==rd_meta$SampleID[i]]
  rd_meta$SiteWaypoint[i] <- rd$SiteWaypoint[rd$SampleID==rd_meta$SampleID[i]]
}
rd_meta$Valley <- factor(rd_meta$Valley)
rd_meta$Site <- factor(rd_meta$Site,
                       levels=c("RBFL","RBFR","RBRL","BTOC","BTQP","BTST","GSTS","GSJS",
                                "STCO","STPR","STMR","STRT","STWA","SPSH","SPKR"))
rd_meta$SiteWaypoint <- factor(rd_meta$SiteWaypoint,
                               levels=c("RBFL","RBFR","RBRL_P226","RBRL_P246","BTOC","BTQP","BTST",
                                        "BTST_NSR","GSTS","GSJS_BS","GSJS_BSN","GSJS_CS02","GSJS_CS03",
                                        "GSJS_CS04","GSJS_CS11","GSJS_CS16","GSJS_CS22A","GSJS_D01",
                                        "GSJS_D09","GSJS_NS03","GSJS_NS05","GSJS_SS19","STCO",
                                        "STPR_BT2","STPR_FT1","STMR","STRT","STWA","SPSH","SPKR"))

##### PC1 x PC2 ####

#Import covariance matrix
rd_mat <- as.matrix(read.table("https://raw.githubusercontent.com/meflacava/RelictDace/main/results/west_PCA.covMat"))

#prep matrix for plotting
eig <- eigen(rd_mat, symm=T)
eig$val <- eig$val/sum(eig$val)
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))

#adding metadata parameters
PC$Valley <- rd_meta$Valley
PC$Site <- rd_meta$Site
PC$SiteWaypoint <- rd_meta$SiteWaypoint

#plot PC1 x PC 2
pdf("west_pca1x2.pdf",
    width=5.5,height=4.5) #default=7
ggplot2::ggplot(PC, ggplot2::aes(x = PC1,y = PC2, color = SiteWaypoint)) + 
  ggplot2::geom_point(cex = 4, alpha = 0.75) + 
  ggplot2::theme_classic() + 
  ggplot2::xlab(paste0("PC1, ", round(eig$val[1] * 100, 2), "% variance explained")) + 
  ggplot2::ylab(paste0("PC2, ", round(eig$val[2] * 100, 2), "% variance explained")) + 
  ggplot2::guides(color = ggplot2::guide_legend(title = "Site"))
dev.off()


##### PC2 x PC3 ####
pdf("west_pca2x3.pdf",
    width=5.5,height=4.5) #default=7
ggplot2::ggplot(PC, ggplot2::aes(x = PC2,y = PC3, color = SiteWaypoint)) + 
  ggplot2::geom_point(cex = 4, alpha = 0.75) + 
  ggplot2::theme_classic() + 
  ggplot2::xlab(paste0("PC2, ", round(eig$val[2] * 100, 2), "% variance explained")) + 
  ggplot2::ylab(paste0("PC3, ", round(eig$val[3] * 100, 2), "% variance explained")) + 
  ggplot2::guides(color = ggplot2::guide_legend(title = "Site"))
dev.off()




##### Admix K=2 ####

i<-2 #set k value you are plotting
pdf(paste0("west_k",i,".pdf"),
    width=9,height=3) #default=7
q <- as.matrix(read.table(paste0("https://raw.githubusercontent.com/meflacava/RelictDace/main/results/admix/west_allsamples_2017snps_k",i,"_run1.qopt")))
group <- rd_meta$SiteWaypoint
grp_name <- levels(factor(group))
idx <- list() #separate samples into groups
for (n in grp_name) {
  x <- which(group == n)
  x <- x[order(q[x, 1L], decreasing = TRUE)] #order q matrix by decreasing ancestry proportions
  idx <- c(idx, list(x))
}
grp_len <- lengths(idx, use.names = FALSE)
xl <- c(0, cumsum(grp_len))
x <- xl[-1L] - 0.5 * grp_len
q <- q[unlist(idx),] #not sure what this does, but without it, site dividing lines don't work
#barplot(t(q),col=hcl.colors(i,palette="Spectral"),space=0,border=NA,yaxt="n",
#        xlab="",ylab="Ancestry proportion") #as many colors in spectrum as there are clusters
barplot(t(q),col=hcl.colors(3,palette="Spectral")[c(1,3)],space=0,border=NA,yaxt="n", #spectral palette messed up for n=2, so subset from n=3
        xlab="",ylab="Ancestry proportion") #as many colors in spectrum as there are clusters
abline(v=xl,col="black")
axis(1,x,grp_name,las=3,cex.axis=0.75,tick=F,hadj=0.8)
axis(2,at=seq(0,1,0.25),las=1,cex.axis=0.74) 
dev.off()




##### Admix K=6 ####

i<-6 #set k value you are plotting
pdf(paste0("west_k",i,".pdf"),
    width=9,height=3) #default=7
q <- as.matrix(read.table(paste0("https://raw.githubusercontent.com/meflacava/RelictDace/main/results/admix/west_allsamples_2017snps_k",i,"_run1.qopt")))
group <- rd_meta$SiteWaypoint
grp_name <- levels(factor(group))
idx <- list() #separate samples into groups
for (n in grp_name) {
  x <- which(group == n)
  x <- x[order(q[x, 1L], decreasing = TRUE)] #order q matrix by decreasing ancestry proportions
  idx <- c(idx, list(x))
}
grp_len <- lengths(idx, use.names = FALSE)
xl <- c(0, cumsum(grp_len))
x <- xl[-1L] - 0.5 * grp_len
q <- q[unlist(idx),] #not sure what this does, but without it, site dividing lines don't work
barplot(t(q),col=hcl.colors(i,palette="Spectral"),space=0,border=NA,yaxt="n",
        #main=paste0("k = ",i),
        xlab="",ylab="Ancestry proportion") #as many colors in spectrum as there are clusters
abline(v=xl,col="black")
axis(1,x,grp_name,las=3,cex.axis=0.75,tick=F,hadj=0.8)
axis(2,at=seq(0,1,0.25),las=1,cex.axis=0.74) 
dev.off()




##### Diversity (theta pi) ####

#List of sites
pops <- read.table("https://raw.githubusercontent.com/meflacava/RelictDace/main/bamlists/west.pop.list")
pops <- pops$V1

#For each site, import diversity stats and calculate standardized theta pi and watterson's theta
theta <- data.frame(site=NA,tPstd=NA,tWstd=NA)
for (i in pops){
  x <- read.table(paste0("https://raw.githubusercontent.com/meflacava/RelictDace/main/results/diversity/",i,".thetas.idx.pestPG"))
  names(x) <- c("index","chr","wincenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites") #got from ANGSD website
  x$tPstd <- x$tP/x$nSites #standardize by number of sites
  x$tWstd <- x$tW/x$nSites #standardize by number of sites
  theta <- rbind(theta,data.frame(site=rep(i,length(x$tPstd)),tPstd=x$tPstd,tWstd=x$tWstd))
}
theta <- theta[-1,]
  
#boxplot
pdf(paste0("west_div.pdf"),
    width=5,height=5) #default=7
boxplot(tPstd~site,data=theta,xlab="",ylab="Theta Pi",las=1,cex.axis=0.8,outline=F,xaxt="n")
axis(side=1,at=1:length(unique(theta$site)),labels=unique(theta$site),las=2,tick=F,cex.axis=0.8,line=-0.8)
dev.off()

#mean theta pi across all east sites
mean(theta$tPstd) #0.006613939
sd(theta$tPstd) #0.003421593
median(theta$tPstd) #0.006048985

#highest/lowest mean site
tP.ave <- tapply(theta$tPstd,theta$site,mean)
sort(tP.ave)
mean(theta$tPstd[theta$site=="BTST_NSR"]) #0.007146443
mean(theta$tPstd[theta$site=="RBRL_P226"]) #0.005862045


##### Fst ####

#import global Fsts
fst <- read.csv("https://raw.githubusercontent.com/meflacava/RelictDace/main/results/compiledFsts.csv")
head(fst)
names(fst) <- c("pop1","pop2","UnweightFst","Fst") #change weight Fst name for legend plotting

#select subset of site pairs to plot
west <- read.table("https://raw.githubusercontent.com/meflacava/RelictDace/main/bamlists/west.pop.list")
west <- west$V1
wfst <- fst[fst$pop1 %in% west,]

#west tile plot
pdf(paste0("west_fst.pdf"),
    width=5.5,height=5) #default=7
ggplot(wfst, aes(pop2, pop1)) + 
  geom_tile(aes(fill = Fst), color="white") + 
  theme_minimal() + 
  scale_fill_gradient(low="#98afd4",high="#04255c") +
  #scale_fill_viridis_c(direction=-1,option = "cividis") + #Finger 2022 manuscript colors
  theme_bw() + 
  theme(axis.text.x=element_text(angle=45,vjust=1,size=10,hjust=1), panel.border=element_blank(),
        panel.grid.major=element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank()) + scale_y_discrete(position = "left") + 
  geom_text(aes(label = round(Fst, 2)), size=3.5, colour = "white") + 
  #ggtitle("D") + 
  theme(plot.title = element_text(face="bold", size=16))
dev.off()




#### Eastern valleys ####

## Sample metadata
#Import bamlist
rd_meta <- read.table("https://raw.githubusercontent.com/meflacava/RelictDace/main/bamlists/bam.east.txt",header=F)

#add metadata to rd_meta
names(rd_meta) <- "SampleID"
rd_meta$SampleID <- sub("_sorted_proper_rmdup.bam.*","",rd_meta$SampleID)
rd_meta$SampleID <- sub(".*/home/mlacava/align/RAD_alignments/","",rd_meta$SampleID)
for (i in 1:nrow(rd_meta)){
  rd_meta$Valley[i] <- rd$Valley[rd$SampleID==rd_meta$SampleID[i]]
  rd_meta$Site[i] <- rd$Site[rd$SampleID==rd_meta$SampleID[i]]
  rd_meta$SiteWaypoint[i] <- rd$SiteWaypoint[rd$SampleID==rd_meta$SampleID[i]]
}
rd_meta$Valley <- factor(rd_meta$Valley)
rd_meta$Site <- factor(rd_meta$Site,
                       levels=c("RBFL","RBFR","RBRL","BTOC","BTQP","BTST","GSTS","GSJS",
                                "STCO","STPR","STMR","STRT","STWA","SPSH","SPKR"))
rd_meta$SiteWaypoint <- factor(rd_meta$SiteWaypoint,
                               levels=c("RBFL","RBFR","RBRL_P226","RBRL_P246","BTOC","BTQP","BTST",
                                        "BTST_NSR","GSTS","GSJS_BS","GSJS_BSN","GSJS_CS02","GSJS_CS03",
                                        "GSJS_CS04","GSJS_CS11","GSJS_CS16","GSJS_CS22A","GSJS_D01",
                                        "GSJS_D09","GSJS_NS03","GSJS_NS05","GSJS_SS19","STCO",
                                        "STPR_BT2","STPR_FT1","STMR","STRT","STWA","SPSH","SPKR"))


##### PC1 x PC2 ####

#Import covariance matrix
rd_mat <- as.matrix(read.table("https://raw.githubusercontent.com/meflacava/RelictDace/main/results/east_PCA.covMat"))

#prep matrix for plotting
eig <- eigen(rd_mat, symm=T)
eig$val <- eig$val/sum(eig$val)
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))

#adding metadata parameters
PC$Valley <- rd_meta$Valley
PC$Site <- rd_meta$Site
PC$SiteWaypoint <- rd_meta$SiteWaypoint

#plot PC1 x PC 2
pdf("east_pca1x2.pdf",
    width=6.2,height=4) #default=7
ggplot2::ggplot(PC, ggplot2::aes(x = PC1,y = PC2, color = SiteWaypoint)) + 
  ggplot2::geom_point(cex = 4, alpha = 0.75) + 
  ggplot2::theme_classic() + 
  ggplot2::xlab(paste0("PC1, ", round(eig$val[1] * 100, 2), "% variance explained")) + 
  ggplot2::ylab(paste0("PC2, ", round(eig$val[2] * 100, 2), "% variance explained")) + 
  ggplot2::guides(color = ggplot2::guide_legend(title = "Site"))
dev.off()


##### PC3 x PC4 ####
pdf("east_pca3x4.pdf",
    width=6.2,height=4) #default=7
ggplot2::ggplot(PC, ggplot2::aes(x = PC3,y = PC4, color = SiteWaypoint)) + 
  ggplot2::geom_point(cex = 4, alpha = 0.75) + 
  ggplot2::theme_classic() + 
  ggplot2::xlab(paste0("PC3, ", round(eig$val[3] * 100, 2), "% variance explained")) + 
  ggplot2::ylab(paste0("PC4, ", round(eig$val[4] * 100, 2), "% variance explained")) + 
  ggplot2::guides(color = ggplot2::guide_legend(title = "Site"))
dev.off()



##### Admix K=16 ####

i<-16 #set k value you are plotting
pdf(paste0("east_k",i,".pdf"),
    width=9,height=3) #default=7
q <- as.matrix(read.table(paste0("https://raw.githubusercontent.com/meflacava/RelictDace/main/results/admix/east_allsamples_2017snps_k",i,"_run1.qopt")))
group <- rd_meta$SiteWaypoint
grp_name <- levels(factor(group))
idx <- list() #separate samples into groups
for (n in grp_name) {
  x <- which(group == n)
  x <- x[order(q[x, 1L], decreasing = TRUE)] #order q matrix by decreasing ancestry proportions
  idx <- c(idx, list(x))
}
grp_len <- lengths(idx, use.names = FALSE)
xl <- c(0, cumsum(grp_len))
x <- xl[-1L] - 0.5 * grp_len
q <- q[unlist(idx),] #not sure what this does, but without it, site dividing lines don't work
barplot(t(q),col=hcl.colors(i,palette="Spectral"),space=0,border=NA,yaxt="n",
        #main=paste0("k = ",i),
        xlab="",ylab="Ancestry proportion") #as many colors in spectrum as there are clusters
abline(v=xl,col="black")
axis(1,x,grp_name,las=3,cex.axis=0.75,tick=F,hadj=0.8)
axis(2,at=seq(0,1,0.25),las=1,cex.axis=0.74) 
dev.off()






##### Diversity (theta pi) ####

#List of sites
pops <- read.table("https://raw.githubusercontent.com/meflacava/RelictDace/main/bamlists/east.pop.list")
pops <- pops$V1

#For each site, import diversity stats and calculate standardized theta pi and watterson's theta
theta <- data.frame(site=NA,tPstd=NA,tWstd=NA)
for (i in pops){
  x <- read.table(paste0("https://raw.githubusercontent.com/meflacava/RelictDace/main/results/diversity/",i,".thetas.idx.pestPG"))
  names(x) <- c("index","chr","wincenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites") #got from ANGSD website
  x$tPstd <- x$tP/x$nSites #standardize by number of sites
  x$tWstd <- x$tW/x$nSites #standardize by number of sites
  theta <- rbind(theta,data.frame(site=rep(i,length(x$tPstd)),tPstd=x$tPstd,tWstd=x$tWstd))
}
theta <- theta[-1,]

#boxplot
pdf(paste0("east_div.pdf"),
    width=10,height=5) #default=7
boxplot(tPstd~site,data=theta,xlab="",ylab="Theta Pi",las=1,cex.axis=0.8,outline=F,xaxt="n")
axis(side=1,at=1:length(unique(theta$site)),labels=unique(theta$site),las=2,tick=F,cex.axis=0.8,line=-0.8)
dev.off()

#mean theta pi across all east sites
mean(theta$tPstd) #0.007840398
sd(theta$tPstd) #0.003763411
median(theta$tPstd) #0.007261283

#highest/lowest mean site
tP.ave <- tapply(theta$tPstd,theta$site,mean)
sort(tP.ave)
mean(theta$tPstd[theta$site=="STWA"]) #0.008513493
sd(theta$tPstd[theta$site=="STWA"]) #0.003535682
mean(theta$tPstd[theta$site=="GSTS"]) #0.006590664
sd(theta$tPstd[theta$site=="GSTS"]) #0.004479207


##### Fst ####

#import global Fsts
fst <- read.csv("https://raw.githubusercontent.com/meflacava/RelictDace/main/results/compiledFsts.csv")
head(fst)
names(fst) <- c("pop1","pop2","UnweightFst","Fst") #change weight Fst name for legend plotting

#select subset of site pairs to plot
east <- read.table("https://raw.githubusercontent.com/meflacava/RelictDace/main/bamlists/east.pop.list")
east <- east$V1
efst <- fst[fst$pop1 %in% east,]

#reorder so GSTS is on the bottom of the y axis
efst$pop1 <- factor(efst$pop1,levels=c("GSTS","GSJS_BS","GSJS_BSN","GSJS_CS02","GSJS_CS03","GSJS_CS04","GSJS_CS11",
                                       "GSJS_CS16","GSJS_CS22A","GSJS_D01","GSJS_D09","GSJS_NS03","GSJS_NS05",
                                       "GSJS_SS19","SPKR","SPSH","STCO","STMR","STPR_BT2","STPR_FT1","STRT"))

#east tile plot
pdf(paste0("east_fst.pdf"),
    width=10,height=5) #default=7
ggplot(efst, aes(pop2, pop1)) + 
  geom_tile(aes(fill = Fst), color="white") + 
  theme_minimal() + 
  scale_fill_gradient(low="#98afd4",high="#04255c") +
  #scale_fill_viridis_c(direction=-1,option = "cividis") + #Finger 2022 manuscript colors
  theme_bw() + 
  theme(axis.text.x=element_text(angle=45,vjust=1,size=10,hjust=1), panel.border=element_blank(),
        panel.grid.major=element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank()) + scale_y_discrete(position = "left") + 
  geom_text(aes(label = round(Fst, 2)), size=3.5, colour = "white") + 
  #ggtitle("D") + 
  theme(plot.title = element_text(face="bold", size=16))
dev.off()

