rm(list = ls())
library(corrplot)
library(mashr)
setwd("./GTEX/08.analysis/MASH/0.rawdata/tissue/")
tissue_list <- list.files()
for (i in 1:48) {
  data <- read.delim(paste0("./GTEX/08.analysis/MASH/0.rawdata/tissue/",tissue_list[i]),header = F)
  data <- data[,c(1,2,3,5)]
  colnames(data) <- c("SNP","gene","beta","FDR")
  data$se <- sqrt(((data$beta)^2)/qchisq(data$FDR,1,lower.tail=F))
  data <- data[!is.na(data$se),]
  data <- data[data$SNP != "---",]
  data$tss <- "tss"
  data$sample <- "sample"
  data$count <- "count"
  data$maf <- "maf"
  data <- data[,c(2,1,6,7,8,9,4,3,5)]
  data$gene <- gsub(data$gene,pattern = "prmtr\\.",replacement = "")
  data$gene <- gsub(data$gene,pattern = "\\.ENSG.*\\.",replacement = "\\.")
  data$gene <- paste0("ENSG00",data$gene)
  colnames(data) <- c("gene_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se")
  write.table(data,paste0("./GTEX/08.analysis/MASH/0.rawdata/fin_tissue/",tissue_list[i]),col.names = T,row.names = F,sep = "\t",quote = F)
}

#Linux
#transform summary data to mash format
#!/bin/bash
#sos run workflows/fastqtl_to_mash --data-list ./tissue.txt --gene-list ./gene.txt -j 20
#echo "finish"

#/bin/bash
#sos run workflows/mashr_flashr_workflow.ipynb mash -j 30 --data ./fastqtl_to_mash_output/tissue.mash.rds

#magnitude
rm(list = ls())
out      <- readRDS("./gtexresults-master/fastqtl_to_mash_output/tissue.mash.rds")
maxb     <- out$strong.b
maxz     <- out$strong.z
out      <-readRDS("./gtexresults-master/mashr_flashr_workflow_output/tissue.mash.EZ.posterior.rds")
pm.mash        <- out$PosteriorMean
lfsr.all       <- out$lfsr
standard.error <- maxb/maxz
pm.mash.beta   <- pm.mash*standard.error
thresh       <- 0.05
pm.mash.beta <- pm.mash.beta[rowSums(lfsr.all<0.05)>0,]
lfsr.mash    <- lfsr.all[rowSums(lfsr.all<0.05)>0,]
shared.fold.size <- matrix(NA,nrow = ncol(lfsr.mash),ncol=ncol(lfsr.mash))
colnames(shared.fold.size) <- rownames(shared.fold.size) <- colnames(maxz)
for (i in 1:ncol(lfsr.mash))
  for (j in 1:ncol(lfsr.mash)) {
    sig.row=which(lfsr.mash[,i]<thresh)
    sig.col=which(lfsr.mash[,j]<thresh)
    a=(union(sig.row,sig.col))
    quotient=(pm.mash.beta[a,i]/pm.mash.beta[a,j])
    shared.fold.size[i,j] = mean(quotient > 0.5 & quotient < 2)
  }

clrs <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                               "#E0F3F8","#91BFDB","#4575B4")))(64)

all.tissue.order <- c("Colon−Sigmoid.txt","Colon−Transverse.txt","Stomach.txt","Esophagus−Gastroesophageal−Junction.txt","Esophagus−Muscularis.txt","Uterus.txt",
                      "Prostate.txt","Vagina.txt","Nerve−Tibial.txt","Breast−Mammary−Tissue.txt","Adipose−Subcutaneous.txt",
                      "Adipose−Visceral−−Omentum.txt","Artery−Coronary.txt","Artery−Aorta.txt","Artery−Tibial.txt",
                      "Small-Intestine-Terminal-Ileum.txt","Ovary.txt","Adrenal-Gland.txt","Pituitary.txt","Spleen.txt",
                      "Esophagus−Mucosa.txt","Liver.txt","Minor−Salivary−Gland.txt",
                      "Heart−Atrial−Appendage.txt","Heart−Left−Ventricle.txt","Pancreas.txt","Cells−Cultured−fibroblasts.txt","Lung.txt",
                      "Skin−Not−Sun−Exposed−−Suprapubic.txt","Skin−Sun−Exposed−−Lower−leg.txt","Thyroid.txt","Cells−EBV−transformed−lymphocytes.txt",
                      "Muscle−Skeletal.txt",
                      "Whole−Blood.txt","Testis.txt","Brain-Cerebellum.txt","Brain-Cerebellar-Hemisphere.txt",
                      "Brain-Spinal-cord--cervical-c-1.txt","Brain-Substantia-nigra.txt","Brain-Frontal-Cortex--BA9.txt",
                      "Brain-Cortex.txt","Brain-Hippocampus.txt","Brain-Hypothalamus.txt","Brain-Amygdala.txt",
                      "Brain-Putamen--basal-ganglia.txt","Brain-Anterior-cingulate-cortex--BA24.txt",
                      "Brain-Nucleus-accumbens--basal-ganglia.txt","Brain-Caudate--basal-ganglia.txt")
all.tissue.order <- gsub(all.tissue.order,pattern = "−",replacement = "-")
all.tissue.order <- gsub(all.tissue.order,pattern = "−−",replacement = "--")
lat=shared.fold.size[all.tissue.order,all.tissue.order]
colnames(lat) <- gsub(colnames(lat),pattern = "\\.txt",replacement = "")
colnames(lat) <- gsub(colnames(lat),pattern = "-",replacement = " ")
rownames(lat) <- gsub(rownames(lat),pattern = "\\.txt",replacement = "")
rownames(lat) <- gsub(rownames(lat),pattern = "-",replacement = " ")
pdf("./GTEX/08.analysis/MASH/heatmap.pdf",width=15,height =15)

p <- corrplot(lat,
              type = "upper",
              title = "aa",
              is.corr = FALSE,
              order = "original",
              col.lim = c(0.2,1), 
              col = colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                                           "#E0F3F8","#91BFDB","#4575B4")))(64),
              method = "shade",
              tl.col = "black",
              diag = TRUE)
print(p)
dev.off()


#sign
rm(list = ls())
out      <- readRDS("./GTEX/08.analysis/MASH/test/gtexresults-master/fastqtl_to_mash_output/tissue.mash.rds")
maxb     <- out$strong.b
maxz     <- out$strong.z
out      <-readRDS("./GTEX/08.analysis/MASH/test/gtexresults-master/mashr_flashr_workflow_output/tissue.mash.EZ.posterior.rds")
pm.mash  <- out$PosteriorMean
lfsr.all <- out$lfsr
standard.error <- maxb/maxz
pm.mash.beta <- pm.mash*standard.error

thresh=0.05
pm.mash.beta=pm.mash.beta[rowSums(lfsr.all<0.05)>0,]
lfsr.mash=lfsr.all[rowSums(lfsr.all<0.05)>0,]
shared.fold.size=matrix(NA,nrow = ncol(lfsr.mash),ncol=ncol(lfsr.mash))
colnames(shared.fold.size)=rownames(shared.fold.size)=colnames(maxz)
for(i in 1:ncol(lfsr.mash)){
  for(j in 1:ncol(lfsr.mash)){
    sig.row=which(lfsr.mash[,i]<thresh)
    sig.col=which(lfsr.mash[,j]<thresh)
    a=(union(sig.row,sig.col))
    quotient=(pm.mash.beta[a,i]/pm.mash.beta[a,j])
    shared.fold.size[i,j]=mean(quotient > 0)
  }
}

all.tissue.order <- c("Colon−Sigmoid.txt","Esophagus−Gastroesophageal−Junction.txt","Esophagus−Muscularis.txt","Uterus.txt",
                      "Prostate.txt","Vagina.txt","Nerve−Tibial.txt","Breast−Mammary−Tissue.txt","Adipose−Subcutaneous.txt",
                      "Adipose−Visceral−−Omentum.txt","Artery−Coronary.txt","Artery−Aorta.txt","Artery−Tibial.txt",
                      "Stomach.txt","Colon−Transverse.txt","Small-Intestine-Terminal-Ileum.txt","Ovary.txt",
                      "Adrenal-Gland.txt","Pituitary.txt","Spleen.txt",
                      "Esophagus−Mucosa.txt","Liver.txt","Minor−Salivary−Gland.txt",
                      "Pancreas.txt","Heart−Atrial−Appendage.txt","Heart−Left−Ventricle.txt","Cells−Cultured−fibroblasts.txt","Lung.txt",
                      "Skin−Not−Sun−Exposed−−Suprapubic.txt","Skin−Sun−Exposed−−Lower−leg.txt","Thyroid.txt","Cells−EBV−transformed−lymphocytes.txt",
                      "Muscle−Skeletal.txt",
                      "Whole−Blood.txt","Testis.txt","Brain-Cerebellum.txt","Brain-Cerebellar-Hemisphere.txt",
                      "Brain-Spinal-cord--cervical-c-1.txt","Brain-Substantia-nigra.txt","Brain-Frontal-Cortex--BA9.txt",
                      "Brain-Cortex.txt","Brain-Hippocampus.txt","Brain-Hypothalamus.txt","Brain-Amygdala.txt",
                      "Brain-Putamen--basal-ganglia.txt","Brain-Anterior-cingulate-cortex--BA24.txt",
                      "Brain-Nucleus-accumbens--basal-ganglia.txt","Brain-Caudate--basal-ganglia.txt")
all.tissue.order <- gsub(all.tissue.order,pattern = "−",replacement = "-")
all.tissue.order <- gsub(all.tissue.order,pattern = "−−",replacement = "--")
lat=shared.fold.size[all.tissue.order,all.tissue.order]

clrs <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                               "#E0F3F8","#91BFDB","#4575B4")))(64)

colnames(lat) <- gsub(colnames(lat),pattern = "\\.txt",replacement = "")
colnames(lat) <- gsub(colnames(lat),pattern = "-",replacement = " ")
rownames(lat) <- gsub(rownames(lat),pattern = "\\.txt",replacement = " ")
rownames(lat) <- gsub(rownames(lat),pattern = "-",replacement = "")
pdf("./sign_heatmap.pdf",width=15,height = 15)
p <- corrplot(lat,
              type = "upper",
              title = "\n\n\n",
              is.corr = FALSE,
              order = "original",
              col.lim = c(0.94,1), 
              col = colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                                           "#E0F3F8","#91BFDB","#4575B4")))(64),
              method = "shade",
              tl.col = "black",
              diag = TRUE)
print(p)
dev.off()

range(lat)
#Number of tissue-specific eQTLs in each tissue
source("./gtexresults-master/code/normfuncs.R")
thresh <- 0.05
out      <- readRDS("./gtexresults-master/fastqtl_to_mash_output/tissue.mash.rds")
maxbeta <- out$strong.b
maxz    <- out$strong.z
standard.error <- out$strong.s
out      <-readRDS("./gtexresults-master/mashr_flashr_workflow_output/tissue.mash.EZ.posterior.rds")
pm.mash        <- out$PosteriorMean
pm.mash.beta   <- pm.mash*standard.error
lfsr           <- out$lfsr
lfsr[lfsr < 0] <- 0

sigmat <- (lfsr <= thresh)
nsig   <- rowSums(sigmat)
pdf("./mag.pdf",width=8,height = 5)
hist((het.func(het.norm(effectsize=pm.mash.beta[nsig>0,]),threshold=0.5)),
     main="",xlab="",breaks=0.5:48.5,col="grey",freq=FALSE,ylim=c(0,0.2),
     xaxt="n")
axis(1,at = seq(1, 48, by=1),labels = c(1:48))
mtext("All Tissues")
dev.off()


sign.func <- function (normeffectsize)
  apply(normeffectsize,1,function(x)(sum(x>0)))
sigmat <- (lfsr<=thresh)
nsig   <- rowSums(sigmat)
pdf("./sign.pdf",width=8,height = 5)
hist(sign.func(het.norm(effectsize=pm.mash.beta[nsig>0,])),main="",xlab="",
     breaks=0.5:48.5,col="grey",freq=FALSE,xaxt="n",ylim=c(0,0.9))
axis(1, at=seq(1, 48, by=1), labels=c(1:48))
mtext("Number of tissues shared by sign")
dev.off()


rm(list = ls())
source("./gtexresults-master/code/normfuncs.R")
thresh <- 0.05

out            <- readRDS("./gtexresults-master/fastqtl_to_mash_output/tissue.mash.rds")
maxbeta        <- out$strong.b
maxz           <- out$strong.z
standard.error <- out$strong.s
out            <- readRDS("./gtexresults-master/mashr_flashr_workflow_output/tissue.mash.EZ.posterior.rds")
pm.mash        <- out$PosteriorMean
pm.mash.beta   <- pm.mash*standard.error
lfsr           <- out$lfsr
lfsr[lfsr < 0] <- 0
tissue.names   <- colnames(lfsr)
colnames(lfsr) <- tissue.names

sigmat <- (lfsr<=thresh)
nsig   <- rowSums(sigmat)
sigall <- mean(het.norm(pm.mash.beta[nsig>0,])>0)

sigmat <- (lfsr[,-c(7:19)]<=thresh)
nsig   <- rowSums(sigmat)
sigall.nobrain <- mean(het.norm(pm.mash.beta[nsig>0,-c(7:19)]) > 0)

sigmat <- (lfsr[,c(7:19)]<=thresh)
nsig   <- rowSums(sigmat)
sigall.brainonly <- mean(het.norm(pm.mash.beta[nsig>0,c(7:19)]) > 0)

out_nobrain            <- readRDS("./gtexresults-master/fastqtl_to_mash_output/no_brain.mash.rds")
maxbeta_nobrain        <- out_nobrain$strong.b
maxz_nobrain           <- out_nobrain$strong.z
standard.errorno_brain <- out_nobrain$strong.s
out_nobrain            <- readRDS("./gtexresults-master/mashr_flashr_workflow_output/no_brain.mash.EZ.posterior.rds")
pm.mash_nobrain         <- out_nobrain$PosteriorMean
pm.mash.beta_nobrain    <- pm.mash_nobrain *standard.errorno_brain

lfsr.nobrain           <- out_nobrain$lfsr
lfsr.nobrain[lfsr.nobrain < 0] <- 0
tissue.names_nobrain   <- colnames(lfsr.nobrain)
colnames(lfsr.nobrain) <- tissue.names_nobrain

sigmat     <- (lfsr.nobrain<=thresh)
nsig       <- rowSums(sigmat)
signobrain <- mean(het.norm(pm.mash.beta_nobrain[nsig>0,]) > 0)


out_brain            <- readRDS("./gtexresults-master/fastqtl_to_mash_output/brain.mash.rds")
maxbeta_brain        <- out_brain$strong.b
maxz_brain           <- out_brain$strong.z
standard.error_brain <- out_brain$strong.s
out_brain            <- readRDS("./gtexresults-master/mashr_flashr_workflow_output/brain.mash.EZ.posterior.rds")
pm.mash_brain         <- out_brain$PosteriorMean
pm.mash.brain.only    <- pm.mash_brain *standard.error_brain

lfsr.brain.only           <- out_nobrain$lfsr
lfsr.brain.only[lfsr.brain.only < 0] <- 0
tissue.names_brain   <- colnames(lfsr.brain.only)
colnames(lfsr.brain.only) <- tissue.names_brain

sigmat   <- (lfsr.brain.only<=thresh)
nsig     <- rowSums(sigmat)
sigbrainonly <- mean(het.norm(pm.mash.brain.only[nsig>0,]) > 0)

sigmat <- (lfsr<=thresh)
nsig   <- rowSums(sigmat)
magall <- mean(het.norm(pm.mash.beta[nsig>0,])>0.5)

sigmat <- (lfsr[,-c(7:19)]<=thresh)
nsig   <- rowSums(sigmat)
magall.excludingbrain <- mean(het.norm(pm.mash.beta[nsig>0,-c(7:19)]) > 0.5)

sigmat <- (lfsr[,c(7:19)]<=thresh)
nsig   <- rowSums(sigmat)
magall.brainonly <- mean(het.norm(pm.mash.beta[nsig>0,c(7:19)]) > 0.5)

sigmat     <- (lfsr.nobrain<=thresh)
nsig       <- rowSums(sigmat)
magnobrain <- mean(het.norm(pm.mash.beta_nobrain[nsig>0,]) > 0.5)

sigmat   <- (lfsr.brain.only<=thresh)
nsig     <- rowSums(sigmat)
magbrain <- mean(het.norm(pm.mash.brain.only[nsig>0,]) > 0.5)



round(matrix(rbind(c(sigall,sigall.nobrain,signobrain,
                     sigall.brainonly,sigbrainonly),
                   c(magall,magall.excludingbrain,magnobrain,
                     magall.brainonly,magbrain)),
             nrow = 2,ncol = 5,
             dimnames = list(c("shared by sign","shared by magnitude"),
                             c("all tissues","non-brain","(non-brain)",
                               "brain","(brain)"))),
      digits = 3)


rm(list = ls())
source("./gtexresults-master/code/normfuncs.R")
thresh <- 0.05
out            <- readRDS("./gtexresults-master/fastqtl_to_mash_output/tissue.mash.rds")
maxbeta        <- out$strong.b
maxz           <- out$strong.z
standard.error <- out$strong.s
out            <- readRDS("./gtexresults-master/mashr_flashr_workflow_output/tissue.mash.EZ.posterior.rds")
pm.mash        <- out$PosteriorMean
lfsr.mash    <- out$lfsr
pm.mash.beta   <- pm.mash*standard.error
nsig                <- rowSums(lfsr.mash < thresh)
pm.mash.beta.norm   <- het.norm(effectsize = pm.mash.beta)
pm.mash.beta.norm   <- pm.mash.beta.norm[nsig > 0,]
lfsr.mash           <- as.matrix(lfsr.mash[nsig > 0,])
colnames(lfsr.mash) <- colnames(maxz)
a         <- which(rowSums(pm.mash.beta.norm > 0.5) == 1)
lfsr.fold <- as.matrix(lfsr.mash[a,])
pm        <- as.matrix(pm.mash.beta.norm[a,])
tspec     <- NULL
for(i in 1:ncol(pm))
  tspec[i] <- sum(pm[,i] > 0.5)
tspec           <- as.matrix(tspec)
rownames(tspec) <- colnames(maxz)
par(mfrow = c(2,1))
barplot(as.numeric(t(tspec)),las = 2,cex.names = 0.75,col = col,
        names = colnames(lfsr.fold))



