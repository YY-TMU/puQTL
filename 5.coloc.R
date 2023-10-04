rm(list = ls())
setwd("./GTEX/06.puQTL/all_data_coloc")
tissue_list <- list.files(pattern = "txt")

anno_snp <- read.delim(paste0("./GTEX/01.preprocess/coloc.txt"),header = F)
colnames(anno_snp) <- c("snp","ID")

#format summary data 
for (k in 39:length(tissue_list)) {
  tissue <- gsub(tissue_list[k],pattern = "\\.txt",replacement = "")
  QTL_data <- read.delim(paste0("./GTEX/06.puQTL/significant/",tissue_list[k]),header = F)
  QTL_data <- QTL_data[,c(1,2,5)]
  colnames(QTL_data) <- c("snp","promoter","P-value")
  QTL_data <- dplyr::inner_join(QTL_data,anno_snp,by="snp")
  maf <- read.delim(paste0("./GTEX/01.imputation/genotype/tissue_maf/",tissue,".maf"),header = T)
  colnames(maf) <- c("snp","MAF")
  QTL_data <- dplyr::inner_join(QTL_data,maf,by="snp")
  write.table(QTL_data,paste0("./GTEX/08.analysis/GWAS/coloc/0.rawdata/significant_puQTL/",tissue_list[k]),sep = "\t",quote = F,col.names = F,row.names = F)
}

#coloc analysis in biobank data
library(coloc)
canner_data <- read.delim("./GTEX/08.analysis/GWAS/coloc/0.rawdata/sample_number",header = F)
colnames(canner_data) <- c("tissue","number")
trait <- read.delim("./coloc/biobank.txt",header = F)
colnames(trait) <- c("traits","gwas")
trait$traits <- paste0(trait$traits,".gwas.imputed_v3.both_sexes")
setwd("./GTEX/08.analysis/GWAS/coloc/0.rawdata/significant_puQTL")
tissue_list <- list.files()
for (i in 1:length(tissue_list)) {
  qtl_data <- read.delim(paste0("./GTEX/08.analysis/GWAS/coloc/0.rawdata/significant_puQTL/",tissue_list[i]),header = F)
  qtl_data <- qtl_data[,c(4,2,3,5)]
  colnames(qtl_data) <- c("rs_id","promoter","pval_nominal_eqtl","maf")
  setwd("./coloc/biobank")
  gwas_list <- list.files(pattern = "tsv")
  for (k in 1:length(gwas_list)) {
    gwas_data <- read.delim(paste0("./coloc/biobank/sig/",gwas_list[k]),header = T)
    colnames(gwas_data) <- c("rs_id","beta","se","pval_nominal_gwas")
    gwas_data$se <- as.numeric(gwas_data$se)
    gwas_data$varbeta <- (gwas_data$se)^2
    gwas_data <- gwas_data[,c("rs_id","pval_nominal_gwas","beta","varbeta")]
    tissue_size <- canner_data[canner_data$tissue == tissue_list[i],]$number
    promoter_data <- read.delim(paste0("./GTEX/06.puQTL/significant/",tissue_list[i]),header = F)
    qtl_data <- qtl_data[qtl_data$promoter %in%  promoter_data$V2,]
    coloc_data <- dplyr::inner_join(gwas_data,qtl_data,by="rs_id")
    if (nrow(coloc_data) > 0) {
      coloc_list <- split(coloc_data,f=coloc_data$promoter)
      for (g in 1:length(coloc_list)) {
        input <- coloc_list[[g]]
        input$beta <- as.numeric(input$beta)
        input$pval_nominal_eqtl <- as.numeric(input$pval_nominal_eqtl)
        input$pval_nominal_gwas <- as.numeric(input$pval_nominal_gwas)
        input$maf <- as.numeric(input$maf)
        input$snp <- c(1:nrow(input))
        input$snp <- paste0("SNP.",input$snp)
        result <- coloc.abf(dataset1=list(pvalues=input$pval_nominal_eqtl, type="quant", N=tissue_size),
                            dataset2=list(beta=input$beta,varbeta=input$varbeta, type="cc", N=337119),
                            MAF=input$maf)
        summary <- data.frame(result$summary)
        summary <- data.frame(t(summary))
        summary$traits <- gsub(gwas_list[k],pattern = "\\.tsv",replacement = "")
        summary <- dplyr::inner_join(summary,trait,by="traits")
        summary$promoter <- names(coloc_list[g])
        write.table(summary,paste0("./08.analysis/GWAS/coloc/1.rawresult/significat_puQTL/biobank/sum_",tissue_list[i]),sep = "\t",quote = F,col.names = F,row.names = F,append = T)
        result <- result$results
        part_result <- result[result$SNP.PP.H4 > 0,]
        if (nrow(part_result) > 0) {
          part_result <- part_result[,c(1,14)]
          part_result <- dplyr::inner_join(part_result,input,by="snp")
          part_result <- part_result[,-1]
          part_result$traits <- gsub(gwas_list[k],pattern = "\\.tsv",replacement = "")
          part_result$source <- "GWAS_catelog"
          part_result <- dplyr::inner_join(part_result,trait,by="traits")
          
          write.table(part_result,paste0("./GTEX/08.analysis/GWAS/coloc/1.rawresult/significat_puQTL/biobank/",tissue_list[i]),sep = "\t",quote = F,col.names = F,row.names = F,append = T)
        }
      }
    }
    
  }}


#coloc analysis in add GWAS data
library(coloc)

canner_data <- read.delim("./08.analysis/GWAS/coloc/0.rawdata/sample_number",header = F)
colnames(canner_data) <- c("tissue","number")
setwd("./08.analysis/GWAS/coloc/0.rawdata/significant_puQTL")
tissue_list <- list.files()
sample_size <- read.delim("./coloc/bb.txt",header = F)
for (i in 1:length(tissue_list)) {
  qtl_data <- read.delim(paste0("./GTEX/08.analysis/GWAS/coloc/0.rawdata/significant_puQTL/",tissue_list[i]),header = F)
  qtl_data <- qtl_data[,c(1,2,3,5)]
  colnames(qtl_data) <- c("rs_id","promoter","pval_nominal_eqtl","maf")
  setwd("./coloc/add_GWAS_data/sig_add_version")
  gwas_list <- list.files(pattern = "txt")
  for (k in 1:length(gwas_list)) {
    gwas_data <- read.delim(paste0("./coloc/add_GWAS_data/sig_add_version/",gwas_list[k]),header = T)
    colnames(gwas_data) <- c("rs_id","beta","se","pval_nominal_gwas")
    gwas_data$se <- as.numeric(gwas_data$se)
    gwas_data$varbeta <- (gwas_data$se)^2
    gwas_data <- gwas_data[,c("rs_id","pval_nominal_gwas","beta","varbeta")]
    tissue_size <- canner_data[canner_data$tissue == tissue_list[i],]$number
    tmp_name <- gsub(gwas_list[k],pattern = ".txt",replacement = "")
    tmp_name <- gsub(tmp_name,pattern = "count",replacement = "Count")
    
    GWAS_size <- sample_size[sample_size$V2 == tmp_name,]$V1
    
    promoter_data <- read.delim(paste0("./GTEX/06.puQTL/significant/",tissue_list[i]),header = F)
    qtl_data <- qtl_data[qtl_data$promoter %in%  promoter_data$V2,]
    coloc_data <- dplyr::inner_join(gwas_data,qtl_data,by="rs_id")
    if (nrow(coloc_data) > 0) {
      coloc_list <- split(coloc_data,f=coloc_data$promoter)
      for (g in 1:length(coloc_list)) {
        input <- coloc_list[[g]]
        input$beta <- as.numeric(input$beta)
        input$pval_nominal_eqtl <- as.numeric(input$pval_nominal_eqtl)
        input$pval_nominal_gwas <- as.numeric(input$pval_nominal_gwas)
        input$maf <- as.numeric(input$maf)
        input$snp <- c(1:nrow(input))
        input$snp <- paste0("SNP.",input$snp)
        result <- coloc.abf(dataset1=list(pvalues=input$pval_nominal_eqtl, type="quant", N=tissue_size),
                            dataset2=list(beta=input$beta,varbeta=input$varbeta, type="cc", N=GWAS_size),
                            MAF=input$maf)
        summary <- data.frame(result$summary)
        summary <- data.frame(t(summary))
        summary$traits <- gsub(gwas_list[k],pattern = "\\.txt",replacement = "")
        summary$promoter <- names(coloc_list[g])
        write.table(summary,paste0("./GTEX/08.analysis/GWAS/coloc/1.rawresult/significat_puQTL_version/add/sum_",tissue_list[i]),sep = "\t",quote = F,col.names = F,row.names = F,append = T)
        result <- result$results
        part_result <- result[result$SNP.PP.H4 > 0,]
        if (nrow(part_result) > 0) {
          part_result <- part_result[,c(1,14)]
          part_result <- dplyr::inner_join(part_result,input,by="snp")
          part_result <- part_result[,-1]
          part_result$traits <- gsub(gwas_list[k],pattern = "\\.txt",replacement = "")
          part_result$source <- "GWAS_catelog"
          
          write.table(part_result,paste0("./GTEX/08.analysis/GWAS/coloc/1.rawresult/significat_puQTL_version/add/",tissue_list[i]),sep = "\t",quote = F,col.names = F,row.names = F,append = T)
        }
      }
    }
    
  }
  
  
}

gene_data <- read.delim("./ref_seq/Human/gene_trans_anno.txt",header = F)
colnames(gene_data) <- c("gene","type","name")
setwd("./GTEX/08.analysis/GWAS/coloc/1.rawresult/significat_puQTL_version/add/")
tissue_list <- list.files(pattern = "sum")
for (i in 1:length(tissue_list)) {
  qtl_data <- read.delim(paste0("./GTEX/08.analysis/GWAS/coloc/1.rawresult/significat_puQTL_version/add/",tissue_list[i]),header = F)
  qtl_data$gene <- gsub(qtl_data$V8,pattern = ".*ENSG",replacement = "ENSG")
  qtl_data$gene <- gsub(qtl_data$gene,pattern = "\\..*",replacement = "")
  qtl_data <- dplyr::inner_join(qtl_data,gene_data)
  write.table(qtl_data,paste0("./GTEX/08.analysis/GWAS/coloc/1.rawresult/significat_puQTL_version/add/",tissue_list[i]),sep = "\t",quote = F,col.names = F,row.names = F)
}