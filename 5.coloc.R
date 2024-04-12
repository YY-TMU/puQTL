.libPaths("./R/x86_64-redhat-linux-gnu-library/4.0")
library(coloc)

rm(list = ls())
canner_data <- read.delim("./GTEX_new_sample/5.analysis/puQTL_summary.txt",header = T)
canner_data <- canner_data[,c(1,4)]
colnames(canner_data) <- c("cancer","number")
canner_data$cancer <- gsub(canner_data$cancer,pattern = " ",replacement = "")
setwd("./GTEX_new_sample/5.analysis/Coloc/0.rawdata/puQTL")

cancer_list <- list.files()
setwd("./coloc/add_GWAS_data/sig_add_version")
gwas_list <- list.files(pattern = "txt")

sample_size <- read.delim("./coloc/sample.txt",header = F)
for (k in 1:length(gwas_list)) {
  #format full summaty of GWAS 
  gwas_data <- read.delim(paste0("./coloc/add_GWAS_data/fin_sum/",gwas_list[k]),header = T)
  gwas_data <- gwas_data[,c(1,5,6,7)]
  colnames(gwas_data) <- c("snp","beta","se","pval_nominal_gwas")
  if (nrow(gwas_data) >2) {
    gwas_data$se <- as.numeric(gwas_data$se)
    gwas_data$varbeta <- (gwas_data$se)^2
    gwas_data <- gwas_data[!is.na(gwas_data$pval_nominal_gwas),]
    gwas_data <- gwas_data[,c("snp","pval_nominal_gwas","beta","varbeta")]
    tmp_name <- gsub(gwas_list[k],pattern = ".txt",replacement = "")
    tmp_name <- gsub(tmp_name,pattern = "count",replacement = "Count")
    GWAS_size <- sample_size[sample_size$V2 == tmp_name,]$V1
    cancer_list <- tissue
    cancer_size <- canner_data[canner_data$cancer == cancer_list,]$number
    setwd(paste0("./GTEX_new_sample/5.analysis/Coloc/0.rawdata/puQTL/",cancer_list))
    promoter_list <- list.files()
    for (m in 1:length(promoter_list)) {
      #full summaty of puQTL 
      promoter_data <- read.delim(paste0("./GTEX_new_sample/5.analysis/Coloc/0.rawdata/puQTL/",cancer_list,"/",promoter_list[m]),header = F)
      promoter_FDR <- promoter_data[,c(1,6)]
      promoter_FDR <- promoter_FDR[!duplicated(promoter_FDR$V1),]
      colnames(promoter_FDR) <- c("snp","FDR")
      promoter_data <- promoter_data[,c(1,2,3,5,7)]
      colnames(promoter_data) <- c("snp","promoter","eqtl_beta","pval_nominal_eqtl","MAF")
      promoter_data <- promoter_data[!duplicated(promoter_data$snp),]
      colnames(promoter_data) <- c("snp","promoter","eqtl_beta","pval_nominal_eqtl","maf")
      input <- dplyr::inner_join(gwas_data,promoter_data,by="snp")
      input <- input[!duplicated(input$snp),]
      if (nrow(input) > 0) {
        colnames(input)[1] <- "rs_id"
        #coloc analysis
        result <- coloc.abf(dataset1=list(snp = input$rs_id,pvalues=input$pval_nominal_eqtl, type="quant", N=cancer_size),
                            dataset2=list(snp = input$rs_id,pvalues=input$pval_nominal_gwas,beta=input$beta,varbeta=input$varbeta, type="cc", N=GWAS_size),
                            MAF=input$maf)
        summary <- data.frame(result$summary)
        summary <- data.frame(t(summary))
        summary$traits <- gsub(gwas_list[k],pattern = "\\.txt",replacement = "")
        summary$promoter <- promoter_list[m]
        write.table(summary,paste0("./GTEX_new_sample/5.analysis/Coloc/full_result/0.rawdata/add/sum_",cancer_list),sep = "\t",quote = F,col.names = F,row.names = F,append = T)
        result <- result$results
        part_result <- result[result$SNP.PP.H4 > 0,]
        if (nrow(part_result) > 0) {
          part_result <- part_result[,c(1,14)]
          part_result$traits <- gsub(gwas_list[k],pattern =  "\\.txt",replacement = "")
          part_result$source <- "GWAS_catelog"
          part_result$first_p4 <- summary$PP.H4.abf
          part_result$promoter <- promoter_list[m]
          part_result <- dplyr::inner_join(part_result,promoter_FDR)
          write.table(part_result,paste0("./GTEX_new_sample/5.analysis/Coloc/full_result/0.rawdata/add/",cancer_list),sep = "\t",quote = F,col.names = F,row.names = F,append = T)
        }
      }
      
    } 
  }
}


