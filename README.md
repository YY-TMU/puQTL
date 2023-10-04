# puQTL

## Description
Promoters play a crucial role in regulating gene transcription. However, our understanding of how genetic variants influence alternative promoter selection is still incomplete. In this study, we implemented a framework to identify genetic variants that affect the relative usage of alternative promoters, known as promoter usage quantitative trait loci (puQTLs). By constructing an atlas of human puQTLs across 48 different tissues from 503 individuals, we have identified approximately 303,025 genetic variants associated with promoter usage. These scripts were used for investigating puQTLs across 48 different tissues from GTEx Consortium.

## Data preprocess for puQTL analysis
The shell script 1.preprocess.sh was used to preprocess individual genotype data and promoter usage profiles for puQTL analysis. The principle components within the genotype, PEER factors, genderm and age were all prepared for identifying puQTLs.

## puQTL mapping for each tissue
The R script 2.MatrixeQTL.R was used to identify puQTLs for each tissue type. The R package Matrix eQTL was utlilized to test for association between genotype and promoter usage.

## Enrichment analysis for puQTLs
The shell script 3.enrichment.sh was used to evaluate enrichment of puQTLs in experimentally annotated epigenomic regulatory features.

## Tissue specificity analysis for puQTLs
The R script 3.mash.R utilized the mash method to elucidate the heterogeneity of puQTL effect sizes across different tissue types. 

## Disease/traits colocalization analysis
The R script 5.coloc.R was used for colocalization analysis of puQTLs with GWAS associated diseases or traits.

## Paritioned heritability
The shell script 6.LDSC.sh was used to compute partitioned heritability.


