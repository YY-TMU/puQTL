source activate ldsc


############ using S-LDSC for GWAS enrichment analysis
## make annotation

for i in cat tissue.txt
do
k=${i%%.*}
for chr in $(seq 1 22)
do
python ./ldsc/make_annot.py 
       --bed-file ./GTEX/06.puQTL/0.1_q_vaue/bed/$i 
       --bimfile ./1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr.bim 
       --annot-file ./GTEX/08.analysis/GWAS/LDSC_enrichment/FDR_puQTL_anno/$k.chr$chr.annot.gz
done
done


## estimate LD score

for i in cat tissue.txt
do
k=${i%%.*}
for chr in $(seq 1 22)
do
python ./ldsc/ldsc.py
      --l2 
      --bfile ./1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr 
      --ld-wind-cm 1 
      --annot ./$k.chr$chr.annot.gz 
      --thin-annot --out ./$k.chr$chr
      --print-snps ./ldsc/listHM3.txt
done
done

## Partition heritability

for trait in cat trait.txt
mkdir ./LDSC_enrichment/result_v1/${trait}
for i in `cat ./LDSC_enrichment/annot/tissue.txt`
do
python ./ldsc/ldsc.py 
       --h2 ./sum_LDSC_GWAS/${trait}.sumstats 
       --ref-ld-chr ./ldsc/baseline_gene_MAF_LD/baseline_gene_MAF_LD.,./LDSC_enrichment/annot/${i}.chr 
       --w-ld-chr ./ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. 
       --frqfile-chr ./ldsc/1000G_Phase3_frq/1000G.EUR.QC. 
       --overlap-annot 
       --print-cov --print-coefficients --print-delete-vals 
       --out ./LDSC_enrichment/result_v1/${trait}/${i}
done
done