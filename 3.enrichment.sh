#snpeff analysis
#!/bin/bash 
for i in `ls ./GTEX/06.puQTL/significant`
do 
	        k=${i%%.*}
		java -Xmx8g -jar ./snpEff/snpEff.jar -c ./snpEff/snpEff.config GRCh37.75 ./GTEX/06.puQTL/vcf/${k}.vcf -verbose -stats ${k}.html -csvStats ${k}.csv > ./GTEX/08.analysis/snpEff/0.rawdata/${i}.vcf
	done
	
#shuffle analysis	
rm(list = ls())
library(stringr)
setwd("./GTEX/06.puQTL/vcf/")
tissue_list <- list.files(pattern = "")
for (m in 1:length(tissue_list)) {
  tissue_puQTL <- read.delim(paste0("./GTEX/06.puQTL/vcf/",tissue_list[m]),header = F)
  tissue_puQTL <- tissue_puQTL[,c(1,2,2,3)]
  tissue <- gsub(tissue_list[m],pattern = "\\.vcf",replacement = "")
  write.table(tissue_puQTL,paste0("./GTEX/08.analysis/shuffle/GTEx/",tissue),quote = F,row.names = F,col.names = F,sep = "\t")

}


for k in `ls ./GTEX/08.analysis/shuffle/GTEx`;
do
echo $k
echo -e "#! /bin/bash
    for ((i=1; i<=10000; i++))
    do
    bedtools shuffle -i ./shuffle/0.rawdata/promoter.bed -g ./ref_seq/Human/hg19.fa.fai -chrom > ./shuffle/0.rawdata/${k}_tmp.bed
    intersectBed -a ./GTEX/08.analysis/shuffle/GTEx/${k} -b ./shuffle/0.rawdata/${k}_tmp.bed -wa  | cut -f4 |sort | uniq | wc -l >> ./shuffle/0.rawdata/${k}.all.txt
    done
" > scripts/GTEx-${k}.sh
chmod a+x scripts/GTEx-${k}.sh
nohup ./scripts/GTEx-${k}.sh &>logs/GTEx-${k}.logs&
  done
  
  
#GREGOR analysis	
perl script/GREGOR.pl --conf example/Brain-Hippocampus.conf.conf
perl script/GREGOR.pl --conf example/Adrenal-Gland.conf
perl script/GREGOR.pl --conf example/Colon-Sigmoid.conf
perl script/GREGOR.pl --conf example/Esophagus-Gastroesophageal-Junction.conf
perl script/GREGOR.pl --conf example/Esophagus-Mucosa.conf
perl script/GREGOR.pl --conf example/Esophagus-Muscularis.conf
perl script/GREGOR.pl --conf example/Heart-Left-Ventricle.conf
perl script/GREGOR.pl --conf example/Liver.conf
perl script/GREGOR.pl --conf example/Lung.conf
perl script/GREGOR.pl --conf example/Ovary.conf
perl script/GREGOR.pl --conf example/Pancreas.conf
perl script/GREGOR.pl --conf example/Small-Intestine-Terminal-Ileum.conf
perl script/GREGOR.pl --conf example/Spleen.conf
perl script/GREGOR.pl --conf example/Stomach.conf

