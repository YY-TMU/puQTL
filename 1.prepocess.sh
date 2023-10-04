#the quality control step initially defined by the GTEx Consortium
./GTEx_v8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf
vcftools --vcf ./GTEx_v8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf --out ./GTEX/01.preprocess/1.rawdata/0.all --plink
./plink --file ./GTEX/01.preprocess/1.rawdata/0.all --geno 0.05 --maf 0.01 --hwe 0.000001 --recode --out ./GTEX/01.preprocess/1.rawdata/prefilter
./plink --file ./GTEX/01.preprocess/1.rawdata/prefilter --recode vcf --out ./GTEX/01.preprocess/1.rawdata/2.filter

grep "rs" ./GTEX/01.imputation/promoter_activity/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt | awk '{split($8,a,"_"); print $1"\t"a[1]"\t"a[2]"\t"$2"\t"$7}'  |sed 1d | grep -v "NA" | grep -v "GL" > ./GTEX/01.preprocess/1.rawdata/snp_38to37 
grep -v "#" ./GTEX/01.preprocess/1.rawdata/2.filter.vcf > ./GTEX/01.preprocess/1.rawdata/2.filter_tmp.vcf
python ./GTEX/01.preprocess/1.rawdata/bin/coord_trans.py /./GTEX/01.preprocess/1.rawdata/2.filter_tmp.vcf ./GTEX/01.preprocess/1.rawdata/snp_38to37 ./GTEX/01.preprocess/1.rawdata/3.filter_37.vcf

grep "#" 2.filter.vcf > 2.filter.vcf_header
cat 2.filter.vcf_header 3.filter_37.vcf > 4.all.vcf

for i in `ls ./GTEX/02.sample/tissue`;
do
echo $i
awk '{print $1"_"$1"_"$1"_"$1}' ./GTEX/01.imputation/overlap/${i}.samples > ./GTEX/01.preprocess/1.rawdata/5.tissue_label/${i}
done


vcftools --vcf ./GTEX/01.preprocess/1.rawdata/4.all.vcf --out ./GTEX/01.preprocess/1.rawdata/5.all --plink
./plink --file ./GTEX/01.preprocess/1.rawdata/5.all --recode vcf --out ./GTEX/01.preprocess/1.rawdata/6.all


for i in `ls ./GTEX/01.preprocess/1.rawdata/5.tissue_label`
do
echo $i
echo "#! /bin/bash
bcftools view -S ./GTEX/01.preprocess/1.rawdata/5.tissue_label/${i} ./GTEX/01.preprocess/1.rawdata/6.all.vcf > /./GTEX/01.preprocess/2.tissue_vcf/${i}.vcf
bcftools view -i \"INFO/AC>10\" ./GTEX/01.preprocess/2.tissue_vcf/${i}.vcf > ./GTEX/01.preprocess/3.filter_vcf/${i}.vcf 
bcftools view --types snps -m 2 -M 2 -q 0.01:minor ./GTEX/01.preprocess/3.filter_vcf/${i}.vcf > ./GTEX/01.preprocess/4.PCA_vcf/${i}.vcf
vcftools --vcf ./GTEX/01.preprocess/4.PCA_vcf/${i}.vcf --out ./GTEX/01.preprocess/5.plink_format/${i} --plink" > ./GTEX/01.preprocess/2.tissue_vcf/script/$i.sh
chmod a+x ./GTEX/01.preprocess/2.tissue_vcf/script/$i.sh
nohup ./GTEX/01.preprocess/2.tissue_vcf/script/$i.sh &> ./GTEX/01.preprocess/2.tissue_vcf/script/$i.log&
  done


##
for i in `ls ./GTEX/01.preprocess/1.rawdata/5.tissue_label`
do
echo $i
echo "#! /bin/bash
bcftools view -i \"INFO/AC>10\" ./GTEX/01.preprocess/2.tissue_vcf/${i}.vcf > ./GTEX/01.preprocess/3.filter_vcf/${i}.vcf 
bcftools view --types snps -m 2 -M 2 -q 0.01:minor ./GTEX/01.preprocess/3.filter_vcf/${i}.vcf > ./GTEX/01.preprocess/4.PCA_vcf/${i}.vcf
vcftools --vcf ./GTEX/01.preprocess/4.PCA_vcf/${i}.vcf --out ./GTEX/01.preprocess/5.plink_format/${i} --plink" > ./GTEX/01.preprocess/2.tissue_vcf/script/filter_${i}.sh 
chmod a+x ./GTEX/01.preprocess/2.tissue_vcf/script/filter_${i}.sh 
nohup ./GTEX/01.preprocess/2.tissue_vcf/script/filter_${i}.sh  &> ./GTEX/01.preprocess/2.tissue_vcf/script/filter_${i}.log&
  done


for i in `ls ./GTEX/01.preprocess/1.rawdata/5.tissue_label`
do
echo $i
echo "#! /bin/bash
bcftools view --types snps -m 2 -M 2 -q 0.01:minor .GTEX/01.preprocess/3.filter_vcf/${i}.vcf > ./GTEX/01.preprocess/4.PCA_vcf/${i}.vcf " > ./GTEX/01.preprocess/2.tissue_vcf/script/pca_filter_${i}.sh 
chmod a+x ./GTEX/01.preprocess/2.tissue_vcf/script/pca_filter_${i}.sh  
nohup ./GTEX/01.preprocess/2.tissue_vcf/script/pca_filter_${i}.sh   &> ./GTEX/01.preprocess/2.tissue_vcf/script/pca_filter_${i}.log&
  done

#convect vcf to plink format
for i in `ls ./GTEX/01.preprocess/1.rawdata/5.tissue_label`
do
echo $i
echo -e "#! /bin/bash
    vcftools --vcf ./GTEX/01.preprocess/4.PCA_vcf/${i}.vcf --out ./GTEX/01.preprocess/5.plink_format/${i} --plink" > scripts/${i}.format.sh
chmod a+x scripts/${i}.format.sh
nohup ./scripts/${i}.format.sh &>logs/${i}.format.logs&
  done

#caculate PCA
for i in `ls ./GTEX/01.preprocess/1.rawdata/5.tissue_label`
do
./plink --noweb --file ./GTEX/01.preprocess/5.plink_format/${i} --pca 5 --out ./GTEX/01.preprocess/5.eigen/${i}
echo ${i}
done

for i in `ls ./GTEX/01.preprocess/1.rawdata/5.tissue_label`
do
echo $i
echo -e "#! /bin/bash
    grep -v \"##\" ./01.preprocess/4.PCA_vcf/${i}.vcf > ./GTEX/01.preprocess/6.snp_matrix/${i}_tmp.vcf
    python ./pu_act/3.sample_name/bin/parse.py -v ./GTEX/01.preprocess/6.snp_matrix/${i}_tmp.vcf -o ./GTEX/01.preprocess/6.snp_matrix/${i}.snpmatrix
    rm ./01.preprocess/6.snp_matrix/${i}_tmp.vcf" >./GTEX/01.preprocess/6.snp_matrix/scripts/$i.sh
chmod a+x ./GTEX/01.preprocess/6.snp_matrix/scripts/$i.sh
nohup ./GTEX/01.preprocess/6.snp_matrix/scripts/$i.sh &> ./GTEX/01.preprocess/6.snp_matrix/scripts/$i.log&
  done


for i in `ls ./GTEX/02.sample/tissue`;
do
echo $i
head -1 ./GTEX/01.preprocess/6.snp_matrix/${i}.snpmatrix | sed 's/\t/\n/g' | sed '1d' > ./GTEX/03.combine_data/matrix_sampleIDs/${i}.samples
done





# make covariates
mkdir covariates
for i in `ls ./GTEX/02.sample/tissue`;
do
echo $i
Rscript ./1.combine_data/bin/format.covariates.R ./GTEX/03.combine_data/covariates/${i}.covariates ./GTEX/03.combine_data/matrix_sampleIDs/${i}.samples ./GTEX/05.analysis_data/covariates/${i}.tmp
done

for i in `ls ./GTEX/02.sample/tissue`;
do
echo $i
sed 's/female/0/g' ./GTEX/05.analysis_data/covariates/${i}.tmp | sed 's/male/1/g' | sed 's/\[Not Available\]/NA/g' > ./GTEX/05.analysis_data/covariates/${i}.covariates
rm -rf ./GTEX/05.analysis_data/covariates/${i}.tmp
done

#make  expression
mkdir expression
for i in `ls ./GTEX/02.sample/tissue`;
do
echo $i
Rscript ./1.combine_data/bin/format.exp.R ./GTEX/03.combine_data/promoter_usage/${i}.txt ./GTEX/03.combine_data/matrix_sampleIDs/${i}.samples ./GTEX/05.analysis_data/promoter_usage/${i}
done

# make eigen

for i in `ls ./GTEX/02.sample/tissue`;
do
echo $i
cut -f1,3- ./GTEX/01.preprocess/5.eigen/${i}.eigenvec > ./GTEX/01.preprocess/5.eigen/${i}.tmp
done

for i in `ls | grep tmp`; do awk -F " " '{split($1,a,"_"); print a[1]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' /./GTEX/01.preprocess/5.eigen/$i > ${i}_tmp; rm $i; mv ${i}_tmp ${i}; done

for i in `ls ./GTEX/02.sample/tissue`;
do
echo $i
Rscript ./1.combine_data/bin/format.eigen.R ./GTEX/01.preprocess/5.eigen/${i}.tmp ./GTEX/03.combine_data/matrix_sampleIDs/${i}.samples ./GTEX/05.analysis_data/eigen/${i}.eigen
done

rm ./GTEX/01.preprocess/5.eigen/*.tmp



# make peer

for i in `ls ./GTEX/02.sample/tissue`;
do
echo $i
Rscript ./1.combine_data/bin/format.peer.R ./GTEX/01.imputation/promoter_activity/4.peer/${i}.txt ./GTEX/03.combine_data/matrix_sampleIDs/${i}.samples ./GTEX/05.analysis_data/peer/${i}.txt
done
# concate covariates

mkdir combine
for i in `ls ./GTEX/02.sample/tissue`;
do
echo $i
cp ./GTEX/05.analysis_data/eigen/${i}.eigen > ./GTEX/05.analysis_data/combine/${i}.eigen.tmp
sed '1d' ./GTEX/05.analysis_data/peer/${i}.txt> ./GTEX/05.analysis_data/combine/${i}.peer.tmp
cat ./GTEX/05.analysis_data/combine/${i}.covariates.tmp ./GTEX/05.analysis_data/combine/${i}.eigen.tmp ./GTEX/05.analysis_data/combine/${i}.peer.tmp > ./GTEX/05.analysis_data/combine/${i}.covariates
done

for i in Testis Ovary Uterus Vagina Prostate
do
grep -v "V2" ./GTEX/05.analysis_data/combine/${i}.covariates > ./GTEX/05.analysis_data/combine/${i}_tmp.covariates
rm ./GTEX/05.analysis_data/combine/${i}.covariates
mv ./GTEX/05.analysis_data/combine/${i}_tmp.covariates ./GTEX/05.analysis_data/combine/${i}.covariates
done



for i in Brain-Putamen--basal-ganglia Brain-Spinal-cord--cervical-c-1 Brain-Substantia-nigra Cells-Cultured-fibroblasts Cells-EBV-transformed-lymphocytes Esophagus-Gastroesophageal-Junction Esophagus-Mucosa Esophagus-Muscularis Heart-Atrial-Appendage Heart-Left-Ventricle Kidney-Cortex Liver Minor-Salivary-Gland Muscle-Skeletal Nerve-Tibial Pituitary Small-Intestine-Terminal-Ileum Spleen
do
echo $i
echo "#! /bin/bash
Rscript ./GTEX/peer/15.R ${i}" > scripts/${i}.sh
chmod a+x scripts/$i.sh
nohup ./scripts/$i.sh &> scripts/$i.log&
  done




for i in Breast-Mammary-Tissue Colon-Sigmoid Colon-Transverse Pancreas Stomach Whole-Blood 
do
echo $i
echo "#! /bin/bash
Rscript ./GTEX/peer/30.R ${i}" > scripts/${i}.sh
chmod a+x scripts/$i.sh
nohup ./scripts/$i.sh &> scripts/$i.log&
  done


for i in Lung Skin-Not-Sun-Exposed--Suprapubic Skin-Sun-Exposed--Lower-leg Thyroid
do
echo $i
echo "#! /bin/bash
Rscript ./GTEX/peer/35.R ${i}" > scripts/${i}.sh
chmod a+x scripts/$i.sh
nohup ./scripts/$i.sh &> scripts/$i.log&
  done

for i in  Ovary Uterus Vagina Prostate
do
echo $i
echo "#! /bin/bash
Rscript ./GTEX/peer/gender.R ${i}" > scripts/${i}.sh
chmod a+x scripts/$i.sh
nohup ./scripts/$i.sh &> scripts/$i.log&
  done


for i in Testis
do
echo $i
echo "#! /bin/bash
Rscript ./GTEX/peer/gender_30.R ${i}" > scripts/${i}.sh
chmod a+x scripts/$i.sh
nohup ./scripts/$i.sh &> scripts/$i.log&
  done


for i in Adipose-Subcutaneous Adipose-Visceral--Omentum Adrenal-Gland Artery-Aorta Artery-Coronary Artery-Tibial Brain-Amygdala Brain-Anterior-cingulate-cortex--BA24 Brain-Caudate--basal-ganglia Brain-Cerebellar-Hemisphere Brain-Cerebellum Brain-Cortex Brain-Frontal-Cortex--BA9 Brain-Hippocampus Brain-Hypothalamus Brain-Nucleus-accumbens--basal-ganglia
do 
echo $i
echo "#! /bin/bash
Rscript ./GTEX/peer/15.R ${i}" > scripts/${i}.sh
chmod a+x scripts/$i.sh
nohup ./scripts/$i.sh &> scripts/$i.log&
  done