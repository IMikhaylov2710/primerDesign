import os

#installing prerequesit programms
os.system('conda install -c bioconda ncbi-blast vcftools')

#wget hg38.fasta
os.system('cd REFs/ && wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz')

#gunzipping and making blast DB 
os.system('cd REFs/ && gunzip hg38.fa && makeblastdb -dbtype nucl -in hg38.fa')

#wget dbsnp most common alleles for allele dropout info
os.system('cd DB/ && wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz && gunzip 00-common_all.vcf.gz')

