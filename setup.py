import os
import sqlite3
import sqlalchemy
from sqlalchemy import create_engine, ForeignKey
from sqlalchemy import Column, Date, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from helpers.setupHelpers import add_rs

#installing prerequesit programms
os.system('conda install -c bioconda blast vcftools')
os.system('conda install biopython primer3-py colorama pandas Bio')

#wget hg38.fasta
if not os.path.isfile('REFs/hg38.fa'):
    os.system('cd REFs/ && wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz')

#gunzipping and making blast DB 
if not os.path.isfile('REFs/hg38.fa.nin'):
    os.system('cd REFs/ && gunzip hg38.fa && makeblastdb -dbtype nucl -in hg38.fa')

#wget dbsnp most common alleles for allele dropout info
if not os.path.isfile('DB/00-common_all.vcf'):
    os.system('cd DB/ && wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz && gunzip 00-common_all.vcf.gz')

#creating initial database
print('===CREATING DATABASE===')
engine = create_engine('sqlite:///DB/dbsnp.db', echo=True)
Base = declarative_base()

class VariantInfo(Base):

    __tablename__ = "dbsnp"

    id = Column(Integer, primary_key=True)
    chromosome = Column(String)  
    coordinate = Column(Integer)
    rs = Column(String)
    refAllele = Column(String)
    altAllele = Column(String)
    refAlleleFrequency = Column(String)
    altAlleleFrequencySum = Column(String)


    def __init__(self, name):

        self.name = name    

Base.metadata.create_all(engine)

#populating database
print('===POPULATING DATABASE===')
e = 1
with sqlite3.connect('DB/dbsnp.db') as conn:
    with open('DB/00-common_all.vcf', 'r') as handle:
        for lin in handle:
            if not lin.startswith('#'):
                ls = lin.split('\t')
                info = ls[-1].split(';')
                for unit in info:
                    if 'CAF' in unit:
                        if len(unit.split(',')) > 2:
                            ref = float(unit.split('=')[1].split(',')[0])
                            maf = [a for a in unit.split('=')[1].split(',')[1:]]
                            totalMaf = 0
                            for i in maf:
                                if i == '.':
                                    totalMaf+=0
                                else:
                                    totalMaf+=float(i)

                        else:
                            ref = float(unit.split('=')[1].split(',')[0])
                            alt = float(unit.split('=')[1].split(',')[1])
                row = [e, 'chr'+ls[0], ls[1], ls[2], ls[3], ls[4], ref, alt]
                e+=1
                add_rs(conn, row)