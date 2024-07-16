import os
import pandas as pd
import sqlite3
from sqlalchemy import create_engine
from sqlalchemy import Column, DateTime, Integer, String, Boolean, ForeignKey
from sqlalchemy.sql import func
from sqlalchemy.ext.declarative import declarative_base

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

class ProgressStatus(Base):

    __tablename__ = "status"

    id = Column(Integer, primary_key=True)
    taskName = Column(String)
    status = Column(String, default = 'Started')
    ready = Column(Boolean, default = False) 
    timeCreated = Column(DateTime(timezone = True), server_default=func.now())
    timeUpdated = Column(DateTime(timezone = True), onupdate=func.now())

    def __init__(self, name):
        
        self.name = name

class RsStatus(Base):

    __tableName__ = "rsstatus"

    id = Column(Integer, primary_key=True)
    taskId = Column(ForeignKey("status.id"))
    rs = Column(String)
    finished = Column(Boolean, default = False)
    timeCreated = Column(DateTime(timezone = True), server_default=func.now())
    timeUpdated = Column(DateTime(timezone = True), onupdate=func.now())

    def __init__(self, name):
        
        self.name = name

Base.metadata.create_all(engine)

#populating database
print('===POPULATING DATABASE===')
e = 1
counter = 0 
rows = []
batchCounter = 0
with sqlite3.connect('DB/dbsnp.db') as conn:
    with open('DB/00-common_all.vcf', 'r') as handle:
        for lin in handle:
            if not lin.startswith('#'):
                ls = lin.split('\t')
                info = ls[-1].split(';')
                for unit in info:
                    if unit.startswith('CAF'):
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
                            alt = unit.split('=')[1].split(',')[1]
                            if alt == '.':
                                alt = 0
                            else:
                                alt = float(alt)
                row = [e, 'chr'+ls[0], ls[1], ls[2], ls[3], ls[4], ref, alt]
                e+=1
                if counter == 10000:
                    counter = 0
                    batchCounter +=1
                    print('added 10000 rows, batch number is', batchCounter)
                    rowsDf = pd.DataFrame(rows, columns = ['id', 'chromosome', 'coordinate', 'rs', 'refAllele', 'altAllele', 'refAlleleFrequency', 'altAlleleFrequencySum'])
                    rowsDf.to_sql(name='dbsnp', con = conn, if_exists = 'append', index = False)
                    rows = []
                else:
                    rows.append(row)
                    counter+=1
