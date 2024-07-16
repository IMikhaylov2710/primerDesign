import sqlalchemy
import pandas as pd

def connectToDBSNP(address):
    dbEngine=sqlalchemy.create_engine(address)
    return dbEngine

def getByRs(rs, engine):
    df = pd.read_sql(f'SELECT * FROM dbsnp WHERE rs = {rs}', engine)
    return df

def gerRsByCoord(chromosome, start, end, engine):
    df = pd.read_sql(f'SELECT * FROM dbsnp WHERE chromosome = {chromosome} AND coordinate >= {start} AND coordinate <= {end}', engine)
    return df