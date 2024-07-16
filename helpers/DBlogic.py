import sqlalchemy
import pandas as pd
from sqlalchemy import insert

def connectToDBSNP(address):
    dbEngine=sqlalchemy.create_engine(address)
    return dbEngine

def getByRs(rs, engine):
    df = pd.read_sql(f'SELECT * FROM dbsnp WHERE rs = {rs}', engine)
    return df

def getRsByCoord(chromosome, start, end, engine):
    df = pd.read_sql(f'SELECT * FROM dbsnp WHERE chromosome = {chromosome} AND coordinate >= {start} AND coordinate <= {end}', engine)
    return df

def createTask(engine, table, taskName, rsList):
    with engine.connect() as conn:
        result = conn.execute(
            insert(table), 
                [
                    {"taskName" : taskName, 
                     "rsList" : rsList}
                ],
                )
        conn.commit()
    return True

def runTask(engine, table, taskName, newStatus):
    with engine.connect() as conn:
        conn = engine.connect()
        stmt = table.update().values(status=newStatus).where(table.taskName == taskName)
        conn.execute(stmt)
    return True