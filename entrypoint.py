from helpers.classes import Variant, Primer, ProbeBatch
from helpers.dbLogic import connectToDBSNP, getByRs, getRsByCoord
from helpers.slurmHelpers import addSlurmTask
import argparse
import os

parser = argparse.ArgumentParser(description='Script for primer design')
parser.add_argument("-b", "--BatchRS", help = "list of rs, separated by ',' to use for batch analysis")
parser.add_argument("--test", action = "store_true", help = "use this flag to only run test case")
parser.add_argument("--batch", action = "store_true", help = "use this flag to perform batch primer analysis")
parser.add_argument("--slurm", action = "store_true", help = "use multiprocessing using SLURM cluster")
args = parser.parse_args()

jobID = 1
resultingRsList = [str(rs.strip())for rs in args.BatchRS.split(',')]
if not args.test:
    for rs in resultingRsList:
        addSlurmTask(f'job{jobID}', 
                     f'job{jobID}', 
                     'chr1', 
                     '1000000', 
                     'testGene', 
                     '~/out/fastaPath/', 
                     '~/out/tsvPath')
        jobID+=1
else:
    rsList = getByRs('rs1274', conn)
    danger = getRsByCoord('chr1', '100000', '200000', conn)