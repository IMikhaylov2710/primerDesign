from helpers.classes import Variant, Primer, ProbeBatch
from helpers.DBlogic import connectToDBSNP, getByRs, getRsByCoord
from helpers.slurmHelpers import addSlurmTask
import argparse
import os

parser = argparse.ArgumentParser(description='Script for primer design')
parser.add_argument("-b", "--BatchRS", help = "list of rs, separated by ',' to use for batch analysis")
parser.add_argument("--test", action = "store_true", help = "use this flag to only run test case")
parser.add_argument("--batch", action = "store_true", help = "use this flag to perform batch primer analysis")
parser.add_argument("--slurm", action = "store_true", help = "use multiprocessing using SLURM cluster")
args = parser.parse_args()

if not args.test:
    addSlurmTask()
else:
    conn = connectToDBSNP()
    rsList = getByRs('rs1274', conn)
    danger = getRsByCoord('chr1', '100000', '200000', conn)