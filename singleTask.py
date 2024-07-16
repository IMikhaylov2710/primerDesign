from helpers.classes import Variant, Primer, ProbeBatch
import argparse

parser = argparse.ArgumentParser(description='Script for primer design')
parser.add_argument("-chr", "--Chromosome", help = "select chromosome number for rs",
                   nargs = '?',
                   type = int)
parser.add_argument("-coord", "--Coordinate", help = "select coordinate number for rs",
                   nargs = '?',
                   type = int)
parser.add_argument("-gn", "--GeneName", help = "select gene name for design",
                    nargs = '?', 
                    type = str, 
                    default = "-")
parser.add_argument("-fasta", "--FastaPath", help = "select path where to write .fasta with primers",
                    nargs = '?', 
                    type = str)
parser.add_argument("-tsv", "--TsvPath", help = "select path where to write .tsv after primer blast",
                    nargs = '?', 
                    type = str)
args = parser.parse_args()

newVar = Variant(args.Coordinate, args.Chromosome, args.GeneName, delta=150)
newVar.getRegion('hg38')
newVar.findPrimers()
newVar.getPairs()
newVar.printPairs()
newVar.reducePairs()
newVar.writeToFasta(args.FastaPath)
newVar.runBlast(args.TsvPath)