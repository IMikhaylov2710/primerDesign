from helpers.classes import Variant, Primer, ProbeBatch

newVar = Variant(155183132, 1, 'TRIM46', delta=150)
newVar.getRegion('hg38')
newVar.findPrimers()
newVar.getPairs()
newVar.printPairs()
newVar.reducePairs()
newVar.writeToFasta('/home/ivan/Downloads/trim46.fasta')
newVar.runBlast('/home/ivan/Downloads/trim46.tsv')