from gff2repeats import *

parseRepOut("SRR8879299.fasta.boxA.rep")
# all chromosomes getting ouptut..to correct this error
chrname, start, end, strand = getGeneLoc(
    gff_file="1453.gff", genename=["psaA"])
print(chrname, start, end, strand)
# print(len(geneIDS))
genetest = gene(name="psaA", start=start, end=end,
                strand=strand, chrName=chrname)
nups, ndowns = genetest.extract_reps(repsFile="SRR8879299.fasta.boxA.rep")
print(genetest.name, "boxA", nups, ndowns)
nups, ndowns = genetest.extract_reps(repsFile="SRR8879299.fasta.boxB.rep")
print(genetest.name, "boxB", nups, ndowns)
nups, ndowns = genetest.extract_reps(repsFile="SRR8879299.fasta.boxC.rep")
print(genetest.name, "boxC", nups, ndowns)
nups, ndowns = genetest.extract_reps(repsFile="SRR8879299.fasta.RUP.rep")
print(genetest.name, "RUP", nups, ndowns)
nups, ndowns = genetest.extract_reps(repsFile="SRR8879299.fasta.SPRITE.rep")
print(genetest.name, "SPRITE", nups, ndowns)
