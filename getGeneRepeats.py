import argparse
import gff2repeats as repcode

parser = argparse.ArgumentParser()
parser.add_argument("--gff", help = "GFF file of the sample")
parser.add_argument("--repsfilePrefix", help = "repeats file Prefix file of the sample")
parser.add_argument("--gene", help = "Name of the gene")
args = parser.parse_args()

reps = ["boxA", "boxB", "boxC", "RUP", "SPRITE"]

for reptype in reps:
    repfile = ".".join([args.repsfilePrefix, reptype, "rep"])
    print(f"searching in {repfile}")
                
    repcode.parseRepOut(repfile)
# all chromosomes getting ouptut..to correct this error
    chrname, start, end, strand = repcode.getGeneLoc(
        gff_file=args.gff, genename=[args.gene])
    print(chrname, start, end, strand)
# print(len(geneIDS))
    genetest = repcode.gene(name=args.gene, start=start, end=end,
                strand=strand, chrName=chrname)
    nups, ndowns = genetest.extract_reps(repsFile=repfile)
    print(genetest.name, reptype, nups, ndowns)
