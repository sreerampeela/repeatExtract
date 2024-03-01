from Bio import SeqIO
import pandas as pd
import re


def readGFF(gffFile):
    with open(gffFile, 'r', encoding="utf-8", newline="\n") as gffInput:
        data = [i.split("\t") for i in gffInput.readlines()
                if i.startswith("#") is False]
        df = pd.DataFrame(data, columns=[
                          "Chr", "method", "type", "start", "end", "d", "strand", 
                          "p", "info"])
        df = df.dropna()
        # print(df.head())
        return df


def get_gene_names(gff_file):
    """From gff files, extract gene names. In case of hypothetical proteins or 
    proteins with no name, the RefSeq protein id or product name is extracted 
    as gene name"""
    genenames = []
    df = readGFF(gffFile=gff_file)
    datacol = list(df["info"])
    # for i in range(df.shape[0]):
    for j in datacol:
        # print(j)
        geneid = [k for k in j.split(";") if k.startswith("gene")]

        if len(geneid) == 0:
            geneid = [k for k in j.split(";") if k.startswith("product")]
        # geneName = geneid[0].split("=")[1]
        # else:
        if len(geneid) == 0:
            geneid = [k for k in j.split(";") if k.startswith("protein_id")]
        # print(geneid)
        geneName = geneid[0].split("=")[1]

        # print(geneName)
        genenames.append(geneName)
    return genenames

def getGeneLoc(gff_file, genename=["hypothetical protein"]):
    gff_file_df = readGFF(gffFile=gff_file)
    geneStarts = gff_file_df["start"]
    geneends = gff_file_df["end"]
    genestrand = gff_file_df["strand"]
    chrnames = gff_file_df["Chr"]
    geneids = get_gene_names(gff_file=gff_file)
    # print(geneids)
    for genesearch in genename:
        if genesearch in geneids:
            print(f"{genesearch} located")
            
            indval = geneids.index(genesearch)
            genechr, start, end, strand = chrnames[indval], geneStarts[indval], geneends[indval], genestrand[indval]
        else:
            print(f"{genesearch} not located")
            genechr, start, end, strand = -1,-1,-1, -1
        return genechr, start,end,strand



class gene:
    def __init__(self, name, start, end, strand, chrName):
        self.name = name
        # self.sequence = sequence
        self.strand = strand
        self.start = start
        self.end = end
        self.chrName = chrName
    
    def extract_reps(self, repsFile):
      nups, ndowns = 0, 0
      chrid = self.chrName
      repsh = open(repsFile, 'r', encoding="utf-8", newline="\n")
      repsOut = repsh.readlines()[13:]
      # print(repsOut)
      for k in repsOut:
        j = k.rstrip().strip().split(" ")
        for i in range(len(j)):
          x = re.search(j[i], chrid)
          if x:
            # print(k.split(" "))
            kdata = k.split(" ")
            startpos = [i.replace(" ", "") for i in kdata if i.startswith("f:") is True]
            # print(startpos)
            startposdata = [i.split(":") for i in startpos]
            # print(startposdata)
            if startposdata[0][1] == "":
              indval = kdata.index("f:")+1
              while kdata[indval] == "":
                indval += 1
              repstart = int(kdata[indval])
            else:
              repstart = int(startposdata[0][1])
            # repstart = int("".join([i.split(":") for i in kdata if i.startswith("f") is True]))
            # # repstop = [i.split(":") for i in kdata if i.startswith("t") is True]
            # print(repstart)
            if int(self.start) > repstart:
              # print("rep located upstream")
              nups += 1
            elif int(self.start) < repstart:
              # print("rep located downstream")
              ndowns += 1
            else:
              print("rep at gene start")
      return nups, ndowns
          


chrname, start,end, strand = getGeneLoc(gff_file="1453.gff", genename=["psaA"])
print(chrname, start, end, strand)
# print(len(geneIDS))
genetest = gene(name="psaA", start=start, end=end, strand=strand, chrName=chrname)
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


