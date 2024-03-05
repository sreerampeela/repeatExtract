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
            genechr, start, end, strand = -1, -1, -1, -1
        return genechr, start, end, strand


def parseRepOut(repoutfile):
    """To read a repeats output file and create a dataframe from the output"""
    repsh = open(repoutfile, 'r', encoding="utf-8", newline="\n")
    repsOut = repsh.readlines()[13:]
    pattern = re.compile(r'([\d.]+)\s+\(bits\) f:(\d+) t:(\d+) Target:\s+(.+)')
    data = [match.groups() for match in pattern.finditer("\n".join(repsOut))]
    df = pd.DataFrame(data, columns=['Bits', 'Start', 'Stop', 'Target'])
    # print(df.head())
    return df


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
        repsDF = parseRepOut(repsFile)
        chrmatch = "_".join(chrid.split("_")[:5])
        x = repsDF["Target"].str.startswith(chrmatch)
        # print(x)
        repsDF['geneloc'] = x
        matchedDF = repsDF[repsDF["geneloc"] == True]
        # print(matchedDF)
        startpos, endpos = self.start, self.end
        nups = matchedDF[matchedDF["Start"] < startpos].shape[0]
        ndowns = matchedDF[matchedDF["Start"] > startpos].shape[0]

        # repsDF.to_csv("tmp_boxA.csv", index=False)
        return nups, ndowns
