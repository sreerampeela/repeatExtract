import pandas as pd
import numpy as np
import os
import re


class Gene:
    def __init__(self, name, start, end, strand, chrName):
        self.name = name
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
        repsDF['geneloc'] = x
        matchedDF = repsDF[repsDF["geneloc"] == True]
        if self.strand == "+":
            startpos, endpos = self.start, self.end
        elif self.strand == "-":
            startpos, endpos = self.end, self.start
        RepstartPos = [int(i) for i in list(matchedDF["Start"])]
        RependPos = [int(i) for i in list(matchedDF["Stop"])]
        nups = 0
        ndowns = 0
        for i in RepstartPos:
            if i < startpos:
                nups += 1
        for j in RependPos:
            if j > endpos:
                ndowns += 1
        return nups, ndowns


def parseRepOut(repoutfile):
    """To read a repeats output file and create a dataframe from the output"""
    repsh = open(repoutfile, 'r', encoding="utf-8", newline="\n")
    repsOut = repsh.readlines()[13:]
    pattern = re.compile(r'([\d.]+)\s+\(bits\) f:(\d+) t:(\d+) Target:\s+(.+)')
    data = [match.groups() for match in pattern.finditer("\n".join(repsOut))]
    df = pd.DataFrame(data, columns=['Bits', 'Start', 'Stop', 'Target'])
    # print(df.head())
    return df


def parseAbricateFile(vfsfile="spn_jip.tab", sampleFasta="GPS_IN_CRL_SPN_259.fasta"):
    vfsDF = pd.read_csv(vfsfile, header=0, sep="\t")
    sampleVFS = vfsDF[vfsDF["#FILE"] == sampleFasta]
    vfgenes = list(sampleVFS["GENE"])
    # print(vfgenes)
    vfgeneOBJS = []
    for vfgene in vfgenes:
        vfgene_chr = str(sampleVFS[sampleVFS["GENE"]
                         == vfgene]["SEQUENCE"].values[0])
        vfgene_start = int(
            sampleVFS[sampleVFS["GENE"] == vfgene]["START"].values[0])
        vfgene_end = int(sampleVFS[sampleVFS["GENE"]
                         == vfgene]["END"].values[0])
        vfgene_strand = str(
            sampleVFS[sampleVFS["GENE"] == vfgene]["STRAND"].values[0])
        vfgeneOBJ = Gene(name=vfgene, start=int(vfgene_start), end=int(
            vfgene_end), strand=vfgene_strand, chrName=vfgene_chr)
        vfgeneOBJS.append(vfgeneOBJ)
    return vfgeneOBJS


fastafiles = [i for i in os.listdir() if i.endswith(".fasta") is True]
sampleIDS = []
geneNames = []
upstreams = []
downstreams = []

for sampleID in fastafiles:
    vfs = parseAbricateFile(sampleFasta=sampleID)
    for vf in vfs:
        repFileName = ".".join([sampleID, "SPRITE", "rep"])
        nups, ndowns = vf.extract_reps(repsFile=repFileName)
        # print(vf.name, vf.end, vf.chrName, vf.start, vf.end, nups, ndowns)
        sampleIDS.append(sampleID)
        geneNames.append(vf.name)
        upstreams.append(nups)
        downstreams.append(ndowns)


dfRes = pd.DataFrame()
dfRes["Sample"] = sampleIDS
dfRes["Gene"] = geneNames
dfRes["Upstream"] = upstreams
dfRes["Downstream"] = downstreams
print(dfRes.head())

dfRes.to_csv("SPRITE_reps_vfs.csv", index=False)
