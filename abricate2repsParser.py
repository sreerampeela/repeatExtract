import os
import argparse
import abricate2reps as ab2rep
import pandas as pd
import numpy as np

parse = argparse.ArgumentParser()
parse.add_argument(
    "--dir", help="directory with abricate output files and repeat files")
parse.add_argument("--out", help="output file prefix")
# parse.add_argument("--dir", help="directory with abricate output files and repeat files")

args = parse.parse_args()
fastafiles = [i for i in os.listdir(args.dir) if i.endswith(".fasta") is True]


for reptype in ["boxA", "boxB", "boxC", "RUP", "SPRITE"]:
    sampleIDS = []
    geneNames = []
    upstreams = []
    downstreams = []
    outfile = "_".join([args.out, reptype]) + ".csv"
    for sampleID in fastafiles:
        vfs = ab2rep.parseAbricateFile(sampleFasta=sampleID)
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

    dfRes.to_csv(outfile, index=False)
