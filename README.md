# repeatExtract
The scripts provided aid in extracting different types of repeats 
located both upstream and downstream of a gene.

Dependencies:
Python - >v3.7

Biopython

Pandas

Numpy

For both the scripts, the help menu can be displayed using "--help" argument 
when running the script


Using GFF files:
Get number of repeats for each gene for each sample
The GFF file is parsed initially.
From the GFF file, the location of gene of interest is determined, and various
repeats present upstream and downstream of this gene is detected and reported

Usage:
python3 getGeneRepeats.py --gff 1453.gff --repsfile SRR8879299.fasta.boxA --gene psaA

Using ABRICATE output:
To use output of ABRICATE and SPN_REPEATS, run abricate2repeatsParser.py
The script accepts the directory path and output file prefix.
To run the script, both ABRICATE output and repeats output file
for each sample should be in the same directory.
BoxA, BoxB, BoxC, RUP and SPRITE Repeats present in all genes listed in the 
ABRICATE tsv file for a sample will be written as separate output files.


usage:
python3 abricate2repsParser.py --dir $PWD --out test

