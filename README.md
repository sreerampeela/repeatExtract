# repeatExtract
Get number of repeats for each gene for each sample
The GFF file is parsed initially.
From the GFF file, the location of gene of interest is determined, and various
repeats present upstream and downstream of this gene is detected and reported

Sample usage


python3 getGeneRepeats.py --gff 1453.gff --repsfile SRR8879299.fasta.boxA --gene psaA
