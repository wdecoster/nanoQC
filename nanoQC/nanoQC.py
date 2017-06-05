# wdecoster

import os
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import argparse
import gzip
from Bio import SeqIO

__version__ = 0.1.0


def getArgs():
	parser = argparse.ArgumentParser(
					description="Investigate overrepresented nucleotides in reads of median length.")
	parser.add_argument("-v", "--version",
					help="Print version and exit.",
					action="version",
					version='NanoQC {}'.format(__version__))
	parser.add_argument("--outdir",
					help="Specify directory in which output has to be created.",
					default=".")
	parser.add_argument("fastq",
					help="Reads data in fastq format.")
	return parser.parse_args()


def main():
	sizeRange = LengthHistogram(gzip.open(args.fastq, 'rt'))
	fq = getBin(gzip.open(args.fastq, 'rt'), sizeRange)
	print("Using {} reads for plotting".format(len(fq)))
	plotNucleotideDiversity([dat[0] for dat in fq])
	plotQual([dat[1] for dat in fq])

def getLengths(fastq):
	'''	Loop over the fastq file, extract length of sequences'''
	return np.array([len(record) for record in SeqIO.parse(fastq, "fastq")])


def LengthHistogram(fqin):
	''' Create a histogram, and return the number of reads per bin and bin edges'''
	lengths = getLengths(fqin)
	plt.hist(lengths, bins='auto')
	plt.savefig(os.path.join(args.outdir, "SequenceLengthDistribution.png"), format='png', dpi=1000)
	plt.close("all")
	hist, bin_edges = np.histogram(lengths, bins='auto')
	maxindex = np.argmax(hist)
	return((bin_edges[maxindex], bin_edges[maxindex + 1]))


def getBin(fq, sizeRange):
	'''
	Loop over the fastq file
	Extract list of nucleotides and list of quality scores in tuples in list
	Only select those reads of which the length is within the size range
	'''
	return [(list(rec.seq), list(rec.letter_annotations["phred_quality"])) for rec in SeqIO.parse(fq, "fastq") if sizeRange[0] < len(rec) < sizeRange[1]]


def plotNucleotideDiversity(fqlists):
	'''
	Create a FastQC-like "￼Per base sequence content" plot
	Plot fraction of nucleotides per position
	zip will stop when shortest read is exhausted
	'''
	numreads = len(fqlists)
	sns.set_style("darkgrid")
	plt.plot(np.array([position.count('A') / numreads for position in zip(*fqlists)]), 'green')
	plt.plot(np.array([position.count('T') / numreads for position in zip(*fqlists)]), 'red')
	plt.plot(np.array([position.count('G') / numreads for position in zip(*fqlists)]), 'black')
	plt.plot(np.array([position.count('C') / numreads for position in zip(*fqlists)]), 'blue')
	plt.legend(['A', 'T', 'G', 'C'], loc='upper right')
	plt.xlabel('Position in read')
	plt.ylabel('Per base sequence content')
	plt.savefig(os.path.join(args.outdir, "PerBaseSequenceContent.png"), format='png', dpi=1000)
	plt.close("all")


def plotQual(quallist):
	'''
	Create a FastQC-like "￼Per base sequence quality￼" plot
	Plot average quality per position
	zip will stop when shortest read is exhausted
	'''
	sns.set_style("darkgrid")
	plt.plot(np.array([np.mean(position) for position in zip(*quallist)]), 'red')
	plt.xlabel('Position in read')
	plt.ylabel('Per base sequence quality')
	plt.savefig(os.path.join(args.outdir, "PerBaseSequenceQuality.png"), format='png', dpi=1000)
	plt.close("all")

if __name__ == "__main__":
	main()
