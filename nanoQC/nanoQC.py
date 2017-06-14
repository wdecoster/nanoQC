# wdecoster

import os
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import argparse
import gzip
from Bio import SeqIO
from .version import __version__



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
	args = getArgs()
	sizeRange = LengthHistogram(gzip.open(args.fastq, 'rt'), os.path.join(args.outdir, "SequenceLengthDistribution.png"))
	fq = getBin(gzip.open(args.fastq, 'rt'), sizeRange)
	print("Using {} reads for plotting".format(len(fq)))
	fqbin = [dat[0] for dat in fq]
	qualbin = [dat[1] for dat in fq]
	perBaseSequenceContentQuality(fqbin, qualbin, args.outdir)

def perBaseSequenceContentQuality(fqbin, qualbin, outdir):
	fig, axs = plt.subplots(2,2, sharex='col', sharey='row')
	lines = plotNucleotideDiversity(axs[0,0], fqbin, os.path.join(outdir, "ForwardPerBaseSequenceContent.png"))
	plotNucleotideDiversity(axs[0,1], fqbin, os.path.join(outdir, "ReversePerBaseSequenceContent.png"), invert=True)
	l_Q = plotQual(axs[1,0], qualbin, os.path.join(outdir, "ForwardPerBaseSequenceQuality.png"))
	plotQual(axs[1,1], qualbin, os.path.join(outdir, "ReversePerBaseSequenceQuality.png"), invert=True)
	plt.setp([a.get_xticklabels() for a in axs[0, :]], visible=False)
	plt.setp([a.get_yticklabels() for a in axs[:, 1]], visible=False)
	for ax in axs[:,1]:
		ax.set_ylabel('', visible=False)
	for ax in axs[0, :]:
		ax.set_xlabel('', visible=False)
	axs[0,1].invert_xaxis()  # Since axes are shared I should only invert once. Twice will restore the original axis order!
	plt.suptitle("Per base sequence content and quality")
	axl = fig.add_axes([0.4,0.4,0.2,0.2])
	pie = ax.plot()
	axl.axis('off')
	lines.append(l_Q)
	plt.legend(lines, ['A', 'T', 'G', 'C', 'Quality'], loc="center", ncol=5)
	plt.savefig(os.path.join(outdir, "PerBaseSequenceContentQuality.png"), format='png', dpi=500)



def getLengths(fastq):
	'''
	Loop over the fastq file, extract length of sequences
	'''
	return np.array([len(record) for record in SeqIO.parse(fastq, "fastq")])


def LengthHistogram(fqin, name):
	'''
	Create a histogram, and return the bin edges of the bin containing the most reads
	'''
	lengths = getLengths(fqin)
	plt.hist(lengths, bins='auto')
	plt.savefig(name, format='png', dpi=100)
	plt.close("all")
	hist, bin_edges = np.histogram(lengths, bins='auto')
	maxindex = np.argmax(hist)
	return (bin_edges[maxindex], bin_edges[maxindex + 1])


def getBin(fq, sizeRange):
	'''
	Loop over the fastq file
	Extract list of nucleotides and list of quality scores in tuples in list
	Only select those reads of which the length is within the size range
	'''
	return [(list(rec.seq), list(rec.letter_annotations["phred_quality"])) for rec in SeqIO.parse(fq, "fastq") if sizeRange[0] < len(rec) < sizeRange[1]]


def plotNucleotideDiversity(ax, fqlists, name, invert=False):
	'''
	Create a FastQC-like "￼Per base sequence content" plot
	Plot fraction of nucleotides per position
	zip will stop when shortest read is exhausted
	'''
	if invert:
		fqlists = [list(reversed(read)) for read in fqlists]
	numreads = len(fqlists)
	sns.set_style("darkgrid")
	l_A, = ax.plot(np.array([position.count('A') / numreads for position in zip(*fqlists)]), 'green', label='A')
	l_T, = ax.plot(np.array([position.count('T') / numreads for position in zip(*fqlists)]), 'red', label='T')
	l_G, = ax.plot(np.array([position.count('G') / numreads for position in zip(*fqlists)]), 'black', label='G')
	l_C, = ax.plot(np.array([position.count('C') / numreads for position in zip(*fqlists)]), 'blue', label='C')
	if invert:
		ax.set_xticklabels(-1*ax.get_xticks().astype(int))
	return [l_A, l_T, l_G, l_C]

def plotQual(ax, quallist, name, invert=False):
	'''
	Create a FastQC-like "￼Per base sequence quality￼" plot
	Plot average quality per position
	zip will stop when shortest read is exhausted
	'''
	sns.set_style("darkgrid")
	if invert:
		l_Q, = ax.plot(np.array([np.mean(position) for position in zip(*[list(reversed(read)) for read in quallist])]), 'orange', label="Quality")
		ax.set_xlabel('Position in read from end')
		ax.set_xticklabels(-1*ax.get_xticks().astype(int))
	else:
		l_Q, = ax.plot(np.array([np.mean(position) for position in zip(*quallist)]), 'orange', label="Quality")
		ax.set_xlabel('Position in read from start')
	return l_Q


if __name__ == "__main__":
	main()
