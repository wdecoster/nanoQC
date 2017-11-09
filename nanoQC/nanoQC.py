# wdecoster

import os
from argparse import ArgumentParser
import gzip
import logging
from Bio import SeqIO
from .version import __version__
import numpy as np
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt


def get_args():
    parser = ArgumentParser(
        description="Investigate nucleotide composition and base quality.")
    parser.add_argument("-v", "--version",
                        help="Print version and exit.",
                        action="version",
                        version='NanoQC {}'.format(__version__))
    parser.add_argument("--outdir",
                        help="Specify directory in which output has to be created.",
                        default=".")
    parser.add_argument("-f", "--format",
                        help="Specify the output format of the plots.",
                        default="png",
                        type=str,
                        choices=['eps', 'jpeg', 'jpg', 'pdf', 'pgf', 'png', 'ps',
                                 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff'])
    parser.add_argument("fastq",
                        help="Reads data in fastq format.")
    return parser.parse_args()


def main():
    args = get_args()
    logging.basicConfig(
        format='%(asctime)s %(message)s',
        filename=os.path.join(args.outdir, "NanoQC.log"),
        level=logging.INFO)
    logging.info("NanoQC started.")
    sizeRange = length_histogram(
        fqin=gzip.open(args.fastq, 'rt'),
        name=os.path.join(args.outdir, "SequenceLengthDistribution.png"))
    fq = get_bin(gzip.open(args.fastq, 'rt'), sizeRange)
    logging.info("Using {} reads for plotting".format(len(fq)))
    fqbin = [dat[0] for dat in fq]
    qualbin = [dat[1] for dat in fq]
    logging.info("Creating plots...")
    per_base_sequence_content_and_quality(fqbin, qualbin, args.outdir, args.format)
    logging.info("per base sequence content and quality completed.")
    logging.info("Finished!")


def per_base_sequence_content_and_quality(fqbin, qualbin, outdir, figformat):
    fig, axs = plt.subplots(2, 2, sharex='col', sharey='row')
    lines = plot_nucleotide_diversity(axs[0, 0], fqbin)
    plot_nucleotide_diversity(axs[0, 1], fqbin, invert=True)
    l_Q = plot_qual(axs[1, 0], qualbin)
    plot_qual(axs[1, 1], qualbin, invert=True)
    plt.setp([a.get_xticklabels() for a in axs[0, :]], visible=False)
    plt.setp([a.get_yticklabels() for a in axs[:, 1]], visible=False)
    for ax in axs[:, 1]:
        ax.set_ylabel('', visible=False)
    for ax in axs[0, :]:
        ax.set_xlabel('', visible=False)
    # Since axes are shared I should only invert once. Twice will restore the original axis order!
    axs[0, 1].invert_xaxis()
    plt.suptitle("Per base sequence content and quality")
    axl = fig.add_axes([0.4, 0.4, 0.2, 0.2])
    ax.plot()
    axl.axis('off')
    lines.append(l_Q)
    plt.legend(lines, ['A', 'T', 'G', 'C', 'Quality'], loc="center", ncol=5)
    plt.savefig(os.path.join(outdir, "PerBaseSequenceContentQuality." +
                             figformat), format=figformat, dpi=500)


def get_lengths(fastq):
    '''
    Loop over the fastq file, extract length of sequences
    '''
    return np.array([len(record) for record in SeqIO.parse(fastq, "fastq")])


def length_histogram(fqin, name):
    '''
    Create a histogram, and return the bin edges of the bin containing the most reads
    '''
    logging.info("Creating length histogram to find bin with most reads.")
    lengths = get_lengths(fqin)
    plt.hist(lengths, bins='auto')
    plt.savefig(name, format='png', dpi=100)
    plt.close("all")
    hist, bin_edges = np.histogram(lengths, bins='auto')
    maxindex = np.argmax(hist)
    return (bin_edges[maxindex], bin_edges[maxindex + 1])


def get_bin(fq, sizeRange):
    '''
    Loop over the fastq file
    Extract list of nucleotides and list of quality scores in tuples in list
    Only select those reads of which the length is within the size range
    '''
    logging.info("Extracting nucleotides and quality scores of selected bin.")
    return [(list(rec.seq), list(rec.letter_annotations["phred_quality"]))
            for rec in SeqIO.parse(fq, "fastq") if sizeRange[0] < len(rec) < sizeRange[1]]


def plot_nucleotide_diversity(ax, fqlists, invert=False):
    '''
    Create a FastQC-like "￼Per base sequence content" plot
    Plot fraction of nucleotides per position
    zip will stop when shortest read is exhausted
    '''
    if invert:
        fqlists = [list(reversed(read)) for read in fqlists]
    numreads = len(fqlists)
    sns.set_style("darkgrid")
    l_A, = ax.plot(
        np.array([pos.count('A') / numreads for pos in zip(*fqlists)]), 'green', label='A')
    l_T, = ax.plot(
        np.array([pos.count('T') / numreads for pos in zip(*fqlists)]), 'red', label='T')
    l_G, = ax.plot(
        np.array([pos.count('G') / numreads for pos in zip(*fqlists)]), 'black', label='G')
    l_C, = ax.plot(
        np.array([pos.count('C') / numreads for pos in zip(*fqlists)]), 'blue', label='C')
    if invert:
        ax.set_xticklabels(-1 * ax.get_xticks().astype(int))
    return [l_A, l_T, l_G, l_C]


def plot_qual(ax, quallist, invert=False):
    '''
    Create a FastQC-like "￼Per base sequence quality￼" plot
    Plot average quality per position
    zip will stop when shortest read is exhausted
    '''
    sns.set_style("darkgrid")
    if invert:
        l_Q, = ax.plot(np.array([np.mean(position) for position in zip(
            *[list(reversed(read)) for read in quallist])]), 'orange', label="Quality")
        ax.set_xlabel('Position in read from end')
        ax.set_xticklabels(-1 * ax.get_xticks().astype(int))
    else:
        l_Q, = ax.plot(np.array([np.mean(position)
                                 for position in zip(*quallist)]), 'orange', label="Quality")
        ax.set_xlabel('Position in read from start')
    return l_Q


if __name__ == "__main__":
    main()
