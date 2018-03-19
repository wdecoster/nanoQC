# wdecoster

import os
import sys
from argparse import ArgumentParser
import gzip
import logging
from Bio import SeqIO
from .version import __version__
import numpy as np
from bokeh.plotting import figure, save, output_file
from bokeh.layouts import gridplot
from bokeh.models import Range1d


def get_args():
    parser = ArgumentParser(
        description="Investigate nucleotide composition and base quality.")
    parser.add_argument("-v", "--version",
                        help="Print version and exit.",
                        action="version",
                        version='NanoQC {}'.format(__version__))
    parser.add_argument("-o", "--outdir",
                        help="Specify directory in which output has to be created.",
                        default=".")
    parser.add_argument("fastq",
                        help="Reads data in fastq.gz format.")
    return parser.parse_args()


def main():
    args = get_args()
    make_output_dir(args.outdir)
    logging.basicConfig(
        format='%(asctime)s %(message)s',
        filename=os.path.join(args.outdir, "NanoQC.log"),
        level=logging.INFO)
    logging.info("NanoQC started.")
    hist, sizeRange = length_histogram(
        fqin=gzip.open(args.fastq, 'rt'),
        name=os.path.join(args.outdir, "SequenceLengthDistribution.png"))
    fq = get_bin(gzip.open(args.fastq, 'rt'), sizeRange)
    logging.info("Using {} reads for plotting".format(len(fq)))
    logging.info("Creating plots...")
    seq_plots, qual_plots = per_base_sequence_content_and_quality(
        fqbin=[dat[0] for dat in fq],
        qualbin=[dat[1] for dat in fq],
        outdir=args.outdir)
    output_file("nanoQC.html", title="nanoQC_report")
    save(
        gridplot(children=[[hist], seq_plots, qual_plots],
                 plot_width=400,
                 plot_height=400,
                 sizing_mode="stretch_both")
    )
    logging.info("Finished!")


def make_output_dir(path):
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    except IOError:
        sys.exit("ERROR: No writing permission to the output directory.")


def per_base_sequence_content_and_quality(fqbin, qualbin, outdir):
    seq_plot_left = plot_nucleotide_diversity_bokeh(fqbin)
    seq_plot_right = plot_nucleotide_diversity_bokeh(
        fqbin, invert=True, y_range=seq_plot_left.y_range)
    qual_plot_left = plot_qual_bokeh(qualbin)
    qual_plot_right = plot_qual_bokeh(qualbin, invert=True, y_range=qual_plot_left.y_range)
    logging.info("Per base sequence content and quality completed.")
    return [seq_plot_left, seq_plot_right], [qual_plot_left, qual_plot_right]


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
    hist, edges = np.histogram(lengths, bins='auto')
    maxindex = np.argmax(hist)

    hist_norm = figure()
    hist_norm.quad(
        top=hist,
        bottom=0,
        left=edges[:-1],
        right=edges[1:],
        fill_color="#036564",
        line_color="#033649")
    hist_norm.xaxis[0].formatter.use_scientific = False

    return (hist_norm, (edges[maxindex], edges[maxindex + 1]))


def get_bin(fq, sizeRange):
    '''
    Loop over the fastq file
    Extract list of nucleotides and list of quality scores in tuples in list
    Only select those reads of which the length is within the size range
    '''
    logging.info("Extracting nucleotides and quality scores of selected bin.")
    return [(list(rec.seq), list(rec.letter_annotations["phred_quality"]))
            for rec in SeqIO.parse(fq, "fastq") if sizeRange[0] < len(rec) < sizeRange[1]]


def plot_nucleotide_diversity_bokeh(seqs, invert=False, y_range=None):
    x_length = len([pos for pos in zip(*seqs)])
    if invert:
        seqs = [list(reversed(read)) for read in seqs]
        p = figure(x_range=Range1d(start=x_length, end=0), y_range=y_range)
    else:
        p = figure()
    p.grid.grid_line_alpha = 0.3
    numreads = len(seqs)
    for nucl, color in zip(['A', 'T', 'G', 'C'], ['green', 'red', 'black', 'blue']):
        p.line(
            x=range(x_length),
            y=np.array([pos.count(nucl) / numreads for pos in zip(*seqs)]),
            color=color,
            legend=nucl)
    if invert:
        p.xaxis.axis_label = 'Position in read from end'
    else:
        p.xaxis.axis_label = 'Position in read from start'
    p.yaxis.axis_label = 'Frequency of nucleotide in read'
    return p


def plot_qual_bokeh(quallist, invert=False, y_range=None):
    '''
    Create a FastQC-like "￼Per base sequence quality￼" plot
    Plot average quality per position
    zip will stop when shortest read is exhausted
    '''
    x_length = len([pos for pos in zip(*quallist)])
    if invert:
        p = figure(x_range=Range1d(start=x_length, end=0), y_range=y_range)
        p.grid.grid_line_alpha = 0.3
        p.line(
            x=range(x_length),
            y=np.array([np.mean(position) for position in zip(
                *[list(reversed(read)) for read in quallist])]),
            color='orange',
            legend="Quality")
        p.xaxis.axis_label = 'Position in read from end'
    else:
        p = figure()
        p.grid.grid_line_alpha = 0.3
        p.line(
            x=range(x_length),
            y=np.array([np.mean(position)
                        for position in zip(*quallist)]),
            color='orange',
            legend="Quality")
        p.xaxis.axis_label = 'Position in read from start'
    p.yaxis.axis_label = 'Mean quality score of base calls'
    return p


if __name__ == "__main__":
    main()
