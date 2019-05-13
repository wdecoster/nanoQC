# wdecoster
# Edited by j.ouwerkerk

import os
import sys
from argparse import ArgumentParser
import logging
from Bio import SeqIO
import numpy as np
from bokeh.plotting import figure, save, output_file
from bokeh.layouts import gridplot
from bokeh.models import Range1d
from .version import __version__


def main():
    args = get_args()
    make_output_dir(args.outdir)
    size_range = int(args.minlen / 2)
    logging.basicConfig(
        format='%(asctime)s %(message)s',
        filename=os.path.join(args.outdir, "NanoQC.log"),
        level=logging.INFO)
    logging.info("NanoQC started.")
    hist = length_histogram(fqin=compressed_input(args.fastq))
    fq = get_bin(fq=compressed_input(args.fastq),
                 size_range=size_range)
    if len(fq) == 0:
        logging.critical(
            "No reads with a higher length of {}.".format(size_range * 2))
        logging.info("Exiting...")
    else:
        logging.info(("Using {} reads with a minimum length of {}bp for plotting")
                     .format(len(fq), size_range * 2))
        logging.info("Creating plots...")
        seq_plots, qual_plots = per_base_sequence_content_and_quality(
            head_seq=[dat[0] for dat in fq],
            head_qual=[dat[1] for dat in fq],
            tail_seq=[dat[2] for dat in fq],
            tail_qual=[dat[3] for dat in fq],
            rna=args.rna)
        output_file(os.path.join(args.outdir, "nanoQC.html"), title="nanoQC_report")
        save(gridplot(children=[[hist], seq_plots, qual_plots],
                      plot_width=400,
                      plot_height=400))
        logging.info("Finished!")


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
    parser.add_argument("--rna",
                        help="Fastq is from direct RNA-seq and contains U nucleotides.",
                        action="store_true")
    parser.add_argument("-l", "--minlen",
                        help=("Filters the reads on a minimal length of the given range.\n"
                              "Also plots the given length/2 of the begin and end of the reads."),
                        default=200)
    parser.add_argument("fastq",
                        help="Reads data in fastq.gz format.")
    return parser.parse_args()


def make_output_dir(path):
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    except IOError:
        sys.exit("ERROR: No writing permission to the output directory.")


def compressed_input(inputfq):
    if inputfq.endswith('.gz'):
        import gzip
        return gzip.open(inputfq, 'rt')
    elif inputfq.endswith('.bz2'):
        import bz2
        return bz2.open(inputfq, 'rt')
    elif inputfq.endswith(('.fastq', '.fq',)):
        return open(inputfq, 'r')
    else:
        sys.exit('INPUT ERROR:\nUnrecognized file extension in {}\n'
                 'Supported are .gz, .bz2, .fastq and .fq'.format(inputfq))


def length_histogram(fqin):
    '''
    Create a histogram, and return the bin edges of the bin containing the most reads
    '''
    lengths = get_lengths(fqin)
    hist, edges = np.histogram(lengths, bins='auto')

    hist_norm = figure()
    hist_norm.quad(
        top=hist,
        bottom=0,
        left=edges[:-1],
        right=edges[1:],
        fill_color="#036564",
        line_color="#033649")
    hist_norm.xaxis[0].formatter.use_scientific = False
    return hist_norm


def get_lengths(fastq):
    '''
    Loop over the fastq file, extract length of sequences
    '''
    return np.array([len(record) for record in SeqIO.parse(fastq, "fastq")])


def per_base_sequence_content_and_quality(head_seq, head_qual, tail_seq, tail_qual, rna=False):
    seq_plot_left = plot_nucleotide_diversity(head_seq, rna=rna)
    seq_plot_right = plot_nucleotide_diversity(tail_seq, invert=True, rna=rna)
    qual_plot_left = plot_qual(head_qual)
    qual_plot_right = plot_qual(tail_qual, invert=True)
    logging.info("Per base sequence content and quality completed.")
    return [seq_plot_left, seq_plot_right], [qual_plot_left, qual_plot_right]


def get_bin(fq, size_range):
    '''
    Loop over the fastq file
    Extract list of nucleotides and list of quality scores in tuples in list
    Only select those reads of which the length is within the size range
    '''
    logging.info("Extracting nucleotides and quality scores.")
    return [(list(rec.seq)[:size_range],
             list(rec.letter_annotations["phred_quality"])[:size_range],
             list(rec.seq[-1 * size_range:]),
             list(rec.letter_annotations["phred_quality"])[-1 * size_range:])
            for rec in SeqIO.parse(fq, "fastq") if len(rec) >= size_range * 2]


def plot_nucleotide_diversity(seqs, invert=False, rna=False):
    x_length = len(seqs[0])
    if invert:
        p = figure(x_range=Range1d(start=x_length, end=0))
    else:
        p = figure()
    p.grid.grid_line_alpha = 0.3
    numreads = len(seqs)
    for nucl, color in zip(['A', 'U' if rna else 'T', 'G', 'C'],
                           ['green', 'pink' if rna else 'red', 'black', 'blue']):
        if invert:
            p.xaxis.axis_label = 'Position in read from end'
            p.line(
                x=range(x_length),
                y=list(reversed(np.array([pos.count(nucl) / numreads for pos in zip(*seqs)]))),
                color=color,
                legend=nucl)
        else:
            p.xaxis.axis_label = 'Position in read from start'
            p.line(
                x=range(x_length),
                y=np.array([pos.count(nucl) / numreads for pos in zip(*seqs)]),
                color=color,
                legend=nucl)
    p.yaxis.axis_label = 'Frequency of nucleotide in read'
    return p


def get_qual_per_pos(quallist):
    '''
    Plot average quality per position
    params:
    -quallist: List of the qualities per read.
    -sizeRange: Range of the first or last bases in the reads.
    '''
    mean_qual_per_pos = []
    # Extract the first positions to the last positions.
    for position in range(len(quallist[0])):
        # Make a list with quality per pos in all the reads in quallist.
        pos_x_list = [quality[position] for quality in quallist]
        # Calculate the mean quality of the position and append it to a list.
        mean_qual_per_pos.append(np.mean(pos_x_list))
    return mean_qual_per_pos


def plot_qual(quallist, invert=False):
    '''
    Create a FastQC-like "￼Per base sequence quality￼" plot
    Plot average quality per position
    zip will stop when shortest read is exhausted
    '''
    mean_quallist = get_qual_per_pos(quallist)
    x_length = len(mean_quallist)
    if invert:
        p = figure(x_range=Range1d(start=x_length, end=0))
        p.grid.grid_line_alpha = 0.3
        p.line(x=range(x_length),
               y=list(reversed(mean_quallist)),
               color='orange')
        p.xaxis.axis_label = 'Position in read from end'
    else:
        p = figure()
        p.grid.grid_line_alpha = 0.3
        p.line(x=range(x_length),
               y=mean_quallist,
               color='orange')
        p.xaxis.axis_label = 'Position in read from start'
    p.yaxis.axis_label = 'Mean quality score of base calls'
    return p


if __name__ == "__main__":
    main()
