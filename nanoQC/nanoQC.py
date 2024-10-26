# wdecoster
# Edited by j.ouwerkerk

import os
import sys
from argparse import ArgumentParser
import logging
from Bio import SeqIO
import numpy as np
from .version import __version__
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def main():
    args = get_args()
    make_output_dir(args.outdir)
    size_range = int(args.minlen / 2)
    logging.basicConfig(
        format="%(asctime)s %(message)s",
        filename=os.path.join(args.outdir, "NanoQC.log"),
        level=logging.INFO,
    )
    logging.info("NanoQC started.")
    hist = length_histogram(fqin=compressed_input(args.fastq))
    fq = get_bin(fq=compressed_input(args.fastq), size_range=size_range)
    if len(fq) == 0:
        logging.critical(f"No reads with a higher length of {size_range * 2}.")
        logging.info("Exiting...")
    else:
        logging.info(f"Using {len(fq)} reads with a minimum length of {size_range * 2}bp for plotting")
        logging.info("Creating plots...")
        seq_plots, qual_plots = per_base_sequence_content_and_quality(
            head_seq=[dat[0] for dat in fq],
            head_qual=[dat[1] for dat in fq],
            tail_seq=[dat[2] for dat in fq],
            tail_qual=[dat[3] for dat in fq],
            rna=args.rna,
        )
        # Create a subplot figure
        fig = make_subplots(
            rows=3,
            cols=2,
            specs=[[{"colspan": 2}, None], [{}, {}], [{}, {}]],
            subplot_titles=(
                "Read length distribution",
                "Nucleotide diversity from start",
                "Nucleotide diversity from end",
                "Quality score from start",
                "Quality score from end",
            ),
            vertical_spacing=0.1,
            shared_yaxes=True,
        )

        # Add the histogram plot
        fig.add_trace(hist.data[0], row=1, col=1)

        # Add the sequence content plots
        for c, seq_plot in enumerate(seq_plots, start=1):
            for trace in seq_plot.data:
                fig.add_trace(trace, row=2, col=c)

        # Add the quality plots
        for c, qual_plot in enumerate(qual_plots, start=1):
            for trace in qual_plot.data:
                fig.add_trace(trace, row=3, col=c)
        # Update layout
        fig.update_layout(height=1600, width=1000, title_text="nanoQC Report")
        # put the legend at the bottom
        fig.update_layout(legend=dict(orientation="h", y=-0.05, x=0.5, xanchor="center"))
        fig.write_html(os.path.join(args.outdir, "nanoQC.html"), include_plotlyjs="cdn")
        logging.info("Finished!")


def get_args():
    parser = ArgumentParser(
        description="Investigate nucleotide composition and base quality."
    )
    parser.add_argument(
        "-v",
        "--version",
        help="Print version and exit.",
        action="version",
        version=f"NanoQC {__version__}",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        help="Specify directory in which output has to be created.",
        default=".",
    )
    parser.add_argument(
        "--rna",
        help="Fastq is from direct RNA-seq and contains U nucleotides.",
        action="store_true",
    )
    parser.add_argument(
        "-l",
        "--minlen",
        help=(
            "Filters the reads on a minimal length of the given range.\n"
            "Also plots the given length/2 of the begin and end of the reads."
        ),
        default=200,
        type=int,
    )
    parser.add_argument("fastq", help="Reads data in fastq.gz format.")
    return parser.parse_args()


def make_output_dir(path):
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    except IOError:
        sys.exit("ERROR: No writing permission to the output directory.")


def compressed_input(inputfq):
    if inputfq.endswith(".gz"):
        import gzip

        return gzip.open(inputfq, "rt")
    elif inputfq.endswith(".bz2"):
        import bz2

        return bz2.open(inputfq, "rt")
    elif inputfq.endswith(
        (
            ".fastq",
            ".fq",
        )
    ):
        return open(inputfq, "r")
    else:
        sys.exit(
            f"INPUT ERROR:\nUnrecognized file extension for {inputfq}\n"
            "Supported are .gz, .bz2, .fastq and .fq"
        )


def length_histogram(fqin):
    """
    Create a histogram
    """
    lengths = get_lengths(fqin)
    hist, edges = np.histogram(lengths, bins="auto")

    hist_norm = go.Figure()
    hist_norm.add_trace(
        go.Bar(
            x=edges[:-1],
            y=hist,
            marker=dict(color="#036564"),
            name="Read Length",
            showlegend=False,
        )
    )
    hist_norm.update_layout(
        xaxis=dict(tickformat=".0f"),
        yaxis_title="Count",
        xaxis_title="Read Length",
        title="Read Length Distribution",
    )
    return hist_norm


def get_lengths(fastq):
    """
    Loop over the fastq file, extract length of sequences
    """
    return np.array([len(record) for record in SeqIO.parse(fastq, "fastq")])


def per_base_sequence_content_and_quality(
    head_seq, head_qual, tail_seq, tail_qual, rna=False
):
    seq_plot_left = plot_nucleotide_diversity(head_seq, rna=rna)
    seq_plot_right = plot_nucleotide_diversity(tail_seq, invert=True, rna=rna)
    qual_plot_left = plot_qual(head_qual)
    qual_plot_right = plot_qual(tail_qual, invert=True)
    logging.info("Per base sequence content and quality completed.")
    return [seq_plot_left, seq_plot_right], [qual_plot_left, qual_plot_right]


def get_bin(fq, size_range):
    """
    Loop over the fastq file
    Extract list of nucleotides and list of quality scores in tuples in list
    Only select those reads of which the length is within the size range
    """
    logging.info("Extracting nucleotides and quality scores.")
    return [
        (
            list(rec.seq)[:size_range],
            list(rec.letter_annotations["phred_quality"])[:size_range],
            list(rec.seq[-1 * size_range :]),
            list(rec.letter_annotations["phred_quality"])[-1 * size_range :],
        )
        for rec in SeqIO.parse(fq, "fastq")
        if len(rec) >= size_range * 2
    ]


def plot_nucleotide_diversity(seqs, invert=False, rna=False):
    x_length = len(seqs[0])
    numreads = len(seqs)
    x = list(range(x_length))
    if invert:
        x = [-1*i for i in list(reversed(x))]

    fig = go.Figure()

    for nucl, color in zip(
        ["A", "U" if rna else "T", "G", "C"],
        ["green", "pink" if rna else "red", "black", "blue"],
    ):
        y = [pos.count(nucl) / numreads for pos in zip(*seqs)]
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                mode="lines",
                name=nucl,
                line=dict(color=color),
                showlegend=True if invert else False,
            )
        )

    fig.update_layout(
        xaxis_title=(
            "Position in read from end" if invert else "Position in read from start"
        ),
        yaxis_title="Frequency of nucleotide in read",
        title="Nucleotide Diversity",
    )

    return fig


def get_qual_per_pos(quallist):
    """
    Plot average quality per position
    params:
    -quallist: List of the qualities per read.
    -sizeRange: Range of the first or last bases in the reads.
    """
    mean_qual_per_pos = []
    # Extract the first positions to the last positions.
    for position in range(len(quallist[0])):
        # Make a list with quality per pos in all the reads in quallist.
        pos_x_list = [quality[position] for quality in quallist]
        # Calculate the mean quality of the position and append it to a list.
        mean_qual_per_pos.append(np.mean(pos_x_list))
    return mean_qual_per_pos


def plot_qual(quallist, invert=False):
    """
    Create a FastQC-like "Per base sequence quality" plot using Plotly
    Plot average quality per position
    zip will stop when shortest read is exhausted
    """
    mean_quallist = get_qual_per_pos(quallist)
    x_length = len(mean_quallist)
    x = list(range(x_length))
    if invert:
        x = [-1*i for i in list(reversed(x))]


    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x,
            y=mean_quallist,
            mode="lines",
            line=dict(color="orange"),
            showlegend=False,
        )
    )

    fig.update_layout(
        xaxis_title=(
            "Position in read from end" if invert else "Position in read from start"
        ),
        yaxis_title="Mean quality score of base calls",
        title="Per Base Sequence Quality",
    )

    return fig


if __name__ == "__main__":
    main()
