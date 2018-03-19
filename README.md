# nanoQC
Quality control tools for long read sequencing data aiming to replicate some of the plots made by fastQC.

[![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/wouter_decoster.svg?style=social&label=Follow%20%40wouter_decoster)](https://twitter.com/wouter_decoster)
[![install with conda](https://anaconda.org/bioconda/nanoqc/badges/installer/conda.svg)](https://anaconda.org/bioconda/nanoqc)
[![Code Health](https://landscape.io/github/wdecoster/nanoQC/master/landscape.svg?style=flat)](https://landscape.io/github/wdecoster/nanoQC/master)

Creates dynamic plots using [bokeh](https://bokeh.pydata.org/en/latest/).
For a static example see my blog [Gigabase or gigabyte](https://gigabaseorgigabyte.wordpress.com/2017/06/15/per-base-sequence-content-and-quality-end-of-reads/)


## INSTALLATION
```bash
pip install nanoQC
```
or  
[![install with conda](https://anaconda.org/bioconda/nanoqc/badges/installer/conda.svg)](https://anaconda.org/bioconda/nanoqc)
```
conda install -c bioconda nanoqc
```


## USAGE
```
nanoQC [-h] [-v] [--outdir OUTDIR] fastq

positional arguments:
  fastq            Reads data in fastq format.

optional arguments:
  -h, --help       show this help message and exit
  -v, --version    Print version and exit.
  -o, --outdir     Specify directory in which output has to be created.
```

## STATUS
[![Code Health](https://landscape.io/github/wdecoster/nanoQC/master/landscape.svg?style=flat)](https://landscape.io/github/wdecoster/nanoQC/master)

## CITATION
If you use this tool, please consider citing our [publication](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty149/4934939).
