# nanoQC
Quality control tools for long read sequencing data aiming to replicate some of the plots made by fastQC.

[![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/wouter_decoster.svg?style=social&label=Follow%20%40wouter_decoster)](https://twitter.com/wouter_decoster)
[![install with conda](https://anaconda.org/bioconda/nanoqc/badges/installer/conda.svg)](https://anaconda.org/bioconda/nanoqc)


Creates dynamic plots using [bokeh](https://bokeh.pydata.org/en/latest/).
For an example see [here](http://decoster.xyz/wouter/)


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
