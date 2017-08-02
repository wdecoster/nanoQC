# nanoQC
Quality control tools for Oxford Nanopore sequencing data aiming to replicate some of the plots made by fastQC.


For an example see my blog [Gigabase or gigabyte](https://gigabaseorgigabyte.wordpress.com/2017/06/15/per-base-sequence-content-and-quality-end-of-reads/)

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
  --outdir OUTDIR  Specify directory in which output has to be created.
```
