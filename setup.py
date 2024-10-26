# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
exec(open("nanoQC/version.py").read())

setup(
    name="nanoQC",
    version=__version__,
    description="Create fastQC-like plots for Oxford Nanopore sequencing data",
    long_description=open(path.join(here, "README.md")).read(),
    long_description_content_type="text/markdown",
    url="https://github.com/wdecoster/nanoQC",
    author="Wouter De Coster",
    author_email="decosterwouter@gmail.com",
    license="GPLv3",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
    ],
    keywords="nanopore sequencing statistics",
    packages=find_packages() + ["scripts"],
    python_requires=">=3",
    install_requires=["plotly>=5.4.0", "numpy", "biopython"],
    package_data={"nanoQC`": []},
    package_dir={"nanoQC": "nanoQC"},
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "nanoQC=nanoQC.nanoQC:main",
        ],
    },
)
