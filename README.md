PHASM: Haplotype-aware *de novo* genome assembly for polyploid organisms
========================================================================

[![Build 
Status](https://travis-ci.org/lrvdijk/phasm.svg?branch=master)](https://travis-ci.org/lrvdijk/phasm)

PHASM is a long read *de novo* genome assembler that phases variants among 
chromosome homologues during the assembly process, and aims to output separate 
contigs for each haplotype. The main idea in PHASM is to build bubble chains: 
consecutive "superbubbles" chained together. While most traditional genome 
assemblers pop these superbubbles by only keeping the best supported path, 
PHASM finds *k* paths through this chain of superbubbles that best represent 
each haplotype.

This program has been created as part of my master thesis project. For now, it 
has only been tested with error free data.

Requirements
------------

* Python >= 3.5
* NumPy >= 1.11
* SciPy >= 0.16
* NetworkX >= 1.9
* (tests) pytest

Installation
------------

PHASM is not on PyPI yet, so for now you'll have to clone this repository and 
run:

    pip install -r requirements.txt
    python setup.py install

Related Repositories
--------------------

* [aneusim][aneusim]: a tool to generate synthetic aneuploid/polyploid genomes
* [phasm-benchmarks][phasm-benchmarks]: a complete [snakemake][snakemake] 
  pipeline that starts with finding pairwise local alignments between reads 
  using DALIGNER, and ends with an assembled polyploid genome.

[aneusim]: https://github.com/lrvdijk/aneusim
[phasm-benchmarks]: https://github.com/lrvdijk/phasm-benchmarks
[snakemake]: https://bitbucket.org/snakemake/snakemake

Documentation
=============

TODO

