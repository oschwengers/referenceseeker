[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/oschwengers/referenceseeker/blob/master/LICENSE)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/referenceseeker.svg)
![GitHub release](https://img.shields.io/github/release/oschwengers/referenceseeker.svg)
![PyPI](https://img.shields.io/pypi/v/referenceseeker.svg)
![PyPI - Status](https://img.shields.io/pypi/status/referenceseeker.svg)
![Conda](https://img.shields.io/conda/v/bioconda/referenceseeker.svg)
![Conda](https://img.shields.io/conda/pn/bioconda/referenceseeker.svg)

# ReferenceSeeker: rapid determination of appropriate reference genomes.

## Contents
- [Description](#description)
- [Input & Output](#input-output)
- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [Databases](#databases)
- [Dependencies](#dependencies)
- [Citation](#citation)

## Description
ReferenceSeeker determines closely related reference genomes from
RefSeq (<https://www.ncbi.nlm.nih.gov/refseq>) following a scalable hierarchical
approach combining an ultra-fast kmer profile-based database lookup of candidate
reference genomes and subsequent computation of specific average nucleotide
identity (ANI) values for the rapid determination of suitable reference genomes.

ReferenceSeeker computes kmer-based genome distances between a query genome and
and a database built on RefSeq genomes via Mash (Ondov et al. 2016). Therefore,
only complete genomes or those stated as 'representative' or 'reference' genome
are included. ReferenceSeeker offers pre-built databases for a broad spectrum of
microbial taxonomic groups, i.e. bacteria, archaea, fungi, protozoa and viruses.
For resulting candidates ReferenceSeeker subsequently computes ANI values picking
genomes meeting community standard thresholds by default (ANI >= 95 % & conserved DNA >= 69 %)
(Goris, Konstantinos et al. 2007) ranked by ANI and conserved DNA.

The reasoning for subsequent calculations of both ANI and conserved DNA values
is that Mash distance values correlate well with ANI values for closely
related genomes but the same is not true for conserved DNA values. A kmer
fingerprint-based comparison alone cannot distinguish if a kmer is missing
due to a SNP, for instance or a lack of the kmer-comprising subsequence.
As DNA conservancy (next to DNA identity) is very important for many kinds of
analyses, e.g. reference based SNP detections, ranking potential reference
genomes based on a mash distance alone is often not sufficient in order to select
the most appropriate reference genomes.

![Mash D vs. ANI / conDNA](mash-ani-cdna.mini.png?raw=true)

## Input & Output
### Input:
Path to a taxon database and a draft or finished genome in fasta format:
```
$ referenceseeker ~/bacteria GCF_000013425.1.fna
```

### Output:
Tab separated lines to STDOUT comprising the following columns:
- RefSeq Assembly ID
- ANI
- Conserved DNA
- Mash Distance
- NCBI Taxonomy ID
- Assembly Status
- Organism

```
#ID	ANI	Con. DNA	Mash Distance	Taxonomy ID	Assembly Status	Organism
GCF_000013425.1	100.00	100.00	0.00000	93061	complete	Staphylococcus aureus subsp. aureus NCTC 8325
GCF_001900185.1	100.00	99.89	0.00002	46170	complete	Staphylococcus aureus subsp. aureus HG001
GCF_900475245.1	100.00	99.57	0.00004	93061	complete	Staphylococcus aureus subsp. aureus NCTC 8325 NCTC8325
GCF_001018725.2	100.00	99.28	0.00016	1280	complete	Staphylococcus aureus FDAARGOS_10
GCF_001018915.2	99.99	96.35	0.00056	1280	complete	Staphylococcus aureus NRS133
...
```

## Installation
Platon can be installed and used in 2 different ways.

In either case, a taxon database must be downloaded which we provide for download at Zenodo:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3562005.svg)](https://doi.org/10.5281/zenodo.3562005)
For more information scroll to [Databases](#databases).

### Conda:
The preferred way to install and run ReferenceSeeker is BioConda:
```
$ conda install -c conda-forge -c bioconda -c defaults referenceseeker
$ referenceseeker --help
```

### GitHub
Alternatively, you can use this raw GitHub repository:
1. install necessary Python dependencies (if necessary)
2. clone the latest version of the repository
3. download and extract a databases

Example:
```
$ pip3 install --user biopython
$ git clone https://github.com/oschwengers/referenceseeker.git
$ ./referenceseeker/bin/referenceseeker --help
```

## Usage
Usage:
```
usage: referenceseeker [-h] [--crg CRG] [--ani ANI]
                       [--conserved-dna CONSERVED_DNA] [--unfiltered]
                       [--verbose] [--threads THREADS] [--version]
                       <database> <genome>

Rapid determination of appropriate reference genomes.

positional arguments:
  <database>            ReferenceSeeker database path
  <genome>              target draft genome in fasta format

optional arguments:
  -h, --help            show this help message and exit
  --crg CRG, -r CRG     max number of candidate reference genomes to pass kmer
                        prefilter (default = 100)
  --ani ANI, -a ANI     ANI threshold value (default = 0.95)
  --conserved-dna CONSERVED_DNA, -c CONSERVED_DNA
                        Conserved DNA threshold value (default = 0.69)
  --unfiltered, -u      set kmer prefilter to extremely conservative values
                        and skip species level ANI cutoffs (ANI >= 0.95 and
                        conserved DNA >= 0.69
  --verbose, -v         print verbose information
  --threads THREADS, -t THREADS
                        number of threads to use (default = number of
                        available CPUs)
  --version, -V         show program's version number and exit
```

## Examples
Simple:
```
$ # referenceseeker <REFERENCE_SEEKER_DB> <GENOME>
$ referenceseeker bacteria/ genome.fasta
```

Expert: verbose output and increased output of candidate reference genomes using a defined number of threads:
```
$ # referenceseeker --crg 500 --verbose --threads 8 <REFERENCE_SEEKER_DB> <GENOME>
$ referenceseeker --crg 500 --verbose --threads 8 bacteria/ genome.fasta
```

## Databases
ReferenceSeeker depends on custom databases based on reference, representative as
well as complete NCBI RefSeq genomes comprising kmer hash profiles and taxonomic information.
We provide the following pre-built databases based on RefSeq 2019-07-02 via [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3562005.svg)](https://doi.org/10.5281/zenodo.3562005) :

| Taxon | URL | # Genomes | Size Zipped | Size Unzipped |
| :---: | --- | ---: | :---: | :---: |
| bacteria | <https://zenodo.org/record/3562005/files/bacteria.tar.gz?download=1> | 18,229 | 22 Gb | 71 Gb |
| archaea | <https://zenodo.org/record/3562005/files/archaea.tar.gz?download=1> | 417 | 364 Mb | 1.2 Gb |
| fungi | <https://zenodo.org/record/3562005/files/fungi.tar.gz?download=1> | 288 | 2.6 Gb | 8 Gb |
| protozoa | <https://zenodo.org/record/3562005/files/protozoa.tar.gz?download=1> | 88 | 1 Gb | 3.4 Gb |
| viral | <https://zenodo.org/record/3562005/files/viral.tar.gz?download=1> | 9,264 | 608 Mb | 835 Mb |

Updated database versions reflecting the latest RefSeq versions can be built
with a shell script and nextflow pipeline.

Download and install Nextflow:
```
$ curl -fsSL get.nextflow.io | bash
```

Build database:
```
$ export REFERENCE_SEEKER_HOME=<REFERENCE_SEEKER_DIR>
$ sh $REFERENCE_SEEKER_HOME/build-db.sh <DB_TYPE_OPTION>
```

List of available database options:
```
$ sh build-db.sh
	-b (bacteria)
	-a (archaea)
	-v (viral)
	-f (fungi)
	-p (protozoa)
```

## Dependencies
ReferenceSeeker needs the following dependencies:
- Python (3.5.2), Biopython (1.71)
- Mash (2.2) <https://github.com/marbl/Mash>
- MUMmer (4.0.0-beta2) <https://github.com/gmarcais/mummer>

ReferenceSeeker has been tested against aforementioned versions.

## Citation

> ReferenceSeeker: rapid determination of appropriate reference genomes. Oliver Schwengers, Torsten Hain, Trinad Chakraborty, Alexander Goesmann. bioRxiv 863621; doi: https://doi.org/10.1101/863621
