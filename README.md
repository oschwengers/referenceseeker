[![DOI](https://joss.theoj.org/papers/10.21105/joss.01994/status.svg)](https://doi.org/10.21105/joss.01994)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/oschwengers/referenceseeker/blob/master/LICENSE)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/referenceseeker.svg)
![GitHub release](https://img.shields.io/github/release/oschwengers/referenceseeker.svg)
![PyPI](https://img.shields.io/pypi/v/referenceseeker.svg)
![PyPI - Status](https://img.shields.io/pypi/status/referenceseeker.svg)
[![Conda](https://img.shields.io/conda/v/bioconda/referenceseeker.svg)](http://bioconda.github.io/recipes/referenceseeker/README.html)
![Python package](https://github.com/oschwengers/referenceseeker/workflows/Python%20package/badge.svg?branch=master)

# ReferenceSeeker: rapid determination of appropriate reference genomes.

## Contents
- [Description](#description)
- [Input & Output](#input-output)
- [Installation](#installation)
  - [BioConda](#bioconda)
  - [GitHub](#github)
- [Usage](#usage)
- [Examples](#examples)
- [Databases](#databases)
  - [RefSeq](#refseq-based)
  - [Custom](#custom-database)
- [Dependencies](#dependencies)
- [Citation](#citation)

## Description
ReferenceSeeker determines closely related reference genomes following a scalable hierarchical
approach combining an fast kmer profile-based database lookup of candidate
reference genomes and subsequent computation of specific average nucleotide
identity (ANI) values for the rapid determination of suitable reference genomes.

ReferenceSeeker computes kmer-based genome distances between a query genome and
potential reference genome candidates via Mash (Ondov et al. 2016).
For resulting candidates ReferenceSeeker subsequently computes (bidirectional) ANI values picking
genomes meeting community standard thresholds by default (ANI >= 95 % & conserved DNA >= 69 %)
(Goris, Konstantinos et al. 2007) ranked by the product of ANI and conserved DNA values
to take into account both genome coverage and identity.

Custom databases can be built with local genomes. For further convenience, we provide pre-built
databases based on RefSeq's (<https://www.ncbi.nlm.nih.gov/refseq>) *complete*,
*reference* and *representative* genomes for the following microbial taxa:
- bacteria
- archaea
- fungi
- protozoa
- viruses

The reasoning for subsequent calculations of both ANI and conserved DNA values
is that Mash distance values correlate well with ANI values for closely
related genomes, however the same is not true for conserved DNA values. A kmer
fingerprint-based comparison alone cannot distinguish if a kmer is missing
due to a SNP, for instance or a lack of the kmer-comprising subsequence.
As DNA conservation (next to DNA identity) is very important for many kinds of
analyses, e.g. reference based SNP detections, ranking potential reference
genomes based on a mash distance alone is often not sufficient in order to select
the most appropriate reference genomes. If desired, ANI and conserved DNA values
can be computed bidirectionally.

![Mash D vs. ANI / conDNA](mash-ani-cdna.mini.png?raw=true)

## Input & Output
### Input:
Path to a taxon database and a draft or finished genome in fasta format:
```
$ referenceseeker ~/bacteria GCF_000013425.1.fna
```

### Output:
Tab separated lines to STDOUT comprising the following columns:

Unidirectionally (query -> references):
- RefSeq Assembly ID
- Mash Distance
- ANI
- Conserved DNA
- NCBI Taxonomy ID
- Assembly Status
- Organism

```
#ID	Mash Distance	ANI	Con. DNA	Taxonomy ID	Assembly Status	Organism
GCF_000013425.1	0.00000	100.00	100.00	93061	complete	Staphylococcus aureus subsp. aureus NCTC 8325
GCF_001900185.1	0.00002	100.00	99.89	46170	complete	Staphylococcus aureus subsp. aureus HG001
GCF_900475245.1	0.00004	100.00	99.57	93061	complete	Staphylococcus aureus subsp. aureus NCTC 8325 NCTC8325
GCF_001018725.2	0.00016	100.00	99.28	1280	complete	Staphylococcus aureus FDAARGOS_10
GCF_003595465.1	0.00185	99.86	96.81	1280	complete	Staphylococcus aureus USA300-SUR6
GCF_003595385.1	0.00180	99.87	96.80	1280	complete	Staphylococcus aureus USA300-SUR2
GCF_003595365.1	0.00180	99.87	96.80	1280	complete	Staphylococcus aureus USA300-SUR1
GCF_001956815.1	0.00180	99.87	96.80	46170	complete	Staphylococcus aureus subsp. aureus USA300_SUR1
...
```
Bidirectionally (query -> references [QR] & references -> query [RQ]):
- RefSeq Assembly ID
- Mash Distance
- QR ANI
- QR Conserved DNA
- RQ ANI
- RQ Conserved DNA
- NCBI Taxonomy ID
- Assembly Status
- Organism

```
#ID	Mash Distance	QR ANI	QR Con. DNA	RQ ANI	RQ Con. DNA	Taxonomy ID	Assembly Status	Organism
GCF_000013425.1	0.00000	100.00	100.00	100.00	100.00	93061	complete	Staphylococcus aureus subsp. aureus NCTC 8325
GCF_001900185.1	0.00002	100.00	99.89	100.00	99.89	46170	complete	Staphylococcus aureus subsp. aureus HG001
GCF_900475245.1	0.00004	100.00	99.57	99.99	99.67	93061	complete	Staphylococcus aureus subsp. aureus NCTC 8325 NCTC8325
GCF_001018725.2	0.00016	100.00	99.28	99.95	98.88	1280	complete	Staphylococcus aureus FDAARGOS_10
GCF_001018915.2	0.00056	99.99	96.35	99.98	99.55	1280	complete	Staphylococcus aureus NRS133
GCF_001019415.2	0.00081	99.99	94.47	99.98	99.36	1280	complete	Staphylococcus aureus NRS146
GCF_001018735.2	0.00096	100.00	94.76	99.98	98.58	1280	complete	Staphylococcus aureus NRS137
GCF_003354885.1	0.00103	99.93	96.63	99.93	96.66	1280	complete	Staphylococcus aureus 164
...
```

## Installation
ReferenceSeeker can be installed via Conda and Git(Hub).

In either case, a taxon database must be downloaded which we provide for download at Zenodo:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3562004.svg)](https://doi.org/10.5281/zenodo.3562004)
For more information scroll to [Databases](#databases).

### BioConda
The preferred way to install and run ReferenceSeeker is [Conda](https://conda.io/docs/install/quick.html) using the [Bioconda](https://bioconda.github.io/) channel:
```
$ conda install -c conda-forge -c bioconda -c defaults referenceseeker
$ referenceseeker --help
```

### GitHub
Alternatively, you can use this raw GitHub repository:
1. install necessary Python dependencies (if necessary)
2. clone the latest version of the repository
3. download and extract a databases

```
$ pip3 install --user biopython
$ git clone https://github.com/oschwengers/referenceseeker.git
$ ./referenceseeker/bin/referenceseeker --help
```

### Test
To test your installation we prepared a tiny mock database comprising 4 `Salmonella spp` genomes
and a query assembly (SRA: SRR498276) in the `tests` directory:
```
$ git clone https://github.com/oschwengers/referenceseeker.git
$
$ # BioConda installation
$ referenceseeker referenceseeker/tests/db referenceseeker/tests/Salmonella_enterica_CFSAN000189.fasta
$
$ # GitHub installation
$ ./referenceseeker/bin/referenceseeker referenceseeker/tests/db referenceseeker/tests/Salmonella_enterica_CFSAN000189.fasta
```

Expected output:
```
#ID	Mash Distance	ANI	Con. DNA	Taxonomy ID	Assembly Status	Organism
GCF_000439415.1	0.00003	100.00	99.55	1173427	complete	Salmonella enterica subsp. enterica serovar Bareilly str. CFSAN000189
GCF_002760915.1	0.01000	99.00	89.86	149539	complete	Salmonella enterica subsp. enterica serovar Enteritidis 56-3991
GCF_900205275.1	0.01522	98.61	83.13	90370	complete	Salmonella enterica subsp. enterica serovar Typhi
```

## Usage
Usage:
```
usage: referenceseeker [--crg CRG] [--ani ANI] [--conserved-dna CONSERVED_DNA]
                       [--unfiltered] [--bidirectional] [--help] [--version]
                       [--verbose] [--threads THREADS]
                       <database> <genome>

Rapid determination of appropriate reference genomes.

positional arguments:
  <database>            ReferenceSeeker database path
  <genome>              target draft genome in fasta format

Filter options / thresholds:
  These options control the filtering and alignment workflow.

  --crg CRG, -r CRG     Max number of candidate reference genomes to pass kmer
                        prefilter (default = 100)
  --ani ANI, -a ANI     ANI threshold (default = 0.95)
  --conserved-dna CONSERVED_DNA, -c CONSERVED_DNA
                        Conserved DNA threshold (default = 0.69)
  --unfiltered, -u      Set kmer prefilter to extremely conservative values
                        and skip species level ANI cutoffs (ANI >= 0.95 and
                        conserved DNA >= 0.69
  --bidirectional, -b   Compute bidirectional ANI/conserved DNA values
                        (default = False)

Runtime & auxiliary options:
  --help, -h            Show this help message and exit
  --version, -V         show program's version number and exit
  --verbose, -v         Print verbose information
  --threads THREADS, -t THREADS
                        Number of used threads (default = number of available
                        CPU cores)
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
ReferenceSeeker depends on custom databases comprising taxonomic genome informations as
well as kmer hash profiles for each entry.

### RefSeq based
We provide the following pre-built databases based on RefSeq release 99 (2020-03-06)
via [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3725706.svg)](https://doi.org/10.5281/zenodo.3725706) :

| Taxon | URL | # Genomes | Size Zipped | Size Unzipped |
| :---: | --- | ---: | :---: | :---: |
| bacteria | <https://zenodo.org/record/3725706/files/bacteria.tar.gz?download=1> | 24,879 | 32 Gb | 101 Gb |
| archaea | <https://zenodo.org/record/3725706/files/archaea.tar.gz?download=1> | 564 | 513 Mb | 1.6 Gb |
| fungi | <https://zenodo.org/record/3725706/files/fungi.tar.gz?download=1> | 299 | 2.7 Gb | 8.4 Gb |
| protozoa | <https://zenodo.org/record/3725706/files/protozoa.tar.gz?download=1> | 88 | 1 Gb | 3.4 Gb |
| viruses | <https://zenodo.org/record/3725706/files/viruses.tar.gz?download=1> | 8,999 | 602 Mb | 898 Mb |

Updated database versions reflecting the latest RefSeq versions can be built
with shell and nextflow scripts contained in this repository.

Download and install Nextflow:
```
$ curl -fsSL get.nextflow.io | bash
```

Clone this repo and build a database:
```
$ git clone git@github.com:oschwengers/referenceseeker.git
$ export REFERENCE_SEEKER_HOME=$(realpath referenceseeker)
$ sh $REFERENCE_SEEKER_HOME/db-scripts/build-db.sh <DB_TYPE_OPTION>
```

List of available database options:
```
$ sh build-db.sh
	-b (bacteria)
	-a (archaea)
	-v (viruses)
	-f (fungi)
	-p (protozoa)
```

### Custom database
If above mentiond RefSeq based databases do not contain sufficiently-close related
genomes or are just too large, ReferenceSeeker provides auxiliary commands
in order to either create databases from scratch or to expand existing ones.
Therefore, a second executable `referenceseeker_db` accepts `init` and `import` subcommands:

Usage:
```
usage: referenceseeker_db [--help] [--version] {init,import} ...

Rapid determination of appropriate reference genomes.

positional arguments:
  {init,import}  sub-command help
    init         Initialize a new database
    import       Add a new genome to database

Runtime & auxiliary options:
  --help, -h     Show this help message and exit
  --version, -V  show program's version number and exit
```

If a new database should be created, use `referenceseeker_db init`:
```
usage: referenceseeker_db init [-h] [--output OUTPUT] --db DB

optional arguments:
  -h, --help            show this help message and exit
  --output OUTPUT, -o OUTPUT
                        output directory (default = current working directory)
  --db DB, -d DB        Name of the new ReferenceSeeker database
```

This new database or an existing one can be used to import genomes in Fasta, GenBank or EMBL format:
```
usage: referenceseeker_db import [-h] --db DB --genome GENOME [--id ID]
                                 [--taxonomy TAXONOMY]
                                 [--status {complete,chromosome,scaffold,contig}]
                                 [--organism ORGANISM]

optional arguments:
  -h, --help            show this help message and exit
  --db DB, -d DB        ReferenceSeeker database path
  --genome GENOME, -g GENOME
                        Genome path [Fasta, GenBank, EMBL]
  --id ID, -i ID        Unique genome identifier (default sequence id of first
                        record)
  --taxonomy TAXONOMY, -t TAXONOMY
                        Taxonomy ID (default = 12908 [unclassified sequences])
  --status {complete,chromosome,scaffold,contig}, -s {complete,chromosome,scaffold,contig}
                        Assembly level (default = contig)
  --organism ORGANISM, -o ORGANISM
                        Organism name (default = "")
```

## Dependencies
ReferenceSeeker needs the following dependencies:
- Python (3.5.2), Biopython (1.71)
- Mash (2.2) <https://github.com/marbl/Mash>
- MUMmer (4.0.0-beta2) <https://github.com/gmarcais/mummer>

ReferenceSeeker has been tested against aforementioned versions.

## Citation

> Schwengers et al., (2020). ReferenceSeeker: rapid determination of appropriate reference genomes. Journal of Open Source Software, 5(46), 1994, https://doi.org/10.21105/joss.01994
