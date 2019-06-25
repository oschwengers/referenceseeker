# ReferenceSeeker: rapid determination of appropriate reference genomes.
Author: Oliver Schwengers (oliver.schwengers@computational.bio.uni-giessen.de)


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
genomes meeting community standard thresholds (ANI >= 95 % & conserved DNA >= 69 %)
(Goris, Konstantinos et al. 2007) ranked by ANI and conserved DNA.
Additionally, ReferenceSeeker can use MeDuSa (Bosi, Donati et al. 2015)
to scaffold contigs based on the 20 most suitable reference genomes.

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
$ referenceseeker.py --db ~/bacteria GCF_000013425.1.fna
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
To setup ReferenceSeeker just do the following:
1. install necessary Python / Java dependencies (if necessary)
2. clone the latest version of the repository
3. set REFERENCE_SEEKER_HOME environment variable pointing to the repository directory
4. download and extract a databases or create one yourself

Example:
```
$ cd ~
$ pip3 install biopython numpy networkx
$ sudo apt-get install openjdk-8-jdk
$ git clone https://github.com/oschwengers/referenceseeker.git
$ export REFERENCE_SEEKER_HOME=~/referenceseeker
$ wget https://s3.computational.bio.uni-giessen.de/swift/v1/referenceseeker/bacteria.tar.gz
$ tar -xzf bacteria.tar.gz
$ rm bacteria.tar.gz
```

Alternatively, use the aforementioned Docker image (oschwengers/referenceseeker)
in order to ease the setup process.


## Usage
Usage:
```
ReferenceSeeker [-h] --db DB [--threads THREADS] [--unfiltered] [--verbose] [--scaffolds] [--output OUTPUT] [--version] <genome>

positional arguments:
  <genome>              Target draft genome in fasta format

optional arguments:
  -h, --help            show this help message and exit
  --db DB, -d DB        ReferenceSeeker database path
  --crg CRG, -c CRG     Max number of candidate reference genomes to assess
                        (default = 100)
  --unfiltered, -u      Set kmer prefilter to extremely conservative values
                        and skip species level ANI cutoffs (ANI >= 0.95 and
                        conserved DNA >= 0.69
  --scaffolds, -s       Build scaffolds via MeDuSa (Bosi, Donati et al. 2015)
                        based on detected references
  --output OUTPUT, -o OUTPUT
                        Output fasta file for built scaffolds
  --threads THREADS, -t THREADS
                        Number of threads to use (default = number of
                        available CPUs)
  --verbose, -v         Print verbose information
  --version             show program's version number and exit
```

## Examples
Simple:
```
$ $REFERENCE_SEEKER_HOME/referenceseeker.py --db <REFERENCE_SEEKER_DB> <GENOME>
```

Expert: creating scaffolds with verbose output using a defined number of threads:
```
$ $REFERENCE_SEEKER_HOME/referenceseeker.py --db <REFERENCE_SEEKER_DB> --crg 500 --scaffolds --output scaffolds.fasta --verbose --threads 8 <GENOME>
```

With Docker:
```
$ sudo docker pull oschwengers/referenceseeker:latest
$ sudo docker run --rm -v <REFERENCE_SEEKER_DB>:/db -v <DATA_DIR>:/data oschwengers/referenceseeker:latest <GENOME>
```

With Docker shell script:
```
$ sudo docker pull oschwengers/referenceseeker:latest
$ referenceseeker.sh <REFERENCE_SEEKER_DB> <GENOME>
```


## Databases
ReferenceSeeker depends on custom databases based on reference, representative as
well as complete NCBI RefSeq genomes comprising kmer hash profiles taxonomic information.
We provide the following pre-built databases based on RefSeq release 90:

| Taxon | URL | # Genomes | Size Zipped | Size Unzipped |
| :---: | --- | ---: | :---: | :---: |
| bacteria | <https://s3.computational.bio.uni-giessen.de/swift/v1/referenceseeker/bacteria.tar.gz> | 14,983 | 18 Gb | 58 Gb |
| archaea | <https://s3.computational.bio.uni-giessen.de/swift/v1/referenceseeker/archaea.tar.gz> | 386 | 335 Mb | 1.1 Gb |
| fungi | <https://s3.computational.bio.uni-giessen.de/swift/v1/referenceseeker/fungi.tar.gz> | 277 | 2.5 Gb | 7.7 Gb |
| protozoa | <https://s3.computational.bio.uni-giessen.de/swift/v1/referenceseeker/protozoa.tar.gz> | 77 | 953 Mb | 3.2 Gb |
| viral | <https://s3.computational.bio.uni-giessen.de/swift/v1/referenceseeker/viral.tar.gz> | 7,855 | 525 Mb | 719 Mb |

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
- Python (3.5.2), Biopython (1.71), NumPy (1.14.5), NetworkX (1.11)
- JAVA (8)
- Mash (2.0) <https://github.com/marbl/Mash>
- MUMmer (4.0.0-beta2) <https://github.com/gmarcais/mummer>
- MeDuSa (1.6) <https://github.com/combogenomics/medusa>

ReferenceSeeker has been tested against aforementioned versions.


## Citation
To temporarily cite our work, please transitionally refer to:
> Schwengers O., Hain T., Chakraborty T., Goesmann A. (2019) ReferenceSeeker: rapid determination of appropriate reference genomes. GitHub https://github.com/oschwengers/referenceseeker
