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
to scaffold contigs based on the 20 closest reference genomes.


## Input & Output
### Input:
Path to a taxon database and a draft or finished genome in fasta format:
```
referenceseeker.py --db ~/bacteria GCF_000013425.1.fna
```

### Output:
Tab separated lines to STDOUT comprising the following columns:
- RefSeq Assembly ID
- ANI
- Conserved DNA
- NCBI Taxonomy ID
- Assembly Status
- Organism

```
#ID	ANI Con. DNA Taxonomy ID Assembly Status Organism
GCF_000013425.1	 100.00	 100.00	93061	complete	Staphylococcus aureus subsp. aureus NCTC 8325
GCF_001900185.1	 100.00	 99.89	46170	complete	Staphylococcus aureus subsp. aureus HG001
GCF_900475245.1	 100.00	 99.57	93061	complete	Staphylococcus aureus subsp. aureus NCTC 8325 NCTC8325
GCF_001018725.2	 100.00	 99.28	1280	complete	Staphylococcus aureus FDAARGOS_10
GCF_003595385.1	 99.87	 96.80	1280	complete	Staphylococcus aureus USA300-SUR2
...
```

## Installation
To setup ReferenceSeeker just do the following:
1. install necessary Python dependencies via pip
2. clone the latest version of the repository
3. set REFERENCE_SEEKER_HOME environment variable pointing to the repository directory
4. download and extract a databases or create one yourself

Example:
```
cd ~
pip3 install biopython numpy networkx
git clone https://github.com/oschwengers/referenceseeker.git
export REFERENCE_SEEKER_HOME=~/referenceseeker
wget https://s3.computational.bio.uni-giessen.de/swift/v1/referenceseeker/bacteria.tar.gz
tar -xzf bacteria.tar.gz
rm bacteria.tar.gz
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
  --threads THREADS, -t THREADS
                        Number of threads to use (default = number of
                        available CPUs)
  --unfiltered, -u      Set kmer prefilter to extremely conservative values
                        and skip species level ANI cutoffs (ANI >= 0.95 and
                        conserved DNA >= 0.69. Use this option for very rare species.
  --verbose, -v         Print verbose information
  --scaffolds, -s       Build scaffolds via MeDuSa (Bosi, Donati et al. 2015)
                        based on detected references
  --output OUTPUT, -o OUTPUT
                        Output fasta file for built scaffolds
  --version             show program's version number and exit
```

## Examples
Simple:
```
$REFERENCE_SEEKER_HOME/referenceseeker.py --db <REFERENCE_SEEKER_DB> <GENOME>
```

Expert: creating scaffolds with verbose output using a defined number of threads:
```
$REFERENCE_SEEKER_HOME/referenceseeker.py --db <REFERENCE_SEEKER_DB> --scaffolds --output scaffolds.fasta --verbose --threads 8 <GENOME>
```

With Docker:
```
sudo docker pull oschwengers/referenceseeker:latest
sudo docker run --rm -v <REFERENCE_SEEKER_DB>:/db -v <DATA_DIR>:/data oschwengers/referenceseeker:latest <GENOME>
```

With Docker shell script:
```
sudo docker pull oschwengers/referenceseeker:latest
referenceseeker.sh <REFERENCE_SEEKER_DB> <GENOME>
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
curl -fsSL get.nextflow.io | bash
```

Build database:
```
export REFERENCE_SEEKER_HOME=<REFERENCE_SEEKER_DIR>
sh $REFERENCE_SEEKER_HOME/build-db.sh <DB_TYPE_OPTION>
```

`build-db.sh -h` prints a list of available database options:
- `-b`: bacteria
- `-a`: archaea
- `-f`: fungi
- `-p`: protozoa
- `-v`: viruses

## Dependencies
ReferenceSeeker depends on the following packages:
- Python (3.5.2), Biopython (1.71), NumPy (1.14.5), NetworkX (1.11)
- Mash (2.0) <https://github.com/marbl/Mash>
- MUMmer (4.0.0-beta2) <https://github.com/gmarcais/mummer>
- MeDuSa (1.6) <https://github.com/combogenomics/medusa>

ReferenceSeeker has been tested against aforementioned versions.
