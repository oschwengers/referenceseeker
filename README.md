# ReferenceSeeker: fast determination of suitable reference genomes.
Author: Oliver Schwengers (oliver.schwengers@computational.bio.uni-giessen.de)


## Contents
- Description
- Input & Output
- Installation
- Usage
- Examples
- Databases
- Dependencies


## Description
ReferenceSeeker determines closely related reference genomes from
RefSeq (<https://www.ncbi.nlm.nih.gov/refseq>) following a scalable hierarchical
approach combining a fast kmer profile based database lookup of candidates and
subsequent calculation of specific average nucleotide identity (ANI) values
for the fast determination of suitable reference genomes.

ReferenceSeeker computes kmer profile based genome distances between a query genome and
and a database built on RefSeq genomes via Mash (Ondov et al. 2016). Hereby, only
complete genomes or those stated as 'representative' or 'reference' genome are included.
Currently, ReferenceSeeker offers pre-built bacterial, archeael, fungi and viral
databases. For resulting candidate reference genomes ReferenceSeeker subsequently
computes ANI and conserved DNA values filtered to community standard thresholds (ANI >= 95 % &
conserved DNA >= 69 %) (Goris, Konstantinos et al. 2007) ranked by the harmonic
mean of ANI and conserved DNA. Optionally, for draft assembly inputs
ReferenceSeeker can use MeDuSa (Bosi, Donati et al. 2015) to scaffold contigs
based on the 20 closest reference genomes.


## Input & Output
Input:
draft or finished genomes in fasta format

Output:
tab separated to STDOUT comprising the following columns:
- RefSeq ID
- ANI
- conserved DNA
- NCBI Taxonomy ID
- Assembly Status
- Organism (genus species strain)


## Installation
To setup ReferenceSeeker just do the following:
1. clone the latest version of the repository
2. set REFERENCE_SEEKER_HOME environment variable pointing to the repository directory
3. download and extract a databases or create one yourself

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

Alternatively, use the aforementioned Docker image (oschwengers/referenceseekr)
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
referenceseeker.py --db <REFERENCE_SEEKER_DB> <GENOME>
```

Expert: creating scaffolds with verbose output using a defined number of threads:
```
referenceseeker.py --db <REFERENCE_SEEKER_DB> --scaffolds --output scaffolds.fasta --verbose --threads 8 <GENOME>
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
ReferenceSeeker depends on custom databases based on reference, representative as well as complete NCBI RefSeq genomes
comprising kmer hash subsets as well as fasta files.
These databases (RefSeq release 90) can be downloaded from the the following list: (type, link, # genomes, size zipped, size unzipped)
- bacteria: <https://s3.computational.bio.uni-giessen.de/swift/v1/referenceseeker/bacteria.tar.gz>, 14,983, 18 Gb, 58 Gb
- archaea: <https://s3.computational.bio.uni-giessen.de/swift/v1/referenceseeker/archaea.tar.gz>, 386, 335 Mb, 1.1 Gb
- viral: <https://s3.computational.bio.uni-giessen.de/swift/v1/referenceseeker/viral.tar.gz>, 7,855, 525 Mb, 719 Mb
- fungi: <https://s3.computational.bio.uni-giessen.de/swift/v1/referenceseeker/fungi.tar.gz>, 277, 2.5 Gb, 7.7 Gb

The latest versions can be built using a custom nextflow pipeline.
Valid values for `DB_TYPE` are:
- 'archaea'
- 'bacteria'
- 'viral'
- 'fungi'

Download and install Nextflow:
```
curl -fsSL get.nextflow.io | bash
```

Build database:
```
export REFERENCE_SEEKER_HOME=<REFERENCE_SEEKER_DIR>
sh build-db.sh <DB_TYPE_OPTION>
```

## Dependencies
ReferenceSeeker depends on the following packages:
- Python (3.5.2), Biopython (1.71), NumPy (1.14.5), NetworkX (1.11)
- Mash (2.0) <https://github.com/marbl/Mash>
- MUMmer (4.0.0-beta2) <https://github.com/gmarcais/mummer>
- MeDuSa (1.6) <https://github.com/combogenomics/medusa>

ReferenceSeeker has been tested against aforementioned versions.
