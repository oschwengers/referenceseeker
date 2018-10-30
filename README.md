# ReferenceSeeker: Fast determination of reference genomes.
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
RefSeq (<https://www.ncbi.nlm.nih.gov/refseq>) following a hierarchical approach
combining a kmer based lookup and subsequent ANI calculations.

ReferenceSeeker computes kmer based genome distances between a query genome and
and a database built on RefSeq genomes via Mash (Ondov et al. 2016). Hereby, only
complete genomes or those stated as 'representative' or 'reference' genome are included.
Currently, ReferenceSeeker offers bacterial, archeael, fungi and viral databases.
For subsequent candidates ReferenceSeeker computes ANI (average nucleotide identity)
values picking genomes meeting community standard thresholds (ANI >= 95 % & conserved DNA >= 69 %)
(Goris, Konstantinos et al. 2007) ranked by ANI and conserved DNA.
Additionally, ReferenceSeeker can use MeDuSa (Bosi, Donati et al. 2015)
to scaffold contigs based on the 20 closest reference genomes.


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
git clone git@github.com:oschwengers/referenceseekr.git
export REFERENCE_SEEKER_HOME=~/referenceseekr
wget https://s3.computational.bio.uni-giessen.de/swift/v1/referenceseeker/bacteria.tar.gz
tar -xzf bacteria.tar.gz
rm bacteria.tar.gz
```

Alternatively, just use the aforementioned Docker image (oschwengers/referenceseekr) in order to ease the setup process.


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
referenceseekr.py --db <REFERENCE_SEEKER_DB> <GENOME>
```

Expert: creating scaffolds with verbose output using a defined number of threads:
```
referenceseekr.py --db <REFERENCE_SEEKER_DB> --scaffolds --output scaffolds.fasta --verbose --threads 8 <GENOME>
```

With Docker:
```
docker pull oschwengers/referenceseekr:latest
docker run --rm -v <REFERENCE_SEEKER_DB>:/db -v <DATA_DIR>:/data oschwengers/referenceseekr:latest <GENOME>
```

With Docker shell script:
```
docker pull oschwengers/referenceseekr:latest
referenceseekr.sh <REFERENCE_SEEKER_DB> <GENOME>
```


## Databases
ReferenceSeeker depends on custom databases based on reference, representative as well as complete NCBI RefSeq genomes
comprising kmer hash subsets as well as fasta files.
These databases (RefSeq release 87) can be downloaded from the the following list: (type, link, # genomes, size zipped, size unzipped)
- bacteria: <https://s3.computational.bio.uni-giessen.de/swift/v1/referenceseeker/bacteria.tar.gz>, 13563, 17 Gb, 53 Gb
- archaea: <https://s3.computational.bio.uni-giessen.de/swift/v1/referenceseeker/archaea.tar.gz>, 374, 324 Mb, 1.1 Gb
- viral: <https://s3.computational.bio.uni-giessen.de/swift/v1/referenceseeker/viral.tar.gz>, 7532, 508 Mb, 761 Mb
- fungi: <https://s3.computational.bio.uni-giessen.de/swift/v1/referenceseeker/fungi.tar.gz>, 250, 2.2 Gb, 6.8 Gb

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
sh build-db.sh <DB_TYPE>
```

## Dependencies
ReferenceSeeker depends on the following packages:
- Python (3.5.2) and BioPython (1.66)
- Mash (2.0) <https://github.com/marbl/Mash>
- MUMmer (4.0.0-beta2) <https://github.com/gmarcais/mummer>
- MeDuSa (1.6) <https://github.com/combogenomics/medusa>

ReferenceSeeker has been tested against aforementioned versions.
