# ReferenceSeeker: Fast determination of reference genomes.
Author: Oliver Schwengers (oliver.schwengers@computatinal.bio.uni-giessen.de)


## Contents
- Description
- Input & Output
- Installation
- Usage
- Examples
- Database
- Dependencies
- Citation


## Description
ReferenceSeeker determines closely related and finished reference genomes from
RefSeq (<https://www.ncbi.nlm.nih.gov/refseq>) following a hierarchical approach
combining a kmer based lookup and ANI calculations.

ReferenceSeeker computes kmer based genome distances between a query genome and
and a database built on finished RefSeq genomes via Mash (Ondov et al. 2016).
Currently, ReferenceSeeker offers bacterial, archeael and viral databases.
For subsequent candidates ReferenceSeeker computes ANI (average nucleotide identity)
values picking genomes meeting community standard thresholds (Goris, Konstantinos et al. 2007)
ranked by ANI and conserved DNA. Additionally, ReferenceSeeker can use MeDuSa
(Bosi, Donati et al. 2015) to build scaffolds based on the 20 closest reference genomes.


## Input & Output
Input:
draft or finished genomes in fasta format

Output:
tab separated to STDOUT comprising the following columns:
- RefSeq ID
- ANI
- conserved DNA
- NCBI Taxonomy ID
- Organism (genus species strain)


## Installation
To setup refseekr just do the followng:
1. clone the latest version of the repository
2. set REFERENCE_SEEKER_HOME environment variable to the repository directory
3. download and extract the database or create it yourself

Example:
```
git clone git@github.com:oschwengers/referenceseekr.git
export REFERENCE_SEEKER_HOME=referenceseekr
wget db...
tar -xzf db.tar.gz
rm db.tar.gz
```

Alternatively, just use the aforementioned Docker image (oschwengers/referenceseekr) in order to ease the setup process.


## Usage
Usage:
    --db          ReferenceSeeker database path
    --genome      Target draft genome
    --scaffolds   Build scaffolds via MeDuSa (Bosi, Donati et al. 2015) based on detected references
    --output      Output fasta file for built scaffolds
    --cpus        Number of CPUs to use (default = all available)
    --help        Show this help message and exit
    --verbose     Print verbose information


## Examples
Simple:
```
refseekr --db <REFERENCE_SEEKER_DB> --genome <GENOME>
```

Expert: creating scaffolds with verbose output using defined # of CPUs:
```
refseekr --db <REFERENCE_SEEKER_DB> --genome <GENOME> --scaffolds --output scaffolds.fasta --verbose --cpus 8
```

With Docker:
```
docker pull oschwengers/refseekr:latest
docker run -v <REFERENCE_SEEKER_DB>:/db -v <DATA_DIR>:/data oschwengers/referenceseekr:latest --genome <GENOME>
```


## Database
ReferenceSeeker depends on custom kmer databases build on NCBI RefSeq as well as
finished reference genomes in fasta format.
These databases can be downloaded HERE: (type, size zipped, size unzipped)
- bacteria, 8.1 Gb, ~27 Gb
- archaea, 184 Mb, 577 Mb
- viral, 490 Mb, 681 Mb

The latest versions can be built using a custom nextflow pipeline.
Valid values for <DB_TYPE> are:
- 'archaea'
- 'bacteria'
- 'viral'

```
export REFERENCE_SEEKER_HOME=<REFERENCE_SEEKER_DIR>
export DB_TYPE=<DB_TYPE>
curl -fsSL get.nextflow.io | bash
nextflow build-db.nf --domain $DB_TYPE
$REFERENCE_SEEKER_HOME/share/mash/mash paste db $DB_TYPE/*.msh
rm -rf work/ .nextflow* $DB_TYPE/*.msh
mv db.msh $DB_TYPE/
```

## Dependencies
ReferenceSeeker depends on the following packages:
- Python (3.5.2) and BioPython (1.66)
- Mash (2.0) <https://github.com/marbl/Mash>
- MUMmer (4.0.0-beta) <https://github.com/gmarcais/mummer>
- MeDuSa (1.6) <https://github.com/combogenomics/medusa>

ReferenceSeeker has been tested against aforementioned versions.


## Citation
Manuscript is submitted for publication.
