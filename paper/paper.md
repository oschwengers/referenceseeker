---
title: 'ReferenceSeeker: rapid determination of appropriate reference genomes'
tags:
  - Python
  - Bioinformatics
  - WGS
  - NGS
  - Microbiology
authors:
  - name: Oliver Schwengers
    orcid: 0000-0003-4216-2721
    affiliation: "1, 2, 3"
  - name: Torsten Hain
    affiliation: "2, 3"
  - name: Trinad Chakraborty
    affiliation: "2, 3"
  - name: Alexander Goesmann
    orcid: 0000-0002-7086-2568
    affiliation: "1, 3"
affiliations:
 - name: Bioinformatics and Systems Biology, Justus Liebig University Giessen, Giessen, 35392, Germany
   index: 1
 - name: Institute of Medical Microbiology, Justus Liebig University Giessen, Giessen, 35392, Germany
   index: 2
 - name: German Centre for Infection Research (DZIF), partner site Giessen-Marburg-Langen, Giessen, Germany
   index: 3
date: 18 December 2019
bibliography: paper.bib
---

# Summary
The enormous success and ubiquitous application of next and third generation sequencing has
led to a large number of available high-quality draft and complete microbial genomes in the
public databases. Today, the NCBI RefSeq database release 90 alone contains 11,060
complete bacterial genomes ​[@Haft:2018​]. Concurrently, selection of appropriate reference
genomes (RGs) is increasingly important as it has enormous implications for routine in-silico
analyses, as for example in detection of single nucleotide polymorphisms, scaffolding of draft
assemblies, comparative genomics and metagenomic tasks. Therefore, a rigorously selected
RG is a prerequisite for the accurate and successful application of the aforementioned
bioinformatic analyses. In order to address this issue several new databases, methods and tools
have been published in recent years *e.g.* RefSeq, DNA-DNA hybridization [@Meier-Kolthoff:2013]​,
average nucleotide identity (ANI) values [@Goris:2007] and Mash ​[@Ondov:2016]​.
Nevertheless, the sheer amount of currently available databases and potential RGs
contained therein, together with the plethora of tools available, often require manual selection of
the most suitable RGs. To the best of the authors’ knowledge, there is currently no such tool
providing both an integrated, highly specific workflow and scalable and rapid implementation.
ReferenceSeeker was designed to overcome this bottleneck. As a novel command line tool, it
combines a fast kmer profile-based lookup of candidate reference genomes (CRGs) from high
quality databases with rapid computation of (mutual) highly specific ANI and conserved DNA values.

# Implementation
ReferenceSeeker is a linux command line tool implemented in Python 3. All necessary external
binaries are bundled with the software. The tool Itself requires no external dependencies other
than Biopython for file input and output.

## Databases
ReferenceSeeker takes advantage of taxon-specific custom databases in order to reduce data
size and overall runtime. Pre-built databases for the taxonomic groups bacteria, archaea, fungi,
protozoa and viruses are provided. Each database integrates genomic as well as taxonomic
information comprising genome sequences of all RefSeq genomes with an assembly level
‘complete’ or whose RefSeq category is either denoted as ‘reference genome’ or ‘representative
genome’, as well as kmer profiles, related species names, NCBI Taxonomy identifiers and
RefSeq assembly identifiers. For convenient and fully automatic updates, we provide locally
executable scripts implemented in bash and Nextflow ​[@Di_Tommaso:2017]​ .

## Database Lookup of CRGs
To reduce the number of necessary ANI calculations a kmer profile-based lookup of CRGs
against custom databases is carried out. This step is implemented via Mash parameterized with
a Mash distance of 0.1, which was shown to correlate well with an ANI of roughly 90% ​
[@Ondov:2016] and thereby establishing a lower limit for reasonably related genomes.
The resulting set of CRGs is subsequently reduced to a configurable number of CRGs (default=100)
with the lowest Mash distances.

## Determination of RG
As a highly specific measure for microbial genome relationships ReferenceSeeker uses ANI and
conserved DNA (conDNA) values. The reasoning for the subsequent calculation of ANI and
conDNA values is that although Mash distance values correlate well with ANI values, the same
cannot be observed for conDNA values (Figure 1).

Calculation of these metrics is implemented vía Nucmer contained in the MUMmer package
[@Marcais:2018] as it was recently shown that Nucmer based implementations (ANIn)
compare favourably against BLAST+ based implementations (ANIb) in terms of runtime. Given
that compared genomes are closely related, i.e. shared ANI is above 90%, it was also shown
that ANIn correlates well with ANIb [@Yoon:2017]. This is ensured by the prior Mash-based
selection of CRGs. As an established threshold for species boundaries [@Goris:2007]​,
results are subsequently filtered by configurable ANI and conDNA values with a default of 95%
and 69%, respectively. Finally, CRGs are sorted according to the harmonic mean of ANI and
conDNA values in order to incorporate both nucleotide identity and genome coverage of query
genomes and resulting CRGs. In this manner, ReferenceSeeker ensures that the resulting RGs
sufficiently reflect the genomic landscape of a query genome.

![Figure 1: Scatter plots showing the correlation between Mash distance, ANI and conDNA
values. ANI and conserved DNA values are plotted against Mash distance values for 500
candidate reference genomes with the lowest Mash distance within the bacterial database for
10 randomly selected *Escherichia coli* genomes from the RefSeq database, each.](mash-ani-cdna.scatter.png)

# Application
ReferenceSeeker takes as input a microbial genome assembly in fasta format and the path to a
taxonomic database of choice. Results are returned as a tabular separated list comprising the
following information: RefSeq assembly identifier, ANI, conDNA, NCBI taxonomy identifier,
assembly status and organism name.
To illustrate the broad applicability at different scales we tested ReferenceSeeker with 12
microbial genomes from different taxonomic groups and measured overall runtimes on a
common consumer laptop providing 4 cores and a server providing 64 cores (Table \ref{Table 1}). For the
tested bacterial genomes, ReferenceSeeker limited the number of resulting RGs to a default
maximum of 100 genomes. Runtimes of archaeal and viral genomes are significantly shorter
due to a small number of available RGs in the database and overall smaller genome sizes,
respectively.

Table: Runtimes and numbers of resulting RG executed locally on a quad-core moderate consumer
laptop and a 64 core server machine. \label{Table 1}

------------------------------------------------------------------------------------
Genome                                           Genome Size   Laptop   Server  # RG
                                                        [kb]  [mm:ss]  [mm:ss]      
----------------------------------------------- ------------ -------- -------- -----
*Escherichia coli* str. K-12 substr. MG1665            4,641     3:24     0:30  100*
(GCF_000005845.2)                                                                   

*Pseudomonas aeruginosa* PAO1                          6,264     5:20     0:44  100*
(GCF_000006765.1)                                                                   

*Listeria monocytogenes* ​EGD-e                         2,944     2:52     0:24  100*
(GCF_000196035.1)                                                                   

*Staphylococcus aureus* subsp aureus NCTC 8325         2,821     2:31     0:21  100*
(GCF_000013425.1)                                                                   

*Halobacterium salinarum* NRC-1                        2,571     0:04     0:03     2
(GCF_000006805.1)                                                                   

*Methanococcus maripaludis* X1                         1,746     0:22     0:09     5
(GCF_000220645.1)                                                                   

*Aspergillus fumigatus* ​Af293                         29,384     3:11     2:07     1
(GCF_000002655.1)                                                                   

*Candida albicans* SC5314                             14,282     0:21     0:19     1
(GCF_000182965.3)                                                                   

*Entamoeba histolytica* HM-1:IMSS                     20,835     6:04     4:41     1
(GCF_000208925.1)                                                                   

*Plasmodium falciparum* ​3D7                           23,326     2:52     1:49     1
(GCF_000002765.4)                                                                   

*Influenza A* virus                                       13     0:03     0:02     1
(GCF_001343785.1)                                                                   

*Human coronavirus* NL63                                  27     0:04     0:02     1
(GCF_000853865.1)                                                                   
------------------------------------------------------------------------------------

# Availability
The source code is available on GitHub under a GPL3 license: https://github.com/oschwengers/referenceseeker​.
The software is packaged and publicly available via BioConda. Pre-built databases for
bacteria, archaea, fungi, protozoa and viruses are hosted at Zenodo: https://doi.org/10.5281/zenodo.3562005​.
All installation instructions, examples and download links are provided on GitHub.

# Funding
This work was supported by the German Center of Infection Research (DZIF) [DZIF grant 8000
 701–3 (HZI), TI06.001 and 8032808811 to T.C.]; the German Network for Bioinformatics
Infrastructure (de.NBI) [BMBF grant FKZ 031A533B to A.G.]; and the German Research
Foundation (DFG) [SFB-TR84 project A04 (TRR84/3 2018) to T.C., KFO309 Z1 (GO 2037/5-1)
to A.G., SFB-TR84 project B08 (TRR84/3 2018) to T.H., SFB1021 Z02 (SFB 1021/2 2017) to
T.H., KFO309 Z1 (HA 5225/1-1) to T.H.].

Authors declare that there are no conflicts of interest.

# Acknowledgement
The authors thank Karina Brinkrolf for valuable discussions, testing and bug reports.

# References
