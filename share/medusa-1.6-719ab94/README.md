Medusa
======

A draft genome scaffolder that uses multiple reference genomes in a
graph-based approach.

Availability and dependencies
-----------------------------

The present document provides a short guide for using the stand-alone
version of the software Medusa. This software has not yet been published.
A web interface is available at <http://combo.dbe.unifi.it/medusa>.
The source code, precompiled version and the present manual are
accessible at <https://github.com/combogenomics/medusa>.

Medusa depends on the following packages being installed on your system
and available in your PATH:

1.  MUMmer: this software is available at
    <http://mummer.sourceforge.net/>.

2.  Python (from 2.6) and BioPython (from 1.61).

3.  Java (from 1.6).

The archive Medusa.tar.gz contains the following files:

1.  A runnable .jar file medusa.jar This is the program you will run.

2.  A sub-folder with python scripts needed to run the program (medusa_scripts). Leave it in the same folder of
    the .jar file.

3.  A sub-folder with a dataset (test) that can be used to test the tool.

4.  A sub-folder with scripts useful for benchmarking the tool.

Input and Output
----------------

The following inputs are required:

-   The *targetGenome* file: a draft genome in fasta format. This is the
    genome you are interested in scaffolding.

-   An arbitrary long list of *auxiliaryDraft* files: other draft
    genomes in fasta format. The closest these organisms are related to
    the target, the better the results will be. These files are expected
    to be collected in a specific directory. It is possible to specify
    the path to the directory, see the command "-f" in the next section.

The following output files will be produced.

-   targetGenome_SUMMARY: a textual file containing information about
    your data. Number of scaffolds, N50 value etc..

-   targetGenomeScaffold.fasta: a fasta file with the sequences grouped
    in scaffolds. Contigs in the same scaffolds are separated by 100 Ns
    by default, or a variable number of Ns (estimate of the distance between
    the contigs), if the option "-d" is used. 
    
The following output files can optionally be produced.

-   targetGenome_distanceTable: a tabular file with the estimation of the
	distance between successive contigs (bp).
	
-   targetGenome_network.gexf: the contig network in gexf format.

-   targetGenome_cover.gexf: the final path cover in gexf format.


Usage
-----

The project folder must contain:

-   the *targetGenome* in fasta format.

-   the medusa.jar file

-   the scripts sub-folder “medusa_scripts”.

-   the comparison genomes sub-folder “drafts”. (In alternative you can
    specify another path for this folder usinf the "-f" option)

Medusa can be run with the following parameters:

1.  The option *-i* is required and indicates the name of the target
    genome file.

2.  The option *-o* is optional and indicates the name of output fasta
    file.

3.  The option *-v* (recommended) print on console the information given
    by the package MUMmer. This option is strongly suggested to
    understand if MUMmer is not running properly.

4.  The option *-f* is optional and indicates the path to the comparison
    drafts folder.

5.  The option *-random* is available (not required). This option allows
    the user to run a given number of cleaning rounds and keep the best
    solution. Since the variability is small, 5 rounds are
    usually sufficient to find the best score.

6.  The option *-w2* is optional and allows for a sequence similarity
    based weighting scheme. Using a different weighting scheme may lead
    to better results.

7. The option *-d* allows for the estimation of the distance between pairs of contigs based on the reference genome(s):
	in this case the scaffolded contigs will be separated by a number of N characters equal to this estimate.
	The estimated distances are also saved in the "*_distanceTable" file.
	By default the scaffolded contigs are separated by 100 Ns.
	
8. The *-gexf* is optional. With this option the gexf format of the contig network and
	the path cover are porvided.

9.  The option *-n50* allows the calculation of the N50 statistic on a FASTA file. 
	In this case the usage is the following: java -jar medusa.jar -n50 <name_of_the_fasta>
	All the other options will be ignored.

10. Finally the *-h* option provides a small recap of the previous ones.

The Medusa archive
------------------

When *medusa* archive is unzipped the following files will be extracted:

-   the medusa.jar file.

-   the scripts sub-folder “medusa_scripts”.

-   the utility test scripts folder "medusa_testing"

-   a folder “test”, containing one test bacterial datasets.

Running an example
-------------------

    java -jar medusa.jar -f test/reference_genomes/ -i test/Rhodobacter_target.fna -v

Additional datasets for benchmarking
------------------------------------

Additional datasets can be retrieved at the medusa_datasets repository <https://github.com/combogenomics/medusa_datasets>.

Just type

    git clone https://github.com/combogenomics/medusa_datasets.git

Compile
-------

The project can be compiled by calling ant in the top-level directory:

    ant
