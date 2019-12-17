
import argparse
import itertools as it
import multiprocessing as mp
import os
import subprocess as sp
import sys

from pathlib import Path

import referenceseeker
import referenceseeker.constants as rc
import referenceseeker.util as util
import referenceseeker.ani as ani


def main():
    # parse options and arguments
    parser = argparse.ArgumentParser(
        prog='ReferenceSeeker',
        description='Rapid determination of appropriate reference genomes.'
    )
    parser.add_argument('genome', metavar='<genome>', help='Target draft genome in fasta format')
    parser.add_argument('--db', '-d', required=True, help='ReferenceSeeker database path')
    parser.add_argument('--crg', '-c', action='store', type=int, default=100, help='Max number of candidate reference genomes to assess (default = 100)')
    parser.add_argument('--unfiltered', '-u', action='store_true', help='Set kmer prefilter to extremely conservative values and skip species level ANI cutoffs (ANI >= 0.95 and conserved DNA >= 0.69')
    parser.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    parser.add_argument('--threads', '-t', action='store', type=int, default=mp.cpu_count(), help='Number of threads to use (default = number of available CPUs)')
    parser.add_argument('--version', '-V', action='version', version='%(prog)s ' + referenceseeker.__version__)
    parser.add_argument('--citation', '-c', action='store_true', help='Print citation')
    args = parser.parse_args()

    # print citation
    if(args.citation):
        print(rc.CITATION)
        sys.exit()

    # setup global configuration
    config = util.setup_configuration()

    # check parameters
    db_path = Path(args.db)
    if(not os.access(str(db_path), os.R_OK)):
        sys.exit('ERROR: database directory not readable!')
    db_path = db_path.resolve()

    genome_path = Path(args.genome)
    if(not os.access(str(genome_path), os.R_OK)):
        sys.exit('ERROR: genome file not readable!')
    genome_path = genome_path.resolve()

    # print verbose information
    if(args.verbose):
        print("ReferenceSeeker v%s" % referenceseeker.__version__)
        print('Options, parameters and arguments:')
        print("\tuse bundled binaries: %s" % str(config['bundled-binaries']))
        print("\tdb path: %s" % str(config['db']))
        print("\tgenome path: %s" % str(genome_path))
        print("\ttmp path: %s" % str(config['tmp']))
        print("\tunfiltered: %s" % str(args.unfiltered))
        print("\t# CRG: %d" % args.crg)
        print("\t# threads: %d" % args.threads)

    # calculate genome distances via Mash
    if(args.verbose):
        print('\nAssess genome distances...')
    mash_result_path = config['tmp'].joinpath('/mash.out')
    with open(mash_result_path, 'w') as fh:
        sp.check_call(
            [
                'mash',
                'dist',
                '-d', rc.UNFILTERED_MASH_DIST if args.unfiltered else rc.MAX_MASH_DIST,
                str(db_path.joinpath('/db.msh')),
                genome_path
            ],
            stdout=fh,
            stderr=sp.STDOUT
        )

    # extract hits and store dist
    accession_ids = []
    mash_distances = {}
    with open(mash_result_path, 'r') as fh:
        for line in fh:
            cols = line.rstrip().split()
            accession_ids.append(cols[0])
            mash_distances[cols[0]] = float(cols[2])
    os.remove(mash_result_path)
    if(args.verbose):
        print("\tscreened %d potential reference genome(s)" % len(accession_ids))

    # reduce Mash output to best args.crg hits
    if len(accession_ids) > args.crg:
        if(args.verbose):
            print("\treduce to best %d hits..." % args.crg)
        tmp_accession_ids = sorted(accession_ids, key=lambda k: mash_distances[k])
        accession_ids = tmp_accession_ids[:args.crg]

    # get assemblies from RefSeq by accessions
    ref_genomes = util.read_reference_genomes(db_path, accession_ids, mash_distances)

    # build dna fragments
    dna_fragments_path = cwd_path + '/dna-fragments.fasta'
    dna_fragments = util.build_dna_fragments(genome_path, dna_fragments_path)

    # copy genomes, extract them and build ANI
    if(args.verbose):
        print('\nCompute ANIs...')
    pool = mp.Pool(args.threads)
    results = pool.starmap(ani.compute_ani, zip(it.repeat(db_path), it.repeat(dna_fragments_path), it.repeat(dna_fragments), ref_genomes))
    pool.close()
    pool.join()
    os.remove(dna_fragments_path)

    # filter and sort results
    tmp_results = []
    for result in results:
        if args.unfiltered or ((result['conserved_dna'] >= 0.69) and (result['ani'] >= 0.95)):
            tmp_results.append(result)
    results = sorted(tmp_results, key=lambda k: -(k['ani'] * k['conserved_dna']))

    # print results to STDOUT
    if(args.verbose):
        print('')
    print('#ID\tANI\tCon. DNA\tMash Distance\tTaxonomy ID\tAssembly Status\tOrganism')
    for result in results:
        print(
            '%s\t%2.2f\t%2.2f\t%1.5f\t%s\t%s\t%s' %
            (
                result['id'],
                result['ani'] * 100,
                result['conserved_dna'] * 100,
                result['mash_dist'],
                result['tax'],
                result['status'],
                result['name']
            )
        )


if __name__ == '__main__':
    main()
