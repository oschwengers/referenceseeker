#!/usr/bin/env python3

import argparse
import itertools as it
import multiprocessing as mp
import os
import subprocess as sp
import sys

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
    parser.add_argument('--threads', '-t', action='store', type=int, default=mp.cpu_count(), help='Number of threads to use (default = number of available CPUs)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    parser.add_argument('--version', action='version', version='%(prog)s ' + referenceseeker.__version__)
    args = parser.parse_args()

    # check parameters & environment variables
    REFERENCE_SEEKER_HOME = os.getenv('REFERENCE_SEEKER_HOME', None)
    if REFERENCE_SEEKER_HOME is None:
        sys.exit('ERROR: REFERENCE_SEEKER_HOME not set!')
    REFERENCE_SEEKER_HOME = os.path.abspath(REFERENCE_SEEKER_HOME)
    if not os.access(REFERENCE_SEEKER_HOME, os.R_OK & os.X_OK):
        sys.exit('ERROR: REFERENCE_SEEKER_HOME (' + REFERENCE_SEEKER_HOME + ') not readable/accessible!')

    db_path = os.path.abspath(args.db)
    if not os.access(db_path, os.R_OK):
        sys.exit('ERROR: database directory not readable!')

    genome_path = os.path.abspath(args.genome)
    if not os.access(genome_path, os.R_OK):
        sys.exit('ERROR: genome file not readable!')

    cwd_path = os.path.abspath(os.getcwd())

    # print verbose information
    if args.verbose:
        print('Options, parameters and arguments:')
        print('\tREFERENCE_SEEKER_HOME: ' + REFERENCE_SEEKER_HOME)
        print('\tdb path: ' + db_path)
        print('\tgenome path: ' + genome_path)
        print('\t# CRG: ' + str(args.crg))
        print('\tunfiltered: ' + str(args.unfiltered))
        print('\tbuild scaffolds: ' + str(args.scaffolds))
        print('\t# threads: ' + str(args.threads))

    # calculate genome distances via Mash
    if args.verbose:
        print('\nAssess genome distances...')
    mash_result_path = cwd_path + '/mash.out'
    with open(mash_result_path, 'w') as fh:
        sp.check_call(
            [
                REFERENCE_SEEKER_HOME + '/share/mash/mash',
                'dist',
                '-d', rc.UNFILTERED_MASH_DIST if args.unfiltered else rc.MAX_MASH_DIST,
                db_path + '/db.msh',
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
    if args.verbose:
        print('\tscreened ' + str(len(accession_ids)) + ' potential reference genome(s)')

    # reduce Mash output to best args.crg hits
    if len(accession_ids) > args.crg:
        if args.verbose:
            print('\treduce to best ' + str(args.crg) + ' hits...')
        tmp_accession_ids = sorted(accession_ids, key=lambda k: mash_distances[k])
        accession_ids = tmp_accession_ids[:args.crg]

    # get assemblies from RefSeq by accessions
    ref_genomes = util.read_reference_genomes(db_path, accession_ids, mash_distances)

    # build dna fragments
    dna_fragments_path = cwd_path + '/dna-fragments.fasta'
    dna_fragments = util.build_dna_fragments(genome_path, dna_fragments_path)

    # copy genomes, extract them and build ANI
    if args.verbose:
        print('\nCompute ANIs...\n\tID\tANI\tConserved DNA\tMash Distance')
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
    if args.verbose:
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
