#!/usr/bin/env python3

import argparse
import itertools as it
import multiprocessing as mp
import os
import re
import shutil
import subprocess as sp
import sys
import tempfile

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
    parser.add_argument('--scaffolds', '-s', action='store_true', help='Build scaffolds via MeDuSa (Bosi, Donati et al. 2015) based on detected references')
    parser.add_argument('--output', '-o', help='Output fasta file for built scaffolds')
    parser.add_argument('--threads', '-t', action='store', type=int, default=mp.cpu_count(), help='Number of threads to use (default = number of available CPUs)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
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
    if args.output:
        scaffolds_path = os.path.abspath(args.output)
        if os.path.isdir(scaffolds_path):
            sys.exit('ERROR: output is a directory! Please, provide a valid path to a fasta file.')
    else:
        scaffolds_path = cwd_path + '/scaffolds.fna'

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
        if args.scaffolds:
            print('\tscaffold path: ' + scaffolds_path)
    fh_def_null = open(os.devnull, 'w')

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

    # create scaffolds based on draft genome and detected references
    if args.scaffolds:
        if args.verbose:
            print('\ncreate scaffolds...')
        # copy first 20 references to tmp file
        tmp_dir = tempfile.mkdtemp()
        if len(results) > 20:
            results = results[:20]
        for result in results:
            os.symlink(db_path + '/' + result['id'] + '.fna', tmp_dir + '/' + result['id'] + '.fna')
        os.putenv('PATH', REFERENCE_SEEKER_HOME + '/share/mummer/:' + os.getenv('PATH'))
        sp.check_call(
            [
                'java', '-jar', REFERENCE_SEEKER_HOME + '/share/medusa/medusa.jar',
                '-f', tmp_dir,
                '-i', genome_path,  # assembled alignments
                '-o', scaffolds_path,  # new name
                '-v',  # verbose flag
                '-threads', str(args.threads),  # threads
                '-random', '1000',  # use 1000 random rounds to find the best scaffolds
                '-scriptPath', REFERENCE_SEEKER_HOME + '/share/medusa/medusa_scripts'  # path to medusa script directory
            ],
            cwd=cwd_path,
            stdout=fh_def_null,
            stderr=sp.STDOUT
        )
        shutil.rmtree(tmp_dir)

        # parse MeDuSa output and print results
        medusa_result_path = genome_path + '_SUMMARY'
        if os.path.isfile(medusa_result_path) and os.access(medusa_result_path, os.R_OK):
            n50_pre_scaffolding, l50_pre_scaffolding = util.compute_n50(genome_path)
            with open(medusa_result_path, 'r') as medusa_result_file:
                medusa_result = medusa_result_file.read()
                mg = re.search('singletons = (\d+), multi-contig scaffold = (\d+)', medusa_result)
                no_contigs = int(mg.group(1))
                no_scaffolds = int(mg.group(2))
                mg = re.search('from (\d+) initial fragments', medusa_result)
                no_init_contigs = int(mg.group(1))
                n50_post_scaffolding, l50_post_scaffolding = util.compute_n50(scaffolds_path)
                if args.verbose:
                    print(
                        (
                            'pre-scaffolding:\n'
                            '\tcontigs: %d\n'
                            '\tN50: %d\n'
                            '\tL50: %d\n'
                            'post-scaffolding:\n'
                            '\tscaffolds: %d\n'
                            '\tcontigs: %d\n'
                            '\tN50: %d\n'
                            '\tL50: %d'
                        ) %
                        (
                            no_init_contigs,
                            n50_pre_scaffolding,
                            l50_pre_scaffolding,
                            no_scaffolds,
                            no_contigs,
                            n50_post_scaffolding,
                            l50_post_scaffolding
                        )
                    )
            os.remove(medusa_result_path)
        else:
            sys.exit(
                (
                    'ERROR: Contigs could not be scaffolded by MeDuSa!\n'
                    'Maybe selected reference genomes do not provide enough common sequence information.'
                )
            )


if __name__ == '__main__':
    main()
