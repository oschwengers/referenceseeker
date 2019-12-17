
import argparse
import itertools as it
import multiprocessing as mp
import os
import sys
import shutil

from pathlib import Path

import referenceseeker
import referenceseeker.constants as rc
import referenceseeker.mash as mash
import referenceseeker.util as util
import referenceseeker.ani as ani


def main():
    # parse options and arguments
    parser = argparse.ArgumentParser(
        prog='referenceseeker',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Rapid determination of appropriate reference genomes.',
        epilog="Citation:\n%s\n\nGitHub:\nhttps://github.com/oschwengers/referenceseeker" % rc.CITATION
    )
    parser.add_argument('db', metavar='<database>', help='ReferenceSeeker database path')
    parser.add_argument('genome', metavar='<genome>', help='target draft genome in fasta format')
    parser.add_argument('--crg', '-r', action='store', type=int, default=100, help='max number of candidate reference genomes to pass kmer prefilter (default = 100)')
    parser.add_argument('--ani', '-a', action='store', type=float, default=0.95, help='ANI threshold value (default = 0.95)')
    parser.add_argument('--conserved-dna', '-c', action='store', dest='conserved_dna', type=float, default=0.69, help='Conserved DNA threshold value (default = 0.69)')
    parser.add_argument('--unfiltered', '-u', action='store_true', help='set kmer prefilter to extremely conservative values and skip species level ANI cutoffs (ANI >= 0.95 and conserved DNA >= 0.69')
    parser.add_argument('--verbose', '-v', action='store_true', help='print verbose information')
    parser.add_argument('--threads', '-t', action='store', type=int, default=mp.cpu_count(), help='number of threads to use (default = number of available CPUs)')
    parser.add_argument('--version', '-V', action='version', version='%(prog)s ' + referenceseeker.__version__)
    args = parser.parse_args()

    # setup global configuration
    config = util.setup_configuration(args)

    # check parameters
    db_path = Path(args.db)
    if(not os.access(str(db_path), os.R_OK)):
        sys.exit('ERROR: database directory not readable!')
    db_path = db_path.resolve()
    config['db_path'] = db_path

    genome_path = Path(args.genome)
    if(not os.access(str(genome_path), os.R_OK)):
        sys.exit('ERROR: genome file not readable!')
    genome_path = genome_path.resolve()
    config['genome_path'] = genome_path

    # print verbose information
    if(args.verbose):
        print("ReferenceSeeker v%s" % referenceseeker.__version__)
        print('Options, parameters and arguments:')
        print("\tuse bundled binaries: %s" % str(config['bundled-binaries']))
        print("\tdb path: %s" % str(config['db_path']))
        print("\tgenome path: %s" % str(config['genome_path']))
        print("\ttmp path: %s" % str(config['tmp']))
        print("\tunfiltered: %s" % str(config['unfiltered']))
        print("\tANI: %0.2f" % config['ani'])
        print("\tconserved DNA: %0.2f" % config['conserved_dna'])
        print("\t# CRG: %d" % config['crg'])
        print("\t# threads: %d" % config['threads'])

    # calculate genome distances via Mash
    if(args.verbose):
        print('\nEstimate genome distances...')
    mash_output_path = config['tmp'].joinpath('mash.out')
    mash.run_mash(config, mash_output_path)

    # extract hits and store dist
    accession_ids, mash_distances = mash.parse_mash_results(config, mash_output_path)
    if(args.verbose):
        print("\tscreened %d potential reference genome(s)" % len(accession_ids))

    # reduce Mash output to best hits (args.crg)
    if(len(accession_ids) > args.crg):
        if(args.verbose):
            print("\treduce to best %d hits..." % args.crg)
        tmp_accession_ids = sorted(accession_ids, key=lambda k: mash_distances[k])
        accession_ids = tmp_accession_ids[:args.crg]

    # get assemblies from RefSeq by accessions
    ref_genomes = util.read_reference_genomes(config, accession_ids, mash_distances)

    # build dna fragments
    dna_fragments_path = config['tmp'].joinpath('dna-fragments.fasta')
    dna_fragments = util.build_dna_fragments(config, dna_fragments_path)

    # copy genomes, extract them and build ANI
    if(args.verbose):
        print('\nCompute ANIs...')
    pool = mp.Pool(args.threads)
    results = pool.starmap(ani.compute_ani, zip(it.repeat(config), it.repeat(dna_fragments_path), it.repeat(dna_fragments), ref_genomes))
    pool.close()
    pool.join()

    # remove tmp dir
    shutil.rmtree(str(config['tmp']))

    # filter and sort results
    tmp_results = []
    for result in results:
        if(args.unfiltered or ((result['conserved_dna'] >= config['conserved_dna']) and (result['ani'] >= config['ani']))):
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
