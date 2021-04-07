
import argparse
import multiprocessing as mp
import os
import sys
import shutil

import concurrent.futures as cf
from pathlib import Path

from xopen import xopen

import referenceseeker
import referenceseeker.constants as rc
import referenceseeker.mash as mash
import referenceseeker.util as util
import referenceseeker.ani as rani


def main():
    # parse options and arguments
    parser = argparse.ArgumentParser(
        prog='referenceseeker',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Rapid determination of appropriate reference genomes.',
        epilog=f'Citation:\n{rc.CITATION}\n\nGitHub:\nhttps://github.com/oschwengers/referenceseeker',
        add_help=False
    )
    parser.add_argument('db', metavar='<database>', help='ReferenceSeeker database path')
    parser.add_argument('genome', metavar='<genome>', help='target draft genome in fasta format')
    group_workflow = parser.add_argument_group('Filter options / thresholds', 'These options control the filtering and alignment workflow.')
    group_workflow.add_argument('--crg', '-r', action='store', type=int, default=100, help='Max number of candidate reference genomes to pass kmer prefilter (default = 100)')
    group_workflow.add_argument('--ani', '-a', action='store', type=float, default=0.95, help='ANI threshold (default = 0.95)')
    group_workflow.add_argument('--conserved-dna', '-c', action='store', dest='conserved_dna', type=float, default=0.69, help='Conserved DNA threshold (default = 0.69)')
    group_workflow.add_argument('--unfiltered', '-u', action='store_true', help='Set kmer prefilter to extremely conservative values and skip species level ANI cutoffs (ANI >= 0.95 and conserved DNA >= 0.69')
    group_workflow.add_argument('--bidirectional', '-b', action='store_true', help='Compute bidirectional ANI/conserved DNA values (default = False)')

    group_runtime = parser.add_argument_group('Runtime & auxiliary options')
    group_runtime.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    group_runtime.add_argument('--version', '-V', action='version', version=f'%(prog)s {referenceseeker.__version__}')
    group_runtime.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    group_runtime.add_argument('--threads', '-t', action='store', type=int, default=mp.cpu_count(), help='Number of used threads (default = number of available CPU cores)')
    args = parser.parse_args()

    # setup global configuration
    config = util.setup_configuration(args)
    util.test_binaries(config)

    # check parameters
    db_path = Path(args.db)
    if(not os.access(str(db_path), os.R_OK)):
        sys.exit('ERROR: database directory not readable!')
    db_path = db_path.resolve()
    config['db_path'] = db_path

    genome_path = Path(args.genome)
    if(not os.access(str(genome_path), os.R_OK)):
        sys.exit('ERROR: genome file not readable!')
    if(genome_path.stat().st_size == 0):
        sys.exit(f'ERROR: genome file ({genome_path}) is empty!')
    genome_path = genome_path.resolve()
    config['genome_path'] = genome_path

    # print verbose information
    if(args.verbose):
        print(f'ReferenceSeeker v{referenceseeker.__version__}')
        print('Options, parameters and arguments:')
        print(f"\tdb path: {config['db_path']}")
        print(f"\tgenome path: {config['genome_path']}")
        print(f"\ttmp path: {config['tmp']}")
        print(f"\tunfiltered: {config['unfiltered']}")
        print(f"\tbidirectional: {config['bidirectional']}")
        print(f"\tANI: {config['ani']:0.2f}")
        print(f"\tconserved DNA: {config['conserved_dna']:0.2f}")
        print(f"\t# CRG: {config['crg']}")
        print(f"\t# threads: {config['threads']}")
    
    # import potentially gzipped query genome
    if(genome_path.suffix == '.gz'):
        genome_unzipped_path = config['tmp'].joinpath(f'{genome_path.stem}.fna')
        with genome_unzipped_path.open(mode='w') as fh_out, xopen(str(genome_path), threads=0) as fh_in:
            for line in fh_in:
                fh_out.write(line)
        genome_path = genome_unzipped_path
        config['genome_path'] = genome_path

    # calculate genome distances via Mash
    if(args.verbose):
        print('\nEstimate genome distances...')
    mash_output_path = config['tmp'].joinpath('mash.out')
    mash.run_mash(config, mash_output_path)

    # extract hits and store dist
    screened_ref_genome_ids, mash_distances = mash.parse_mash_results(config, mash_output_path)
    if(args.verbose):
        print(f'\tscreened {len(screened_ref_genome_ids)} potential reference genome(s)')

    # reduce Mash output to best hits (args.crg)
    if(len(screened_ref_genome_ids) > args.crg):
        if(args.verbose):
            print(f'\treduce to best {args.crg} hits...')
        tmp_screened_ref_genome_ids = sorted(screened_ref_genome_ids, key=lambda k: mash_distances[k])
        screened_ref_genome_ids = tmp_screened_ref_genome_ids[:args.crg]

    # get genomes from RefSeq by accessions
    ref_genomes = util.read_reference_genomes(config)
    screened_ref_genomes = {k: v for k, v in ref_genomes.items() if k in screened_ref_genome_ids}

    # build dna fragments
    dna_fragments_path = config['tmp'].joinpath('dna-fragments.fasta')
    dna_fragments = util.build_dna_fragments(genome_path, dna_fragments_path)

    # align query fragments to reference genomes and compute ANI/conserved DNA
    results = {}
    if(args.verbose):
        print('\nCompute ANIs...')
    with cf.ThreadPoolExecutor(max_workers=args.threads) as tpe:
        futures = []
        for id, ref_genome in screened_ref_genomes.items():
            futures.append(tpe.submit(rani.align_query_genome, config, dna_fragments_path, dna_fragments, id))
        for f in futures:
            ref_genome_id, ani, conserved_dna = f.result()
            results[ref_genome_id] = [(ani, conserved_dna)]

    # align reference genomes fragments to query genome and compute ANI/conserved DNA
    if(args.bidirectional):
        if(args.verbose):
            print('\nCompute reverse ANIs...')
        with cf.ProcessPoolExecutor(args.threads) as ppe:
            futures = []
            for id, ref_genome in screened_ref_genomes.items():
                futures.append(ppe.submit(rani.align_reference_genome, config, genome_path, id))
            for f in futures:
                ref_genome_id, ani, conserved_dna = f.result()
                result = results[ref_genome_id]
                result.append((ani, conserved_dna))

    # remove tmp dir
    shutil.rmtree(str(config['tmp']))

    # filter and sort results
    filtered_reference_ids = []
    for ref_genome_id, result in results.items():
        if(args.unfiltered):
            filtered_reference_ids.append(ref_genome_id)
        else:
            if(args.bidirectional):
                query_ref = result[0]
                ref_query = result[1]
                if((query_ref[0] >= config['ani']) and (query_ref[1] >= config['conserved_dna'])
                        and (ref_query[0] >= config['ani']) and (ref_query[1] >= config['conserved_dna'])):
                    filtered_reference_ids.append(ref_genome_id)
            else:
                (ani, conserved_dna) = result[0]
                if((conserved_dna >= config['conserved_dna']) and (ani >= config['ani'])):
                    filtered_reference_ids.append(ref_genome_id)

    # sort and print results according to ANI * conserved DNA values
    if(args.bidirectional):
        filtered_reference_ids = sorted(filtered_reference_ids, key=lambda k: (results[k][0][0] * results[k][0][1] * results[k][1][0] * results[k][1][1]), reverse=True)
        if(args.verbose):
            print('')
        print('#ID\tMash Distance\tQR ANI\tQR Con. DNA\tRQ ANI\tRQ Con. DNA\tTaxonomy ID\tAssembly Status\tOrganism')
        for id in filtered_reference_ids:  # print results to STDOUT
            ref_genome = ref_genomes[id]
            result = results[id]
            print(f"{id}\t{mash_distances[id]:1.5f}\t{(result[0][0] * 100):2.2f}\t{(result[0][1] * 100):2.2f}\t{(result[1][0] * 100):2.2f}\t{(result[1][1] * 100):2.2f}\t{ref_genome['tax']}\t{ref_genome['status']}\t{ref_genome['name']}")
    else:
        filtered_reference_ids = sorted(filtered_reference_ids, key=lambda k: (results[k][0][0] * results[k][0][1]), reverse=True)
        if(args.verbose):
            print('')
        print('#ID\tMash Distance\tANI\tCon. DNA\tTaxonomy ID\tAssembly Status\tOrganism')
        for id in filtered_reference_ids:  # print results to STDOUT
            ref_genome = ref_genomes[id]
            result = results[id][0]
            print(f"{id}\t{mash_distances[id]:1.5f}\t{(result[0] * 100):2.2f}\t{(result[1] * 100):2.2f}\t{ref_genome['tax']}\t{ref_genome['status']}\t{ref_genome['name']}")


if __name__ == '__main__':
    main()
