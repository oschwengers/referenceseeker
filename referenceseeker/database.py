
import argparse
from pathlib import Path
import os
import shutil
import subprocess as sp
import sys
import tempfile

from Bio import SeqIO
from xopen import xopen

import referenceseeker
import referenceseeker.constants as rc
import referenceseeker.util as util


def check(db_path):
    """Check if database directory exists, is accessible and contains necessary files."""

    if(db_path is None):
        sys.exit('ERROR: database directory not provided nor detected! Please provide a valid path to the database directory.')

    if(not os.access(str(db_path), os.R_OK & os.X_OK)):
        sys.exit(f'ERROR: database directory ({db_path}) not readable/accessible!')

    file_names = [
        'db.tsv',
        'db.msh'
    ]

    for file_name in file_names:
        path = db_path.joinpath(file_name)
        if(not os.access(str(path), os.R_OK) or not path.is_file()):
            sys.exit(f'ERROR: database file ({file_name}) not readable!')


def init(args):
    try:
        cwd_path = Path(args.output).resolve()
        db_path = cwd_path.joinpath(args.db)
        if(db_path.exists()):
            print(f'detected existing database directory: {db_path}')
            if(not db_path.is_dir()):
                sys.exit(f'Database path ({db_path}) exists but is not a directory!')
        else:
            db_path.mkdir(mode=0o770)
            print('created database directory')

        db_tsv_path = db_path.joinpath('db.tsv')
        db_tsv_path.touch(mode=0o660)
        print('created database genome information file (db.tsv)')

        db_sketch_path = db_path.joinpath('db.msh')
        db_sketch_path.touch(mode=0o660)
        print('created database genome kmer sketch file (db.msh)')
    except Exception as e:
        print(f'Error: could not init database ({args.output}) in output directory ({args.db})!', file=sys.stderr)
        raise e
    print(f'\nSuccessfully initialized empty database at {db_path}')
    print("Use 'referenceseeker_db import' to import genomes into database")


def import_genome(args):
    try:
        db_path = Path(args.db).resolve()
        check(db_path)
        tmp_path = Path(tempfile.mkdtemp())
        genome_id = args.id
        genome_path = Path(args.genome).resolve()
        genome_suffix = genome_path.suffix.lower()
        if(genome_suffix == '.gz'):
            genome_suffix = Path(genome_path.stem).suffix
        if(genome_suffix in ['.fasta', '.fas', '.fsa', '.fna', '.fa']):
            # import fasta
            input_path = genome_path
            genome_path = tmp_path.joinpath('genome.fasta')
            with xopen(str(input_path), threads=0) as fh_in, genome_path.open('w') as fh_out:
                sequences_in = SeqIO.parse(fh_in, 'fasta')
                sequences = test_sequences(sequences_in)
                for id, seq in sequences:
                    fh_out.write(f'>{id}\n{seq}\n')
        elif(genome_suffix in ['.genbank', '.gbff', '.gbk', '.gb']):
            # import genbank
            input_path = genome_path
            genome_path = tmp_path.joinpath('genome.fasta')
            with xopen(str(input_path), threads=0) as fh_in, genome_path.open('w') as fh_out:
                sequences_in = SeqIO.parse(fh_in, 'genbank')
                sequences = test_sequences(sequences_in)
                for id, seq in sequences:
                    fh_out.write(f'>{id}\n{seq}\n')
        elif(genome_suffix in ['.embl', '.ebl', '.el']):
            # import embl
            input_path = genome_path
            genome_path = tmp_path.joinpath('genome.fasta')
            with xopen(str(input_path), threads=0) as fh_in, genome_path.open('w') as fh_out:
                sequences_in = SeqIO.parse(fh_in, 'embl')
                sequences = test_sequences(sequences_in)
                for id, seq in sequences:
                    fh_out.write(f'>{id}\n{seq}\n')
        else:
            raise Exception(f'Unknown genome file extension ({genome_suffix})')

        # extract genome id if it is empty
        if(genome_id is None):
            for record in SeqIO.parse(str(genome_path), 'fasta'):
                genome_id = record.id
                break

        # copy genome fasta file to database directory
        genome_database_path = db_path.joinpath(f'{genome_id}.fna.gz')
        with genome_path.open() as fh_in, xopen(str(genome_database_path), mode='wb') as fh_out:
            for line in fh_in:
                fh_out.write(bytes(line, 'utf-8'))

        # copy/rename genome file to have Mash use the right ID in its internal db
        old_genome_path = genome_path
        genome_path = tmp_path.joinpath(genome_id)
        shutil.copyfile(str(old_genome_path), str(genome_path))

        # sketch genome
        cmd = [
            'mash',
            'sketch',
            '-o', 'genome',
            '-k', '32',
            '-s', '10000',
            genome_id
        ]
        proc = sp.run(
            cmd,
            cwd=str(tmp_path),
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            universal_newlines=True
        )
        if(proc.returncode != 0):
            sys.exit(f'ERROR: failed to create genome kmer sketches via Mash!\nexit={proc.returncode}\ncmd={cmd}')

        existing_db_sketch_path = db_path.joinpath('db.msh')
        if(existing_db_sketch_path.stat().st_size == 0):
            # empty db -> replace existing with new sketch file
            shutil.move(str(tmp_path.joinpath('genome.msh')), str(existing_db_sketch_path))
        else:
            # add genome sketch to db sketch file
            cmd = [
                'mash',
                'paste',
                'db',
                str(existing_db_sketch_path),
                str(tmp_path.joinpath('genome.msh'))
            ]
            proc = sp.run(
                cmd,
                cwd=str(tmp_path),
                stdout=sp.PIPE,
                stderr=sp.PIPE,
                universal_newlines=True
            )
            if(proc.returncode != 0):
                sys.exit(f'ERROR: failed to import sketch to database via Mash!\nexit={proc.returncode}\ncmd={cmd}')
            # replace old database sketch file by new db.msh
            shutil.move(str(tmp_path.joinpath('db.msh')), str(existing_db_sketch_path))

        # remove tmp dir
        shutil.rmtree(str(tmp_path))

        # add genome metainformation to db.tsv
        with db_path.joinpath('db.tsv').open(mode='a') as fh:
            fh.write(f'{genome_id}\t{args.taxonomy}\t{args.status}\t{args.organism}\n')
    except Exception as e:
        print(f'ERROR: could not import genome ({args.organism}/{args.genome}) into database ({args.db})!', file=sys.stderr)
        raise e
    print(f'\nSuccessfully imported genome ({args.organism}/{args.genome}) into database ({db_path})')


def test_sequences(sequences_in):
    sequence_ids = set()
    sequences = []
    for record in sequences_in:
        if(len(record.seq) == 0):
            raise Exception(f'Record {record.id} with zero length sequence')
        if(record.id in sequence_ids):
            raise Exception(f'Duplicated record id: {record.id}')
        else:
            sequences.append((record.id, str(record.seq)))
            sequence_ids.add(record.id)
    return sequences


def main():
    #  setup path and test if necessary 3rd party executables are available
    util.test_binaries()

    # parse options and arguments
    parser = argparse.ArgumentParser(
        prog='referenceseeker_db',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Rapid determination of appropriate reference genomes.',
        epilog=f'Citation:\n{rc.CITATION}\n\nGitHub:\nhttps://github.com/oschwengers/referenceseeker',
        add_help=False
    )
    #  add common options
    group_runtime = parser.add_argument_group('Runtime & auxiliary options')
    group_runtime.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    group_runtime.add_argument('--version', '-V', action='version', version=f'%(prog)s {referenceseeker.__version__}')

    subparsers = parser.add_subparsers(dest='subcommand', help='sub-command help')
    #  add init sub-command options
    parser_init = subparsers.add_parser('init', help='Initialize a new database')
    parser_init.add_argument('--output', '-o', action='store', default=Path.cwd(), help='output directory (default = current working directory)')
    parser_init.add_argument('--db', '-d', action='store', required=True, help='Name of the new ReferenceSeeker database')

    #  add import sub-command options
    parser_import = subparsers.add_parser('import', help='Add a new genome to database')
    parser_import.add_argument('--db', '-d', action='store', required=True, help='ReferenceSeeker database path')
    parser_import.add_argument('--genome', '-g', action='store', required=True, help='Genome path [Fasta, GenBank, EMBL]')
    parser_import.add_argument('--id', '-i', action='store', default=None, help='Unique genome identifier (default sequence id of first record)')
    parser_import.add_argument('--taxonomy', '-t', action='store', type=int, default=12908, help='Taxonomy ID (default = 12908 [unclassified sequences])')
    parser_import.add_argument('--status', '-s', action='store', choices=['complete', 'chromosome', 'scaffold', 'contig'], default='contig', help='Assembly level (default = contig)')
    parser_import.add_argument('--organism', '-o', action='store', default='NA', help='Organism name (default = "NA")')

    args = parser.parse_args()

    if(args.subcommand == 'init'):
        init(args)
    elif(args.subcommand == 'import'):
        import_genome(args)
    else:
        parser.print_help()
        sys.exit('Error: no subcommand provided!')


if __name__ == '__main__':
    main()
