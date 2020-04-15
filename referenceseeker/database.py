
import argparse
from pathlib import Path
import shutil
import subprocess as sp
import sys
import tempfile

from Bio import SeqIO

import referenceseeker
import referenceseeker.constants as rc
import referenceseeker.util as util


def init(args):
    try:
        cwd_path = Path(args.output).resolve()
        db_path = cwd_path.joinpath(args.db)
        if(db_path.exists()):
            print("detected existing database directory: %s" % db_path)
            if(not db_path.is_dir()):
                sys.exit("Database path (%s) exists but is not a directory!" % db_path)
        else:
            db_path.mkdir(mode=0o770)
            print("created database directory")

        db_tsv_path = db_path.joinpath('db.tsv')
        db_tsv_path.touch(mode=0o660)
        print("created database genome information file (db.tsv)")

        db_sketch_path = db_path.joinpath('db.msh')
        db_sketch_path.touch(mode=0o660)
        print("created database genome kmer sketch file (db.msh)")
    except:
        print("Error: could not init database (%s) in output directory (%s)!" % (args.output, args.db), file=sys.stderr)
        raise
        sys.exit(-1)
    print("\nSuccessfully initialized empty database at %s" % db_path)
    print("Use 'referenceseeker_db import' to import genomes into database")


def import_genome(config, args):
    try:
        db_path = Path(args.db).resolve()
        genome_path = Path(args.genome).resolve()
        tmp_path = Path(tempfile.mkdtemp())
        genome_suffix = genome_path.suffix.lower()
        genome_id = args.id
        if(genome_suffix in ['.fasta', '.fas', '.fsa', '.fna', '.fa']):
            # import fasta
            with genome_path.open() as fh_in:
                sequences = SeqIO.parse(fh_in, "fasta")
                test_sequences(sequences)
        elif(genome_suffix in ['.genbank', '.gbff', '.gbk', '.gb']):
            # import genbank
            input_path = genome_path
            genome_path = tmp_path.joinpath('genome.fasta')
            with input_path.open() as fh_in, genome_path.open('w') as fh_out:
                sequences = SeqIO.parse(fh_in, "genbank")
                test_sequences(sequences)
                SeqIO.write(sequences, fh_out, "fasta")
        elif(genome_suffix in ['.embl', '.ebl', '.el']):
            # import embl
            input_path = genome_path
            genome_path = tmp_path.joinpath('genome.fasta')
            with input_path.open() as fh_in, genome_path.open('w') as fh_out:
                sequences = SeqIO.parse(fh_in, "embl")
                test_sequences(sequences)
                SeqIO.write(sequences, fh_out, "fasta")
        else:
            raise Exception("Unknown genome file extension (%s)" % genome_suffix)

        # extract genome id if it is empty
        if(genome_id is None):
            for record in SeqIO.parse(str(genome_path), 'fasta'):
                genome_id = record.id
                break

        # copy genome fasta file to database directory
        shutil.copyfile(str(genome_path), str(db_path.joinpath("%s.fna" % genome_id)))

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
            env=config['env'],
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            universal_newlines=True
        )
        if(proc.returncode != 0):
            sys.exit("ERROR: failed to create genome kmer sketches via Mash!\nexit=%d\ncmd=%s" % (proc.returncode, cmd))

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
                env=config['env'],
                stdout=sp.PIPE,
                stderr=sp.PIPE,
                universal_newlines=True
            )
            if(proc.returncode != 0):
                sys.exit("ERROR: failed to import sketch to database via Mash!\nexit=%d\ncmd=%s" % (proc.returncode, cmd))
            # replace old database sketch file by new db.msh
            shutil.move(str(tmp_path.joinpath('db.msh')), str(existing_db_sketch_path))

        # remove tmp dir
        shutil.rmtree(str(tmp_path))

        # add genome metainformation to db.tsv
        with db_path.joinpath('db.tsv').open(mode='a') as fh:
            fh.write(
                "%s\t%i\t%s\t%s\n" %
                (genome_id, args.taxonomy, args.status, args.organism)
            )
    except Exception as e:
        print(e)
        print("ERROR: could not import genome (%s/%s) into database (%s)!" % (args.organism, args.genome, args.db), file=sys.stderr)
        sys.exit(-1)
    print("\nSuccessfully imported genome (%s/%s) into database (%s)" % (args.organism, args.genome, db_path))


def test_sequences(sequences):
    sequence_ids = set()
    for record in sequences:
        if(len(record.seq) == 0):
            raise Exception("Record %s with zero length sequence" % record.id)
        if(record.id in sequence_ids):
            raise Exception("Duplicated record id: %s" % record.id)
        else:
            sequence_ids.add(record.id)


def main():
    #  setup path and test if necessary 3rd party executables are available
    config = {}
    util.set_path(config)
    util.test_binaries(config)

    # parse options and arguments
    parser = argparse.ArgumentParser(
        prog='referenceseeker_db',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Rapid determination of appropriate reference genomes.',
        epilog="Citation:\n%s\n\nGitHub:\nhttps://github.com/oschwengers/referenceseeker" % rc.CITATION,
        add_help=False
    )
    #  add common options
    group_runtime = parser.add_argument_group('Runtime & auxiliary options')
    group_runtime.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    group_runtime.add_argument('--version', '-V', action='version', version='%(prog)s ' + referenceseeker.__version__)

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
    parser_import.add_argument('--organism', '-o', action='store', default='', help='Organism name (default = "")')

    args = parser.parse_args()

    if(args.subcommand == 'init'):
        init(args)
    elif(args.subcommand == 'import'):
        import_genome(config, args)
    else:
        parser.print_help()
        sys.exit("Error: no subcommand provided!")


if __name__ == '__main__':
    main()
