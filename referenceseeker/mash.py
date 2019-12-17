
import subprocess as sp
import sys

import referenceseeker.constants as rc


def run_mash(config, mash_output_path):
    with mash_output_path.open(mode='w') as fh:
        cmd = [
            'mash',
            'dist',
            '-d', rc.UNFILTERED_MASH_DIST if config['unfiltered'] else rc.MAX_MASH_DIST,
            '-p', str(config['threads']),
            str(config['db_path'].joinpath('db.msh')),
            str(config['genome_path'])
        ]
        proc = sp.run(
            cmd,
            cwd=str(config['tmp']),
            env=config['env'],
            stdout=fh,
            stderr=sp.PIPE,
            universal_newlines=True
        )
        if(proc.returncode != 0):
            sys.exit("ERROR: failed to execute Mash!\nexit=%d\ncmd=%s" % (proc.returncode, cmd))


def parse_mash_results(config, mash_output_path):
    accession_ids = []
    mash_distances = {}
    with mash_output_path.open() as fh:
        for line in fh:
            cols = line.rstrip().split()
            accession_ids.append(cols[0])
            mash_distances[cols[0]] = float(cols[2])
    return accession_ids, mash_distances
