
import os
import subprocess as sp
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO

import referenceseeker.constants as rc


def read_reference_genomes(config, accession_ids, mash_distances):
    ref_genomes = []
    with open(str(config['db_path'].joinpath('db.tsv')), 'r') as fh:
        for line in fh:
            if(line[0] != '#'):
                cols = line.strip().split('\t')
                accession_id = cols[0]
                if(accession_id in accession_ids):
                    ref_genomes.append(
                        {
                            'id': accession_id,
                            'tax': cols[1],
                            'status': cols[2],
                            'name': cols[3],
                            'mash_dist': mash_distances[accession_id]
                        }
                    )
    return ref_genomes


def build_dna_fragments(config, dna_fragments_path):
    """Build DNA fragments.

    :param config: a global config object encapsulating global runtime vars
    :param dna_fragments_path: Path to DNA fragments output Fasta file.

    :rtype {idx, dna_fragment}: A dict containing index and DNA fragment objects.
    """

    dna_fragments = {}
    dna_fragment_idx = 1
    genome_path = config['genome_path']
    with dna_fragments_path.open(mode='w') as fh:
        for record in SeqIO.parse(str(genome_path), 'fasta'):
            sequence = record.seq
            while len(sequence) > (rc.FRAGMENT_SIZE + rc.MIN_FRAGMENT_SIZE):  # forestall fragments shorter than MIN_FRAGMENT_SIZE
                dnaFragment = sequence[:rc.FRAGMENT_SIZE]
                fh.write(">%s\n%s\n" % (str(dna_fragment_idx), str(dnaFragment)))
                dna_fragments[dna_fragment_idx] = {
                    'id': dna_fragment_idx,
                    'length': len(dnaFragment)
                }
                sequence = sequence[rc.FRAGMENT_SIZE:]
                dna_fragment_idx += 1
            dnaFragment = sequence
            fh.write(">%s\n%s\n" % (str(dna_fragment_idx), str(dnaFragment)))
            dna_fragments[dna_fragment_idx] = {
                'id': dna_fragment_idx,
                'length': len(dnaFragment)
            }
            sequence = sequence[rc.FRAGMENT_SIZE:]
            dna_fragment_idx += 1
    return dna_fragments


def setup_configuration(args):
    """Test environment and build a runtime configuration."""

    config = {
        'env': os.environ.copy(),
        'tmp': Path(tempfile.mkdtemp()),
        'bundled-binaries': False,
        'threads': args.threads,
        'unfiltered': args.unfiltered,
        'crg': args.crg,
        'ani': args.ani,
        'conserved_dna': args.conserved_dna
    }
    base_dir = Path(__file__).parent.parent
    share_dir = base_dir.joinpath('share')
    if(os.access(str(share_dir), os.R_OK & os.X_OK)):
        config['env']["PATH"] = str(share_dir) + ':' + config['env']["PATH"]
        config['bundled-binaries'] = True

    db_path = base_dir.joinpath('db')
    if(os.access(str(db_path), os.R_OK & os.X_OK)):
        config['db'] = db_path
    return config


def test_binaries():
    """Test the proper installation of necessary 3rd party executables."""

    # test prodigal
    try:
        sp.check_call(
            ['mash', 'dist', '-h'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        sys.exit('ERROR: \'prodigal\' was not found!')
    except:
        sys.exit('ERROR: \'prodigal\' was not exeutable!')

    # test nucmer
    try:
        sp.check_call(
            ['nucmer', '--help'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        sys.exit('ERROR: \'nucmer\' was not found!')
    except:
        sys.exit('ERROR: \'nucmer\' was not exeutable!')

    # test delta-filter
    try:
        sp.check_call(
            ['delta-filter', '-h'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        sys.exit('ERROR: \'delta-filter\' was not found!')
    except:
        sys.exit('ERROR: \'delta-filter\' was not exeutable!')

    return
