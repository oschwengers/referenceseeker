
import os
import subprocess as sp
import sys
import tempfile
from pathlib import Path
from xopen import xopen

from Bio import SeqIO

import referenceseeker.constants as rc


def read_reference_genomes(config):
    ref_genomes = {}
    with open(str(config['db_path'].joinpath('db.tsv')), 'r') as fh:
        for line in fh:
            if(line[0] != '#'):
                cols = line.strip().split('\t')
                accession_id = cols[0]
                ref_genomes[accession_id] = {
                    'id': accession_id,
                    'tax': cols[1],
                    'status': cols[2],
                    'name': cols[3]
                }
    return ref_genomes


def build_dna_fragments(genome_path, dna_fragments_path):
    """Build DNA fragments.

    :param genome_path: Path to source DNA Fasta file.
    :param dna_fragments_path: Path to DNA fragments output Fasta file.

    :rtype {idx, dna_fragment}: A dict containing index and DNA fragment objects.
    """

    dna_fragments = {}
    dna_fragment_idx = 0
    with dna_fragments_path.open(mode='w') as fh_out, xopen(str(genome_path), threads=0) as fh_in:
        for record in SeqIO.parse(fh_in, 'fasta'):
            sequence = str(record.seq)
            while len(sequence) > (rc.FRAGMENT_SIZE + rc.MIN_FRAGMENT_SIZE):  # forestall fragments shorter than MIN_FRAGMENT_SIZE
                dna_fragment = sequence[:rc.FRAGMENT_SIZE]
                dna_fragment_idx += 1
                fh_out.write('>')
                fh_out.write(str(dna_fragment_idx))
                fh_out.write('\n')
                fh_out.write(str(dna_fragment))
                fh_out.write('\n')
                dna_fragments[dna_fragment_idx] = {
                    'id': dna_fragment_idx,
                    'length': len(dna_fragment)
                }
                sequence = sequence[rc.FRAGMENT_SIZE:]
            dna_fragment = sequence
            dna_fragment_idx += 1
            fh_out.write('>')
            fh_out.write(str(dna_fragment_idx))
            fh_out.write('\n')
            fh_out.write(str(dna_fragment))
            fh_out.write('\n')
            dna_fragments[dna_fragment_idx] = {
                'id': dna_fragment_idx,
                'length': len(dna_fragment)
            }
            # sequence = sequence[rc.FRAGMENT_SIZE:]
    return dna_fragments


def setup_configuration(args):
    """Test environment and build a runtime configuration."""

    config = {
        'tmp': Path(tempfile.mkdtemp()),
        'threads': args.threads,
        'unfiltered': args.unfiltered,
        'bidirectional': args.bidirectional,
        'crg': args.crg,
        'ani': args.ani,
        'conserved_dna': args.conserved_dna
    }

    base_dir = Path(__file__).parent.parent
    db_path = base_dir.joinpath('db')
    if(os.access(str(db_path), os.R_OK & os.X_OK)):
        config['db'] = db_path
    return config


def test_binaries(config):
    """Test the proper installation of necessary 3rd party executables."""

    # test Mash
    try:
        sp.check_call(
            ['mash', 'dist', '-h'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        sys.exit('ERROR: \'Mash\' was not found!')
    except sp.CalledProcessError as err:
        sys.exit(f'ERROR: \'Mash\' test execution failed! error-code={err.returncode}')
    except:
        sys.exit('ERROR: \'Mash\' was not executable!')

    # test nucmer
    try:
        sp.check_call(
            ['nucmer', '--help'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        sys.exit('ERROR: \'nucmer\' was not found!')
    except sp.CalledProcessError as err:
        sys.exit(f'ERROR: \'nucmer\' test execution failed! error-code={err.returncode}')
    except:
        sys.exit('ERROR: \'nucmer\' was not executable!')

    # test delta-filter
    try:
        sp.check_call(
            ['delta-filter', '-h'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        sys.exit('ERROR: \'delta-filter\' was not found!')
    except sp.CalledProcessError as err:
        sys.exit(f'ERROR: \'delta-filter\' test execution failed! error-code={err.returncode}')
    except:
        sys.exit('ERROR: \'delta-filter\' was not executable!')

    return
