
import shutil
import sys
import subprocess as sp
import tempfile

from pathlib import Path

import referenceseeker.util as util


def align_query_genome(config, dna_fragments_path, dna_fragments, ref_genome_id):
    """Perform per-genome calculation of ANI/conserved DNA values.

    :param config: a global config object encapsulating global runtime vars
    :param dna_fragments: A dict comprising information on fragments.
    :param ref_genome_id: reference genome id.

    :rtype: A dict representing a reference genome and additionally comprising ANI / conserved DNA values.
    """

    reference_genome_path = config['db_path'].joinpath("%s.fna" % ref_genome_id)
    tmp_dir = Path(tempfile.mkdtemp())

    dna_fragment_matches = execute_nucmer(config, tmp_dir, dna_fragments, dna_fragments_path, reference_genome_path)

    shutil.rmtree(str(tmp_dir))

    ani = calculate_ani(dna_fragment_matches)
    conserved_dna = calculate_conserved_dna(dna_fragments, dna_fragment_matches)

    return (ref_genome_id, ani, conserved_dna)


def align_reference_genome(config, query_genome_path, ref_genome_id):
    """Perform per-genome calculation of ANI/conserved DNA values.

    :param config: a global config object encapsulating global runtime vars
    :param query_genome_path: Path to query genome Fasta file.
    :param ref_genome_id: reference genome id.

    :rtype: A dict representing a reference genome and additionally comprising ANI / conserved DNA values.
    """

    reference_genome_path = config['db_path'].joinpath("%s.fna" % ref_genome_id)
    tmp_dir = Path(tempfile.mkdtemp())

    dna_fragments_path = tmp_dir.joinpath('dna-fragments.fasta')
    dna_fragments = util.build_dna_fragments(reference_genome_path, dna_fragments_path)

    # perform global alignments via nucmer
    dna_fragment_matches = execute_nucmer(config, tmp_dir, dna_fragments, dna_fragments_path, query_genome_path)

    shutil.rmtree(str(tmp_dir))

    ani = calculate_ani(dna_fragment_matches)
    conserved_dna = calculate_conserved_dna(dna_fragments, dna_fragment_matches)

    return (ref_genome_id, ani, conserved_dna)


def execute_nucmer(config, tmp_dir, dna_fragments, query_path, reference_genome_path):
    cmd = [
        'nucmer',
        '--threads=1',
        str(reference_genome_path),
        str(query_path)
    ]
    proc = sp.run(
        cmd,
        cwd=str(tmp_dir),
        env=config['env'],
        stdout=sp.DEVNULL,
        stderr=sp.STDOUT,
        universal_newlines=True
    )
    if(proc.returncode != 0):
        sys.exit("ERROR: failed to execute nucmer!\nexit=%d\ncmd=%s" % (proc.returncode, cmd))

    filtered_delta_path = tmp_dir.joinpath('out-filtered.delta')
    with filtered_delta_path.open(mode='w') as fh:
        cmd = [
            'delta-filter',
            '-q',
            'out.delta'
        ]
        proc = sp.run(
            cmd,
            cwd=str(tmp_dir),
            env=config['env'],
            stdout=fh,
            stderr=sp.STDOUT
        )
        if(proc.returncode != 0):
            sys.exit("ERROR: failed to execute delta-filter!\nexit=%d\ncmd=%s" % (proc.returncode, cmd))

    # parse nucmer output
    dna_fragment = None
    dna_fragment_matches = []
    with filtered_delta_path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if(line[0] == '>'):
                dna_fragment = dna_fragments.get(int(line.split(' ')[1]), None)
            elif(dna_fragment is not None):
                cols = line.split(' ')
                if(len(cols) == 7):
                    dna_fragment['alignment_length'] = abs(int(cols[3]) - int(cols[2])) + 1  # abs( qStop - qStart ) + 1
                    dna_fragment['no_non_identities'] = int(cols[4])  # number of non-identities
                    dna_fragment_matches.append(dna_fragment)

    return dna_fragment_matches


def calculate_conserved_dna(dna_fragments, dna_fragment_matches):
    """Calculate conserved DNA value for a set of DNA fragment matches.

    :param dna_fragments: A dict comprising information on DNA fragments.
    :param dna_fragment_matches: A dict comprising matches of DNA fragments.

    :rtype: conserved DNA value.
    """

    alignment_sum = 0
    for fm in dna_fragment_matches:
        rel_alignment_length = float(fm['alignment_length'] - fm['no_non_identities']) / float(fm['length'])
        if(rel_alignment_length > 0.9):
            alignment_sum += fm['alignment_length']
    genome_length = 0
    for fm in dna_fragments.values():
        genome_length += fm['length']

    return (float(alignment_sum) / float(genome_length)) if genome_length > 0 else 0


def calculate_ani(dna_fragment_matches):
    """Calculate ANI value for a set of DNA fragment matches.

    :param dna_fragment_matches: A dict comprising matches of DNA fragments.

    :rtype: ANI value.
    """

    ani_matches = 0
    ni_sum = 0.0
    for fm in dna_fragment_matches:
        if(((float(fm['alignment_length'] - fm['no_non_identities']) / float(fm['length'])) > 0.3)
                and ((float(fm['alignment_length']) / float(fm['length'])) >= 0.7)):
            ni_sum += float(fm['alignment_length'] - fm['no_non_identities']) / float(fm['alignment_length'])
            ani_matches += 1

    return (ni_sum / float(ani_matches)) if ani_matches > 0 else 0
