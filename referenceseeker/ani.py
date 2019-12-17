
import shutil
import sys
import subprocess as sp
import tempfile

from pathlib import Path


def compute_ani(config, dna_fragments_path, dna_fragments, ref_genome):
    """Perform per-genome calculation of ANI/conserved DNA values.

    :param config: a global config object encapsulating global runtime vars
    :param dna_fragments: A dict comprising information on fragments.
    :param ref_genome: A dict representing a reference genome.

    :rtype: A dict representing a reference genome and additionally comprising ANI / conserved DNA values.
    """

    reference_path = config['db_path'].joinpath("%s.fna" % ref_genome['id'])
    tmp_dir = Path(tempfile.mkdtemp())

    # perform global alignments via nucmer
    cmd = [
        'nucmer',
        '--threads=1',
        str(reference_path),
        str(dna_fragments_path)
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

    # calc % conserved DNA
    alignment_sum = 0
    for fm in dna_fragment_matches:
        rel_alignment_length = float(fm['alignment_length'] - fm['no_non_identities']) / float(fm['length'])
        if(rel_alignment_length > 0.9):
            alignment_sum += fm['alignment_length']
    genome_length = 0
    for fm in dna_fragments.values():
        genome_length += fm['length']
    conserved_dna = (float(alignment_sum) / float(genome_length)) if genome_length > 0 else 0

    # calc average nucleotide identity
    ani_matches = 0
    ni_sum = 0.0
    for fm in dna_fragment_matches:
        if(((float(fm['alignment_length'] - fm['no_non_identities']) / float(fm['length'])) > 0.3)
                and ((float(fm['alignment_length']) / float(fm['length'])) >= 0.7)):
            ni_sum += float(fm['alignment_length'] - fm['no_non_identities']) / float(fm['alignment_length'])
            ani_matches += 1
    ani = (ni_sum / float(ani_matches)) if ani_matches > 0 else 0

    shutil.rmtree(str(tmp_dir))

    # if(args.verbose):
    #     print(
    #         '\t%s\t%2.2f\t%2.2f\t%1.5f' %
    #         (
    #             reference_path.split('/')[-1][:15],
    #             ani * 100,
    #             conserved_dna * 100,
    #             ref_genome['mash_dist'],
    #         )
    #     )

    ref_genome['ani'] = ani
    ref_genome['conserved_dna'] = conserved_dna

    return ref_genome
