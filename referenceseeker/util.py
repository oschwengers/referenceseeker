
from Bio import SeqIO

import referenceseeker.constants as rc


def read_reference_genomes(db_path, accession_ids, mash_distances):
    ref_genomes = []
    with open(db_path + '/db.tsv', 'r') as fh:
        for line in fh:
            if line[0] != '#':
                cols = line.strip().split('\t')
                accession_id = cols[0]
                if accession_id in accession_ids:
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


def build_dna_fragments(genome_path, dna_fragments_path):
    """Build DNA fragments.

    :param genome_path: Path to input sequence Fasta file.
    :param dna_fragments_path: Path to DNA fragments output Fasta file.

    :rtype {idx, dna_fragment}: A dict containing index and DNA fragment objects.
    """

    dna_fragments = {}
    dna_fragment_idx = 1
    with open(dna_fragments_path, 'w') as fh:
        for record in SeqIO.parse(genome_path, 'fasta'):
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
