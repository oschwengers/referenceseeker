
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


def compute_n50(fasta_path):
    """Calculate N50 metric.

    :param fasta_path: Path to sequence Fasta file.

    :rtype (N50,L50): A tuple containing the N50 and L50 metrics for the assembly in 'fasta_path'.
    """

    genome_length = 0
    contig_lengths = []
    for record in SeqIO.parse(fasta_path, 'fasta'):
        length = len(record.seq)
        genome_length += length
        contig_lengths.append(length)
    contig_lengths = sorted(contig_lengths, reverse=True)
    tmp_sum = 0
    i = -1
    while tmp_sum <= 0.5 * genome_length:
        i += 1
        tmp_sum += contig_lengths[i]
    return contig_lengths[i], i + 1
