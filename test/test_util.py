
from pathlib import Path

from referenceseeker import constants as rc
from referenceseeker import util as ru


def test_build_dna_fragments(tmpdir):
    genome_path = Path('test/data/Salmonella_enterica_CFSAN000189.fasta').resolve()
    dna_fragments_path = tmpdir / 'fragments.fna'
    dna_fragments = ru.build_dna_fragments(genome_path, dna_fragments_path)

    # first nucleotide fragment must have standard length
    assert dna_fragments[1]['length'] == rc.FRAGMENT_SIZE

    # last nucleotide fragment should have a non-standard length
    last_fragment_id = sorted(list(dna_fragments.keys()))[-1]
    assert dna_fragments[last_fragment_id]['length'] != rc.FRAGMENT_SIZE
