
import pytest

from referenceseeker import ani as ra


@pytest.fixture(scope='module')
def test_fragment_matches():
    dna_fragment_matches = [
        {  # ANI: valid, conDNA: valid
            'id': 1,
            'alignment_length': 950,
            'no_non_identities': 10,
            'length': 1000
        },
        {  # ANI: valid, conDNA: valid
            'id': 2,
            'alignment_length': 951,
            'no_non_identities': 50,
            'length': 1000
        },
        {  # ANI: not valid as seq cov is < 0.7, conDNA: not valid as cov < 0.9
            'id': 3,
            'alignment_length': 500,
            'no_non_identities': 10,
            'length': 1000
        },
        {  # ANI: not valid as seq identity is < 0.3, conDNA: not valid as cov < 0.9
            'id': 4,
            'alignment_length': 700,
            'no_non_identities': 600,
            'length': 1000
        }
    ]
    return dna_fragment_matches


def test_calculate_conserved_dna_all(test_fragment_matches):
    #  check conDNA over all test fragment matches
    conserved_dna = ra.calculate_conserved_dna({k['id']: k for k in test_fragment_matches}, test_fragment_matches)
    valid_dna_fragment_matches = test_fragment_matches[0:2]
    expected_conserved_dna = sum(map(lambda k: k['alignment_length'], valid_dna_fragment_matches)) / sum(map(lambda k: k['length'], test_fragment_matches))
    assert conserved_dna == expected_conserved_dna


def test_calculate_conserved_dna_valid(test_fragment_matches):
    #  check conDNA for *valid* test fragment matches
    test_fragment_matches = test_fragment_matches[0:2]
    conserved_dna = ra.calculate_conserved_dna({k['id']: k for k in test_fragment_matches}, test_fragment_matches)
    expected_conserved_dna = sum(map(lambda k: k['alignment_length'], test_fragment_matches)) / sum(map(lambda k: k['length'], test_fragment_matches))
    assert conserved_dna == expected_conserved_dna


def test_calculate_conserved_dna_unvalid(test_fragment_matches):
    #  check conDNA for *unvalid* test fragment matches
    test_fragment_matches = test_fragment_matches[2:]
    conserved_dna = ra.calculate_conserved_dna({k['id']: k for k in test_fragment_matches}, test_fragment_matches)
    expected_conserved_dna = 0.0
    assert conserved_dna == expected_conserved_dna


def test_calculate_conserved_dna_empty():
    #  check conDNA for empty test fragment matches
    conserved_dna = ra.calculate_conserved_dna({}, [])
    expected_conserved_dna = 0.0
    assert conserved_dna == expected_conserved_dna


def test_calculate_ani_all(test_fragment_matches):
    #  check ANI over all test fragment matches
    ani = ra.calculate_ani(test_fragment_matches)
    valid_dna_fragment_matches = test_fragment_matches[0:2]
    expected_ani = sum(map(lambda k: (k['alignment_length'] - k['no_non_identities']) / k['alignment_length'], valid_dna_fragment_matches)) / len(valid_dna_fragment_matches)
    assert ani == expected_ani


def test_calculate_ani_valid(test_fragment_matches):
    #  check ANI for *valid* test fragment matches
    test_fragment_matches = test_fragment_matches[0:2]
    ani = ra.calculate_ani(test_fragment_matches)
    expected_ani = sum(map(lambda k: (k['alignment_length'] - k['no_non_identities']) / k['alignment_length'], test_fragment_matches)) / len(test_fragment_matches)
    assert ani == expected_ani


def test_calculate_ani_unvalid(test_fragment_matches):
    #  check ANI for *unvalid* test fragment matches
    test_fragment_matches = test_fragment_matches[2:]
    ani = ra.calculate_ani(test_fragment_matches)
    expected_ani = 0.0
    assert ani == expected_ani


def test_calculate_ani_empty():
    #  check ANI for empty test fragment matches
    ani = ra.calculate_ani([])
    expected_ani = 0.0
    assert ani == expected_ani
