from pathlib import Path
from subprocess import run


def test_referenceseeker_genome(tmpdir):
    # full test on genome
    proc = run(['bin/referenceseeker', 'test/db', 'test/data/Salmonella_enterica_CFSAN000189.fasta'], capture_output=True, text=True)
    assert proc.returncode == 0
    assert proc.stdout != ''

    lines = proc.stdout.splitlines()
    assert len(lines) == 3
    assert lines[0][0] == '#'  # check if first line is header
    assert lines[1][0] != '#'  # check if second line is a non-header line
    
    cols = lines[1].split('\t')
    assert len(cols) == 7  # check number of colums
    for i in [1,2,3]:  # check Mash dist, ANI and conserved DNA column values
        assert float(cols[i]) <= 100
        assert float(cols[i]) >= 0
    
    assert int(cols[4]) >= 2  # check taxonomy column value
    assert cols[5] in ['complete', 'chromosome', 'scaffold', 'contig']  # check assembly status column values


def test_referenceseeker_genome_zipped(tmpdir):
    # full test on genome
    proc = run(['bin/referenceseeker', 'test/db', 'test/data/Salmonella_enterica_CFSAN000189.fasta.gz'], capture_output=True, text=True)
    assert proc.returncode == 0
    assert proc.stdout != ''

    lines = proc.stdout.splitlines()
    assert len(lines) == 3
    assert lines[0][0] == '#'  # check if first line is header
    assert lines[1][0] != '#'  # check if second line is a non-header line
    
    cols = lines[1].split('\t')
    assert len(cols) == 7  # check number of colums
    for i in [1,2,3]:  # check Mash dist, ANI and conserved DNA column values
        assert float(cols[i]) <= 100
        assert float(cols[i]) >= 0
    
    assert int(cols[4]) >= 2  # check taxonomy column value
    assert cols[5] in ['complete', 'chromosome', 'scaffold', 'contig']  # check assembly status column values


def test_referenceseeker_genome_bidirectional(tmpdir):
    # full test on genome
    proc = run(['bin/referenceseeker', '--bidirectional', 'test/db', 'test/data/Salmonella_enterica_CFSAN000189.fasta'], capture_output=True, text=True)
    assert proc.returncode == 0
    assert proc.stdout != ''

    lines = proc.stdout.splitlines()
    assert len(lines) == 3
    assert lines[0][0] == '#'  # check if first line is header
    assert lines[1][0] != '#'  # check if second line is a non-header line
    
    cols = lines[1].split('\t')
    assert len(cols) == 9  # check number of colums
    for i in [1,2,3,4,5]:  # check Mash dist, ANI and conserved DNA column values
        assert float(cols[i]) <= 100
        assert float(cols[i]) >= 0
    
    assert int(cols[6]) >= 2  # check taxonomy column value
    assert cols[7] in ['complete', 'chromosome', 'scaffold', 'contig']  # check assembly status column values


def test_referenceseeker_genome_bidirectional_zipped(tmpdir):
    # full test on genome
    proc = run(['bin/referenceseeker', '--bidirectional', 'test/db', 'test/data/Salmonella_enterica_CFSAN000189.fasta.gz'], capture_output=True, text=True)
    assert proc.returncode == 0
    assert proc.stdout != ''

    lines = proc.stdout.splitlines()
    assert len(lines) == 3
    assert lines[0][0] == '#'  # check if first line is header
    assert lines[1][0] != '#'  # check if second line is a non-header line
    
    cols = lines[1].split('\t')
    assert len(cols) == 9  # check number of colums
    for i in [1,2,3,4,5]:  # check Mash dist, ANI and conserved DNA column values
        assert float(cols[i]) <= 100
        assert float(cols[i]) >= 0
    
    assert int(cols[6]) >= 2  # check taxonomy column value
    assert cols[7] in ['complete', 'chromosome', 'scaffold', 'contig']  # check assembly status column values