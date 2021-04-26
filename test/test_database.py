import pytest

from pathlib import Path
from subprocess import run


def test_database(tmpdir):
    # init test db
    proc = run(['bin/referenceseeker_db', 'init', '--output', tmpdir, '--db', 'test'])
    assert proc.returncode == 0
    
    tmpdir_path = Path(tmpdir)
    db_path = tmpdir_path.joinpath('test')
    assert Path.exists(db_path)
    for file in ['db.tsv', 'db.msh']:  # check file existence
        assert Path.exists(db_path.joinpath(file))
    
    # import genome
    proc = run(['bin/referenceseeker_db', 'import', '--db', f'{tmpdir}/test', '--genome', 'test/data/Salmonella_enterica_CFSAN000189.fasta'])
    assert proc.returncode == 0

    for file in ['db.tsv', 'db.msh']:  # check db files are not empty
        file_path = db_path.joinpath(file)
        assert Path.exists(file_path)
        assert file_path.stat().st_size > 0
    
    # import genome
    proc = run(['bin/referenceseeker_db', 'import', '--db', f'{tmpdir}/test', '--genome', 'test/data/Salmonella_enterica_CFSAN000189.fasta.gz'])
    assert proc.returncode == 0

    for file in ['db.tsv', 'db.msh']:  # check db files are not empty
        file_path = db_path.joinpath(file)
        assert Path.exists(file_path)
        assert file_path.stat().st_size > 0
    
    # test new database
    proc = run(['bin/referenceseeker', f'{tmpdir}/test', 'test/data/Salmonella_enterica_CFSAN000189.fasta.gz'], capture_output=True, text=True)
    assert proc.returncode == 0
    assert proc.stdout != ''

    lines = proc.stdout.splitlines()
    assert len(lines) == 2
    assert lines[0][0] == '#'  # check if first line is header
    assert lines[1][0] != '#'  # check if second line is a non-header line
    
    cols = lines[1].split('\t')
    assert len(cols) == 7  # check number of colums
    for i in [1,2,3]:  # check Mash dist, ANI and conserved DNA column values
        assert float(cols[i]) <= 100
        assert float(cols[i]) >= 0
    
    assert int(cols[4]) >= 2  # check taxonomy column value
    assert cols[5] in ['complete', 'chromosome', 'scaffold', 'contig']  # check assembly status column values