import os
import pytest

from subprocess import run


@pytest.mark.parametrize(
    'parameters',
    [
        ([]),  # no parameter
        (['']),  # empty argument
        (['foo.fasta'])  # argument not existing
    ]
)
def test_genome_failing(parameters, tmpdir):
    # test genome arguments
    cmd_line = ['bin/referenceseeker', '--output', tmpdir, 'test/db'] + parameters
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        ([]),  # not provided
        (['', ]),  # empty
        (['test/foo']),  # not existing
    ]
)
def test_database_failing_parameter(parameters, tmpdir):
    # test database arguments

    cmd_line = ['bin/referenceseeker'] + parameters + ['test/data/Salmonella_enterica_CFSAN000189.fasta']
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'arguments',
    [
        ([]),  # not provided
        (['']),  # empty
        (['foo']),  # not integer
        (['-1']),  # smaller than zero
        (['0']),  # zero
        (['1.1'])  # floating point
    ]
)
def test_crg_failing(arguments, tmpdir):
    # test candidate reference genome arguments
    cmd_line = ['bin/referenceseeker', '--crg'] + arguments + ['test/db', 'test/data/Salmonella_enterica_CFSAN000189.fasta']
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'arguments',
    [
        ([]),  # not provided
        (['']),  # empty
        (['foo']),  # not integer
        (['-1']),  # smaller than zero
        (['0']),  # zero
        (['1.1'])  # larger than one
    ]
)
def test_ani_failing(arguments, tmpdir):
    # test ANI arguments
    cmd_line = ['bin/referenceseeker', '--ani'] + arguments + ['test/db', 'test/data/Salmonella_enterica_CFSAN000189.fasta']
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'arguments',
    [
        ([]),  # not provided
        (['']),  # empty
        (['foo']),  # not integer
        (['-1']),  # smaller than zero
        (['0']),  # zero
        (['1.1'])  # larger than one
    ]
)
def test_cdna_failing(arguments, tmpdir):
    # test cDNA arguments
    cmd_line = ['bin/referenceseeker', '--conserved-dna'] + arguments + ['test/db', 'test/data/Salmonella_enterica_CFSAN000189.fasta']
    proc = run(cmd_line)
    assert proc.returncode != 0