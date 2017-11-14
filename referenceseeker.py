#!/usr/bin/env python3

import argparse
import os
import re
import sys
import shutil
import tempfile
import subprocess as sp
from Bio import SeqIO
import multiprocessing as mp


parser = argparse.ArgumentParser( description='Fast determination of finished reference genomes.' )
parser.add_argument( '--db', '-d', required=True, help='ReferenceSeeker database path' )
parser.add_argument( '--genome', '-g', required=True, help='Target draft genome' )
parser.add_argument( '--scaffolds', '-s', action='store_const', default=False, const=True, help='Build scaffolds via MeDuSa (Bosi, Donati et al. 2015) based on detected references' )
parser.add_argument( '--output', '-o', required=False, help='Output fasta file for built scaffolds' )
parser.add_argument( '--cpus', '-c', required=False, type=int, help='Number of cpus to use (default = all available)' )
parser.add_argument( '--verbose', '-v', action='store_const', default=False, const=True, help='Print verbose information' )
args = parser.parse_args()


__MASH_THRESHOLD__ = '0.01'


# check parameters & environment variables
REFERENCE_SEEKER_HOME = os.getenv( 'REFERENCE_SEEKER_HOME', None )
if( REFERENCE_SEEKER_HOME == None ): sys.exit( 'ERROR: REFERENCE_SEEKER_HOME not set!' )
REFERENCE_SEEKER_HOME = os.path.abspath( REFERENCE_SEEKER_HOME )
if( not os.access( REFERENCE_SEEKER_HOME, os.R_OK & os.X_OK ) ): sys.exit( 'ERROR: REFERENCE_SEEKER_HOME ('+REFERENCE_SEEKER_HOME+') not readable/accessible!' )
if( args.verbose ): print( 'REFERENCE_SEEKER_HOME: ' + REFERENCE_SEEKER_HOME )

dbPath = os.path.abspath( args.db )
if( not os.access( dbPath, os.R_OK ) ): sys.exit( 'ERROR: database directory not readable!' )
if( args.verbose ): print( 'db path: ' + dbPath )

genomePath = os.path.abspath( args.genome )
if( not os.access( genomePath, os.R_OK ) ): sys.exit( 'ERROR: genome file not readable!' )
if( args.verbose ): print( 'genome path: ' + genomePath )

if( args.verbose ): print( 'build scaffolds: ' + str(args.scaffolds) )
cwdPath = os.path.abspath( os.getcwd() )
scaffoldsPath = os.path.abspath(args.output) if args.output else cwdPath + '/scaffolds.fna'
if( args.verbose  and  args.scaffolds ): print( 'scaffold path: ' + scaffoldsPath )

noCpus = int(args.cpus) if args.cpus else mp.cpu_count()
if( args.verbose ): print( '# cpus: ' + str(noCpus) )

fhFNULL = open( os.devnull, 'w' )

def compute_ani( refGenome ):
    reference = dbPath + '/' + refGenome['id'] + '.fna'
    # build and write DNA fragemnts of length 1020
    dnaFragmentIdx = 1
    dnaFragments = {}
    tmpDir = tempfile.mkdtemp()
    dnaFragmentsPath = tmpDir + '/dna-fragments.fasta'
    with open( dnaFragmentsPath, 'w' ) as fhDnaFragmentsPath:
        for record in SeqIO.parse( refGenome['query'], 'fasta' ):
            sequence = record.seq
            while( len(sequence) > 1020 ):
                dnaFragment = sequence[:1020]
                fhDnaFragmentsPath.write( '>' + str(dnaFragmentIdx) + '\n' )
                fhDnaFragmentsPath.write( str(dnaFragment) + '\n' )
                dnaFragments[ dnaFragmentIdx ] = { 'id': dnaFragmentIdx, 'length': len( dnaFragment ) }
                sequence = sequence[1020:]
                dnaFragmentIdx += 1
            dnaFragment = sequence
            fhDnaFragmentsPath.write( '>' + str(dnaFragmentIdx) + '\n' )
            fhDnaFragmentsPath.write( str(dnaFragment) + '\n' )
            dnaFragments[ dnaFragmentIdx ] = { 'id': dnaFragmentIdx, 'length': len( dnaFragment ) }
            sequence = sequence[1020:]
            dnaFragmentIdx += 1

    # perform global alignments via nucmer
    sp.check_call( [ REFERENCE_SEEKER_HOME + '/share/mummer/nucmer',
        '--threads=1',
        reference,
        dnaFragmentsPath ],
        cwd = tmpDir,
        stdout = fhFNULL,
        stderr = sp.STDOUT
    )
    filteredDeltaPath = tmpDir + '/out-filtered.delta'
    with open( filteredDeltaPath, 'w' ) as fhFilteredDeltaPath:
        sp.check_call( [ REFERENCE_SEEKER_HOME+'/share/mummer/delta-filter',
            '-q',
            'out.delta' ],
            cwd = tmpDir,
            stdout = fhFilteredDeltaPath,
            stderr = sp.STDOUT
        )

    # parse nucmer output
    dnaFragment = None
    dnaFragementMatches = []
    with open( filteredDeltaPath, 'r' ) as fhFilteredDeltaPath:
        for line in fhFilteredDeltaPath:
            line = line.rstrip()
            if( line[0] == '>' ):
                dnaFragment = dnaFragments[ int( line.split( ' ' )[1] ) ]
            elif( dnaFragment != None ):
                cols = line.split( ' ' )
                if( len(cols) == 7 ):
                    dnaFragment[ 'alignmentLength' ] = abs( int( cols[3] ) - int( cols[2] ) ) + 1  # abs( qStop - qStart ) + 1
                    dnaFragment[ 'noNonIdentities' ] = int( cols[4] )  # number of non-identities
                    dnaFragementMatches.append( dnaFragment )

    # calc % conserved DNA
    alignmentSum = 0
    for fm in dnaFragementMatches:
        relAlignmentLength = float(fm['alignmentLength'] - fm['noNonIdentities']) / float(fm['length'])
        if( relAlignmentLength > 0.9 ):
            alignmentSum += fm['alignmentLength']
    genomeLength = 0
    for fm in dnaFragments.values():
        genomeLength += fm['length']
    conservedDna = float(alignmentSum) / float(genomeLength)

    # calc average nucleotide identity
    aniMatches = 0
    niSum = 0.0
    for fm in dnaFragementMatches:
        if( ((float(fm['alignmentLength']-fm['noNonIdentities'])/float(fm['length'])) > 0.3 )  and  ( (float(fm['alignmentLength'])/float(fm['length'])) >= 0.7) ):
            niSum += float(fm['alignmentLength']-fm['noNonIdentities'])/float(fm['alignmentLength'])
            aniMatches += 1
    ani = 0 if aniMatches == 0 else niSum / float(aniMatches)

    shutil.rmtree( tmpDir )

    if( args.verbose ): print( '\tref: %s\tANI: % 2.2f\tcons. DNA: % 2.2f' % (reference.split('/')[-1], ani*100, conservedDna*100) )
    refGenome[ 'ani' ] = ani
    refGenome[ 'conservedDna' ] = conservedDna
    return refGenome




# Calculate genome distances via Mash
if( args.verbose ): print( 'calc genome distances...' )
mashResultPath = cwdPath + '/mash.out'
with open( mashResultPath, 'w' ) as fhMashResultPath:
    sp.check_call( [ REFERENCE_SEEKER_HOME + '/share/mash/mash',
        'dist',
        '-d', __MASH_THRESHOLD__,
        dbPath + '/db.msh',
        genomePath ],
        stdout = fhMashResultPath,
        stderr = sp.STDOUT
    )


# Extract hits and store dist
accessionIds = []
mashDistances = {}
with open( mashResultPath, 'r' ) as fhMashResultPath:
    for line in fhMashResultPath:
        cols = line.rstrip().split()
        accessionIds.append( cols[0] )
        mashDistances[ cols[0] ] = float( cols[2] )
os.remove( mashResultPath )


if( args.verbose ): print( '\tfound ' + str(len(accessionIds)) + ' potential reference genome(s)' )
# Reduce Mash output to best 100 hits
if( len( accessionIds) > 100 ):
    if( args.verbose ): print( '\treduce to best 100 hits...' )
    tmpAccessionIds = sorted( accessionIds, key=lambda k: mashDistances[ k ] )
    accessionIds = tmpAccessionIds[:100]


# Get assemblies from RefSeq by accessions
if( args.verbose ): print( 'lookup corresponding reference files...' )
refGenomes = []
with open( dbPath + '/db.tsv', 'r' ) as fhDbPath:
    for line in fhDbPath:
        if( line[0] != '#' ):
            cols = line.strip().split( '\t' )
            accessionId = cols[0]
            if( accessionId in accessionIds ):
                refGenomes.append( { 'id': accessionId, 'tax': cols[1], 'name': cols[2], 'query': genomePath, 'dist': mashDistances[ accessionId ] } )


# Copy genomes, extract them and build ANI
if( args.verbose ): print( 'compute ANIs...' )
pool = mp.Pool( noCpus )
results = pool.map( compute_ani, refGenomes )
pool.close()
pool.join()


# filter and sort results
tmp_results = []
for result in results :
    if( (result['conservedDna'] >= 0.69)  and (result['ani'] >= 0.98) ):
        tmp_results.append( result )
results = sorted( tmp_results, key=lambda k: -(k['ani']*k['conservedDna']) )


# print results to STDOUT
if( args.verbose ): print( 'ANIs:' )
for result in results:
    print( '%s\t% 2.2f\t% 2.2f\t%s\t%s' % (result['id'], result['ani']*100, result['conservedDna']*100, result['tax'], result['name'] ) )
    #print( '%s\t%s\t% 2.2f\t% 2.2f\t%s\t%s' % (result['id'], result['dist'], result['ani']*100, result['conservedDna']*100, result['tax'], result['name'] ) )


# create scaffolds based on draft genome and detected references
if( args.scaffolds ):
    if( args.verbose ): print( 'create scaffolds...' )
    # copy first 20 references to tmp file
    tmpDir = tempfile.mkdtemp()
    if( len(results) > 20 ): results = results[:20]
    for result in results:
        os.symlink( dbPath + '/' + result['id'] + '.fna', tmpDir + '/' + result['id'] + '.fna' )
    os.putenv( 'PATH', REFERENCE_SEEKER_HOME+'/share/mummer/:' + os.getenv( 'PATH' ) )
    sp.check_call( [ 'java', '-jar', REFERENCE_SEEKER_HOME+'/share/medusa/medusa.jar',
        '-f', tmpDir,
        '-i', genomePath, # assembled alignments
        '-o', scaffoldsPath, # new name
        '-v', # verbose flag
        '-random', '1000', # use 1000 random rounds to find the best scaffolds
        '-scriptPath', REFERENCE_SEEKER_HOME+'/share/medusa/medusa_scripts' # path to medusa script directory
        ],
        cwd = cwdPath,
        stdout = fhFNULL,
        stderr = sp.STDOUT
    )
    shutil.rmtree( tmpDir )

    # parse MeDuSa output and print results
    medusaResultPath = genomePath + '_SUMMARY'
    if( os.path.isfile( medusaResultPath ) ):
        with open( medusaResultPath, 'r' ) as medusaResultFile:
            medusaResult = medusaResultFile.read()
            mg = re.search( 'singletons = (\d+), multi-contig scaffold = (\d+)', medusaResult )
            noContigs   = mg.group(1)
            noScaffolds = mg.group(2)
            mg = re.search( 'from (\d+) initial fragments', medusaResult )
            noInitContigs = mg.group(1)
            if( args.verbose ): print( '\tordered ' + noInitContigs + ' initial contigs' )
            if( args.verbose ): print( '\t# new scaffolds:' + noScaffolds )
            if( args.verbose ): print( '\t# new contigs:' + noContigs )
        os.remove( medusaResultPath )
