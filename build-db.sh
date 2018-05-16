nextflow $REFERENCE_SEEKER_HOME/build-db.nf --domain $1  ||  { echo "Nextflow failed!"; exit; }
$REFERENCE_SEEKER_HOME/share/mash/mash paste db $1/*.msh  ||  { echo "Mash failed!"; exit; }
rm -rf work/ .nextflow* $1/*.msh
mv db.msh $1/
