nextflow $REFERENCE_SEEKER_HOME/build-db.nf --domain $1  ||  { echo "Nextflow failed!"; exit; }
$REFERENCE_SEEKER_HOME/share/mash/mash paste db sketches/*.msh  ||  { echo "Mash failed!"; exit; }
rm -rf work/ .nextflow* sketches/
mv db.msh $1/
