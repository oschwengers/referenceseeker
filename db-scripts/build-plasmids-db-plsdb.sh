
echo "Download plasmid sequences..."
wget -O plasmids-plsdb.bz2 https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/plsdb.fna.bz2

echo "Unzip plasmid sequences..."
bunzip2 plsdb.fna.bz2
nextflow run $REFERENCE_SEEKER_HOME/db-scripts/build-plasmids-db-plsdb.nf --plasmids plsdb.fna  ||  { echo "Nextflow failed!"; exit; }

find sketches/ -type f -name '*.msh' -exec realpath {} + > sketches.fof
mash paste -l db sketches.fof  ||  { echo "Mash failed!"; exit; }

rm -rf work/ .nextflow* sketches/ sketches.fof plsdb.fna

mv db.msh plasmids-plsdb/

tar -I pigz -cf plasmids-plsdb.tar.gz plasmids-plsdb
md5sum plasmids-plsdb.tar.gz