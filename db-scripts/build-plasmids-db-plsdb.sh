
echo "Download plasmid sequences..."
wget -O plasmids-plsdb.zip https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/?zip

echo "Unzip plasmid sequences..."
unzip plasmids-plsdb.zip
blastdbcmd -entry all -db plsdb.fna -out plsdb.fna
nextflow run $REFERENCE_SEEKER_HOME/db-scripts/build-plasmids-db-plsdb.nf --plasmids plsdb.fna  ||  { echo "Nextflow failed!"; exit; }

$REFERENCE_SEEKER_HOME/share/mash paste db sketches/*.msh  ||  { echo "Mash failed!"; exit; }

rm -rf work/ .nextflow* sketches/ plsdb* README.md

mv db.msh plasmids-plsdb/

tar -I pigz -cf plasmids-plsdb.tar.gz plasmids-plsdb
md5sum plasmids-plsdb.tar.gz