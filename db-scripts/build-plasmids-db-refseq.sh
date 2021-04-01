
NCBI_PATH='ftp://ftp.ncbi.nlm.nih.gov/genomes'

echo "Download plasmid sequences..."
for i in 1 2 3 4 5 6 7 8 9; do
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.${i}.1.genomic.fna.gz
done

echo "Unzip plasmid sequences..."
zcat plasmid.* > plasmids.fna

nextflow run $REFERENCE_SEEKER_HOME/db-scripts/build-plasmids-db-refseq.nf --plasmids plasmids.fna  ||  { echo "Nextflow failed!"; exit; }

$REFERENCE_SEEKER_HOME/share/mash paste db sketches/*.msh  ||  { echo "Mash failed!"; exit; }

rm -rf work/ .nextflow* sketches/ plasmid.* plasmids.fna

mv db.msh plasmids-refseq/

tar -I pigz -cf plasmids-refseq.tar.gz plasmids-refseq
md5sum plasmids.tar.gz