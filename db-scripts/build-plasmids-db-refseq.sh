
echo "Download plasmid sequences..."
for i in 1 2 3 4 5; do
    for k in 1 2; do
        wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.${i}.${k}.genomic.fna.gz
    done
done

echo "Unzip plasmid sequences..."
zcat plasmid.* > plasmids.fna

nextflow run $REFERENCE_SEEKER_HOME/db-scripts/build-plasmids-db-refseq.nf --plasmids plasmids.fna  ||  { echo "Nextflow failed!"; exit; }

find sketches/ -type f -name '*.msh' -exec realpath {} + > sketches.fof
mash paste -l db sketches.fof  ||  { echo "Mash failed!"; exit; }

rm -rf work/ .nextflow* sketches/ plasmid.* sketches.fof plasmids.fna

mv db.msh plasmids-refseq/

tar -I pigz -cf plasmids-refseq.tar.gz plasmids-refseq
md5sum plasmids.tar.gz