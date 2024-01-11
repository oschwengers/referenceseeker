
echo "Download GTDB representatives genomes..."
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz
tar -xzf gtdb_genomes_reps.tar.gz

DOMAIN=""
METADATA=""
getopts "ba" opt;
case "$opt" in
h|\?)
    echo -e "Select one of the following options:\n\t-b (bacteria)\n\t-a (archaea)\n"
    exit 0
    ;;
b)  DOMAIN="bacteria-gtdb"
    METADATA="https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv.gz"
    ;;
a)  DOMAIN="archaea-gtdb"
    METADATA="https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_metadata.tsv.gz"
    ;;
\?)
    echo "Invalid option: -$OPTARG" >&2
    exit 1
    ;;
esac

echo "Download $DOMAIN database..."
wget $METADATA
mv $METADATA metadata.tsv.gz
gunzip metadata.tsv.gz

DATA_PATH=$(readlink -f gtdb_genomes_reps_r*/database)
nextflow run $REFERENCE_SEEKER_HOME/db-scripts/build-db-gtdb.nf --metadata ./metadata.tsv --representatives $DATA_PATH --domain $DOMAIN  ||  { echo "Nextflow failed!"; exit; }

mash paste db sketches/*.msh  ||  { echo "Mash failed!"; exit; }

rm -rf work/ .nextflow* sketches/ *.tsv *_metadata.tar.gz gtdb_genomes_reps.tar.gz

mv db.msh $DOMAIN/

tar -I pigz -cf "$DOMAIN.tar.gz" $DOMAIN
md5sum "$DOMAIN.tar.gz"