
NCBI_PATH='https://ftp.ncbi.nlm.nih.gov/genomes'

DOMAIN=''
getopts "bavfp" opt;
case "$opt" in
h|\?)
    echo -e "Select one of the following options:\n\t-b (bacteria)\n\t-a (archaea)\n\t-v (viral)\n\t-f (fungi)\n\t-p (protozoa)"
    exit 0
    ;;
b)  DOMAIN='bacteria'
    ;;
a)  DOMAIN='archaea'
    ;;
v)  DOMAIN='viral'
    ;;
f)  DOMAIN='fungi'
    ;;
p)  DOMAIN='protozoa'
    ;;
\?)
    echo "Invalid option: -$OPTARG" >&2
    exit 1
    ;;
esac

echo "Download $DOMAIN database..."
wget -O assembly_summary.txt ${NCBI_PATH}/refseq/${DOMAIN}/assembly_summary.txt

nextflow run $REFERENCE_SEEKER_HOME/db-scripts/build-db-refseq.nf --ass_sum ./assembly_summary.txt --ncbiPath $NCBI_PATH --domain $DOMAIN  ||  { echo "Nextflow failed!"; exit; }

$REFERENCE_SEEKER_HOME/share/mash paste db sketches/*.msh  ||  { echo "Mash failed!"; exit; }

rm -rf work/ .nextflow* sketches/ assembly_summary.txt

mv db.msh "$DOMAIN-refseq"/

tar -I pigz -cf "$DOMAIN-refseq.tar.gz" "$DOMAIN-refseq"
md5sum "$DOMAIN-refseq.tar.gz"
