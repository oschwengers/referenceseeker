
NCBI_PATH='ftp://ftp.ncbi.nlm.nih.gov/genomes'

#OPTIND=1
DOMAIN=''
getopts "bavf" opt;
case "$opt" in
h|\?)
    echo "Select one of the following options:\n\t-b (bacteria)\n\t-a (archaea)\n\t-v (viral)\n\t-f (fungi)"
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
\?)
    echo "Invalid option: -$OPTARG" >&2
    exit 1
    ;;
esac

echo "Download $DOMAIN database..."
wget -O assembly_summary.txt ${NCBI_PATH}/refseq/${DOMAIN}/assembly_summary.txt

nextflow $REFERENCE_SEEKER_HOME/build-db.nf --ass_sum ./assembly_summary.txt --ncbiPath $NCBI_PATH --domain $DOMAIN  ||  { echo "Nextflow failed!"; exit; }

$REFERENCE_SEEKER_HOME/share/mash/mash paste db sketches/*.msh  ||  { echo "Mash failed!"; exit; }

rm -rf work/ .nextflow* sketches/ assembly_summary.txt

mv db.msh $DOMAIN/
