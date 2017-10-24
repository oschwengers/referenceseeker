
import java.nio.file.*


domain   = params.domain
ncbiPath = 'ftp://ftp.ncbi.nlm.nih.gov/genomes'


Channel.fromPath( "${ncbiPath}/refseq/${domain}/assembly_summary.txt" )
    .splitCsv( skip: 2, sep: '\t'  )
    .filter( { it[11] == 'Complete Genome' } )
    .map( {
        def species = it[7]
        def strain  = it[8] - 'strain='
        if( species.contains( strain ) )
            return [ it[0], it[5], species, it[19] - 'ftp://ftp.ncbi.nlm.nih.gov/genomes/' ]
        else
            return [ it[0], it[5], "${species} ${strain}", it[19] - 'ftp://ftp.ncbi.nlm.nih.gov/genomes/' ]
    } )
    .set { validGenomes }


process sketch {

    maxForks 3
    maxRetries 3
    tag { "${acc} - ${orgName}" }

    input:
    set val(acc), val(taxId), val(orgName), val(path) from validGenomes

    output:
    set val(acc), val(taxId), val(orgName) into dbEntries
    file("${acc}.msh") into outMash
    file("${acc}.fna") into outFasta

    publishDir path: "./${domain}/", mode: 'link'

    script:
    """
    wget -O ${acc}.gz ${ncbiPath}/${path}/${path.split('/').last()}_genomic.fna.gz
    gunzip ${acc}.gz
    $REFERENCE_SEEKER_HOME/share/mash/mash sketch -k 32 -s 10000 ${acc}
    mv ${acc} ${acc}.fna
    """
}


dbEntries.map { "${it[0]}\t${it[1]}\t${it[2]}" }
    .collectFile( name: 'db.tsv', storeDir: "./${domain}/", newLine: true )
