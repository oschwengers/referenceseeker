
import java.nio.file.*


metadata        = params.metadata
representatives = params.representatives
domain          = params.domain


Channel.fromPath( metadata )
    .splitCsv( skip: 1, sep: '\t'  )
    .filter( { it[15].toLowerCase() == 't' } )
    .map( {
        def acc = it[0] - 'RS_' - 'GB_'
        def orgName = it[16].split(';').last() - 's__'
        def path = acc.substring(0,3) + '/' + acc.substring(4,7) + '/' + acc.substring(7,10) + '/' + acc.substring(10,13)
        def status = it[45].toLowerCase()
        return [ acc, '-', status, orgName, path ]
    } )
    .dump()
    .set { validGenomes }


process sketch {

    tag { "${acc} - ${orgName}" }

    errorStrategy 'ignore'
    maxRetries 3
    conda 'mash=2.3'

    input:
    tuple val(acc), val(taxId), val(status), val(orgName), val(path) from validGenomes

    output:
    tuple val(acc), val(taxId), val(status), val(orgName) into chDbEntries
    file("${acc}.msh") into outMash
    file("${acc}.fna.gz") into outFasta

    publishDir pattern: '*.fna.gz', path: "./${domain}/", mode: 'move'
    publishDir pattern: '*.msh', path: './sketches/',  mode: 'move'

    script:
    """
    gzip -dc ${params.representatives}/${path}/${acc}_genomic.fna.gz > ${acc}
    mash sketch -k 32 -s 10000 ${acc}
    mv ${acc} ${acc}.fna
    gzip ${acc}.fna
    """
}


chDbEntries.map { "${it[0]}\t${it[1]}\t${it[2]}\t${it[3]}" }
    .collectFile( name: 'db.tsv', storeDir: "./${domain}/", newLine: true )
