
import java.nio.file.*


metadata        = params.metadata
representatives = params.representatives
domain          = params.domain


Channel.fromPath( metadata )
    .splitCsv( skip: 1, sep: '\t'  )
    .filter( { it[15].toLowerCase() == 't' } )
    .map( {
        def acc = it[0] - 'RS_' - 'GB_'
	    def status = it[45].toLowerCase()
        def orgName = it[16].split(';').last() - 's__'
        return [ acc, '-', orgName, status ]
    } )
    .set { validGenomes }


process sketch {

    tag { "${acc} - ${orgName}" }

    maxForks 3
    errorStrategy 'ignore'
    maxRetries 3
    conda 'mash=2.3'

    input:
    set val(acc), val(taxId), val(orgName), val(status) from validGenomes

    output:
    set val(acc), val(taxId), val(status), val(orgName) into dbEntries
    file("${acc}.msh") into outMash
    file("${acc}.fna.gz") into outFasta

    publishDir pattern: '*.fna.gz', path: "./${domain}/", mode: 'move'
    publishDir pattern: '*.msh', path: './sketches/',  mode: 'move'

    script:
    """
    gzip -dc ${params.representatives}/${acc}_genomic.fna.gz > ${acc}
    mash sketch -k 32 -s 10000 ${acc}
    mv ${acc} ${acc}.fna
    gzip ${acc}.fna
    """
}


dbEntries.map { "${it[0]}\t${it[1]}\t${it[2]}\t${it[3]}" }
    .collectFile( name: 'db.tsv', storeDir: "./${domain}/", newLine: true )
