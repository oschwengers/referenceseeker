nextflow.enable.dsl=2

import java.nio.file.*


// params.metadata
// params.representatives
// params.domain


process sketch {

    tag { "${acc} - ${orgName}" }

    errorStrategy 'ignore'
    maxRetries 3
    conda 'mash=2.3'
    container 'quay.io/biocontainers/mash:2.3--hd3113c8_6'

    input:
    tuple val(acc), val(taxId), val(status), val(orgName), val(path)

    output:
    tuple val(acc), val(taxId), val(status), val(orgName), emit: results
    path("${acc}.msh")
    path("${acc}.fna.gz")

    publishDir pattern: "${acc}.msh", path: './sketches/',  mode: 'move'
    publishDir pattern: "${acc}.fna.gz", path: "./${params.domain}/", mode: 'move'
    
    script:
    """
    gzip -dc ${params.representatives}/${path}/${acc}_genomic.fna.gz > ${acc}
    mash sketch -k 32 -s 10000 ${acc}
    mv ${acc} ${acc}.fna
    gzip ${acc}.fna
    """
}


workflow {

    processedGenomes = Channel.fromPath( params.metadata )
        | splitCsv( skip: 1, sep: '\t'  )
        | filter( { it[18].toLowerCase() == 't' } )
        | map( {
            def acc = it[0] - 'RS_' - 'GB_'
            def orgName = it[65].split(';').last() - 's__'
            def path = acc.substring(0,3) + '/' + acc.substring(4,7) + '/' + acc.substring(7,10) + '/' + acc.substring(10,13)
            def status = it[48].toLowerCase()
            return [ acc, '-', status, orgName, path ]
        } )
        | sketch

    processedGenomes.results
        | map( { "${it[0]}\t${it[1]}\t${it[2]}\t${it[3]}\n" } )
        | collectFile( name: 'db.tsv', sort: false, tempDir: "${workDir}/", storeDir: "./${params.domain}/" )

}
