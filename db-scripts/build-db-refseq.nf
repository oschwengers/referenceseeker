nextflow.enable.dsl=2

import java.nio.file.*


// params.ass_sum
// params.ncbiPath
// params.domain


process download {

    tag { "${acc} - ${orgName}" }

    executor 'local'
    maxForks 5
    errorStrategy 'ignore'
    maxRetries 3

    input:
    tuple val(acc), val(taxId), val(orgName), val(status), val(path)

    output:
    tuple val(acc), val(taxId), val(orgName), val(status), path("${acc}.gz")

    script:
    """
    wget -O ${acc}.gz ${params.ncbiPath}/${path}/${path.split('/').last()}_genomic.fna.gz
    """
}


process sketch {

    tag { "${acc} - ${orgName}" }

    errorStrategy 'ignore'
    maxRetries 3
    conda 'mash=2.3'
    container 'quay.io/biocontainers/mash:2.3--hd3113c8_6'

    input:
    tuple val(acc), val(taxId), val(orgName), val(status), path("${acc}.gz")

    output:
    tuple val(acc), val(taxId), val(status), val(orgName), emit: results
    file("${acc}.msh")
    file("${acc}.fna.gz")

    publishDir pattern: "${acc}.msh", path: './sketches/',  mode: 'move'
    publishDir pattern: "${acc}.fna.gz", path: "./${params.omain}-refseq/", mode: 'move'
    
    script:
    """
    gunzip -c ${acc}.gz > ${acc}
    mash sketch -k 32 -s 10000 ${acc}
    mv ${acc} ${acc}.fna
    gzip ${acc}.fna
    """
}


workflow {

    processedGenomes = Channel.fromPath( params.assemblySummary )
        | splitCsv( skip: 2, sep: '\t'  )
        | filter( { (it[11].toLowerCase() == 'complete genome')  ||  (it[4].toLowerCase() == 'representative genome')  ||  (it[4].toLowerCase() == 'reference genome') } )
        | map( {
            def species = it[7]
            def strain  = it[8] - 'strain='
	        def status  = it[11].split(' ')[0].toLowerCase()
            if( species.contains( strain ) )
                return [ it[0], it[5], species, status, it[19] - 'https://ftp.ncbi.nlm.nih.gov/genomes/' ]
            else
                return [ it[0], it[5], "${species} ${strain}", status, it[19] - 'https://ftp.ncbi.nlm.nih.gov/genomes/' ]
        } )
        | download
        | sketch

    processedGenomes.results
        | map { "${it[0]}\t${it[1]}\t${it[2]}\t${it[3]}" }
        | collectFile( name: 'db.tsv', storeDir: "./${params.domain}-refseq/", newLine: true )

}
