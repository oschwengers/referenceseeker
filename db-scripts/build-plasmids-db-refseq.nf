nextflow.enable.dsl=2

import java.nio.file.*


// params.plasmids


Channel.fromPath(plasmids)
    .into { fastaFiles; fastaEntries }

fastaFiles
    .splitFasta( by: 1, file: true)
    .set { fastaFiles }

fastaEntries
    .splitFasta( record: [id: true, desc: true ] )
    .map( {
        def plasmidName = it.desc - ', complete sequence'
	    return [ it.id, '', plasmidName, 'complete' ]
    } )
    .into { validGenomes; dbEntries }


process sketch {

    tag { "${acc} - ${orgName}" }

    errorStrategy 'ignore'
    maxRetries 3
    conda 'mash=2.3'
    container 'quay.io/biocontainers/mash:2.3--hd3113c8_6'

    input:
    file(sequence)
    set val(acc), val(taxId), val(orgName), val(status)

    output:
    file("${acc}.msh")
    file("${acc}.fna.gz")

    publishDir pattern: '*.msh', path: './sketches/', mode: 'move'
    publishDir pattern: '*.fna.gz', path: "./plasmids-refseq/", mode: 'move'
    
    script:
    """
    mv ${sequence} ${acc}
    mash sketch -k 32 -s 1000 ${acc}
    cp -L ${acc} ${acc}.fna
    gzip ${acc}.fna
    """
}

workflow {
    
    plasmidsFasta = Channel.fromPath(params.plasmids)

    plasmidSequences = plasmidsFasta
        | splitFasta( by: 1, file: true)
    
    plasmidRecords = plasmidsFasta
        | splitFasta( record: [id: true, desc: true ] )
        | map( {
            def plasmidName = it.desc - ', complete sequence'
	        return [ it.id, '', plasmidName, 'complete' ]
        } )

    sketch( plasmidSequences, plasmidRecords)

    plasmidRecords
        | map( { "${it[0]}\t${it[1]}\t${it[2]}\t${it[3]}" } )
        | collectFile( name: 'db.tsv', storeDir: './plasmids-refseq/', newLine: true )
}
