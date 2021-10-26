
import java.nio.file.*


plasmids        = params.plasmids


Channel.fromPath(plasmids)
    .into { fastaFiles; fastaEntries }

fastaFiles
    .splitFasta( by: 1, file: true)
    .set { fastaFiles }

fastaEntries
    .splitFasta( record: [id: true, desc: true ] )
    .map( {
	    return [ it.id, '', it.desc, 'complete' ]
    } )
    .into { validGenomes; dbEntries }


process sketch {

    tag { "${acc} - ${orgName}" }

    errorStrategy 'ignore'
    maxRetries 3
    conda 'mash=2.3'

    input:
    file(sequence) from fastaFiles
    set val(acc), val(taxId), val(orgName), val(status) from validGenomes

    output:
    file("${acc}.msh") into outMash
    file("${acc}.fna.gz") into outFasta

    publishDir pattern: '*.fna.gz', path: "./plasmids-plsdb/", mode: 'move'
    publishDir pattern: '*.msh', path: './sketches/', mode: 'move'

    script:
    """
    mv ${sequence} ${acc}
    mash sketch -k 32 -s 1000 ${acc}
    cp -L ${acc} ${acc}.fna
    gzip ${acc}.fna
    """
}


dbEntries.map { "${it[0]}\t${it[1]}\t${it[2]}\t${it[3]}" }
    .collectFile( name: 'db.tsv', storeDir: './plasmids-plsdb/', newLine: true )
