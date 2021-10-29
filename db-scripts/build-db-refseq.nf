
import java.nio.file.*


assemblySummary = params.ass_sum
ncbiPath        = params.ncbiPath
domain          = params.domain


Channel.fromPath( assemblySummary )
    .splitCsv( skip: 2, sep: '\t'  )
    .filter( { (it[11].toLowerCase() == 'complete genome')  ||  (it[4].toLowerCase() == 'representative genome')  ||  (it[4].toLowerCase() == 'reference genome') } )
    .map( {
        def species = it[7]
        def strain  = it[8] - 'strain='
	def status  = it[11].split(' ')[0].toLowerCase()
        if( species.contains( strain ) )
            return [ it[0], it[5], species, status, it[19] - 'https://ftp.ncbi.nlm.nih.gov/genomes/' ]
        else
            return [ it[0], it[5], "${species} ${strain}", status, it[19] - 'https://ftp.ncbi.nlm.nih.gov/genomes/' ]
    } )
    .set { chValidGenomes }


process download {

    tag { "${acc} - ${orgName}" }

    executor 'local'
    maxForks 5
    errorStrategy 'ignore'
    maxRetries 3

    input:
    tuple val(acc), val(taxId), val(orgName), val(status), val(path) from chValidGenomes

    output:
    tuple val(acc), val(taxId), val(orgName), val(status), path("${acc}.gz") into chDownloadedGenomes

    script:
    """
    wget -O ${acc}.gz ${ncbiPath}/${path}/${path.split('/').last()}_genomic.fna.gz
    """
}


process sketch {

    tag { "${acc} - ${orgName}" }

    errorStrategy 'ignore'
    maxRetries 3
    conda 'mash=2.3'

    input:
    tuple val(acc), val(taxId), val(orgName), val(status), path("${acc}.gz") from chDownloadedGenomes

    output:
    tuple val(acc), val(taxId), val(status), val(orgName) into chDbEntries
    file("${acc}.msh") into outMash
    file("${acc}.fna.gz") into outFasta

    publishDir pattern: '*.fna.gz', path: "./${domain}-refseq/", mode: 'move'
    publishDir pattern: '*.msh', path: './sketches/',  mode: 'move'

    script:
    """
    gunzip -c ${acc}.gz > ${acc}
    mash sketch -k 32 -s 10000 ${acc}
    mv ${acc} ${acc}.fna
    gzip ${acc}.fna
    """
}


chDbEntries.map { "${it[0]}\t${it[1]}\t${it[2]}\t${it[3]}" }
    .collectFile( name: 'db.tsv', storeDir: "./${domain}-refseq/", newLine: true )
