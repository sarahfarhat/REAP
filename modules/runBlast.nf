#!/usr/bin/env nextflow

process runBlast {
    tag "$species"

    input:
    tuple path(fasta), val(species)
    path(db_files)
    val blast_program
    val evalue

    publishDir '01_blastResults', mode: 'copy'

    output:
    tuple path("${species}.repeats.blastout6"), val(species)

    script:
    // Strip any BLAST db extension (works for both nucl and prot databases)
    def db_prefix = db_files[0].getName().replaceAll(/\.(phr|pin|psq|pdb|pto|psi|ptf|nhr|nin|nsq|ndb|not|nto|ntf)$/, '')

    """
    ${blast_program} -db ${db_prefix} -query ${fasta} -out ${species}.repeats.blastout6 \
        -evalue ${evalue} -outfmt 6 -num_threads ${task.cpus}
    """
}
