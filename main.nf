#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { sfAssign } from "${baseDir}/modules/sfAssign.nf"
include { runBlast } from "${baseDir}/modules/runBlast.nf"

process makeblastdb {
    tag "Creating BLAST database"

    input:
    path fasta_file
    val  dbtype

    output:
    path("repeats.*")

    script:
    """
    makeblastdb -in ${fasta_file} -dbtype ${dbtype} -out repeats
    """
}

// Derive BLAST program from query/db type combination
def blastProgram(query_type, db_type) {
    if      (query_type == "nucl" && db_type == "nucl") return "blastn"
    else if (query_type == "prot" && db_type == "prot") return "blastp"
    else if (query_type == "prot" && db_type == "nucl") return "tblastn"
    else if (query_type == "nucl" && db_type == "prot") return "blastx"
    else error "Invalid combination: query_type=${query_type}, db_type=${db_type}. " +
               "Accepted values are 'nucl' or 'prot'."
}

def blast_program    = blastProgram(params.query_type, params.db_type)
def makeblastdb_type = (params.db_type == "prot") ? "prot" : "nucl"

log.info """\
        =========================================
        Repeat annotation PIPELINE
        =========================================
        Repeat database   : ${params.repeatdb}
        DB type           : ${params.db_type}
        Query type        : ${params.query_type}
        BLAST program     : ${blast_program}
        Input table       : ${params.inputTable}
        E-value threshold : ${params.blastevalue}
        Separator         : ${params.sep}
        =========================================
        """.stripIndent()

Channel
    .fromPath(params.inputTable)
    .splitText()
    .filter { line -> line.trim() }
    .map { line ->
        def parts = line.trim().split(/\s+/)
        if (parts.size() < 2) {
            error "Malformed line in inputTable: '${line}'"
        }
        def fasta_path = parts[0]
        def species    = parts[1]
        return [file(fasta_path), species]
    }
    .set { sequences_species_ch }

workflow {
    makeblastdb(params.repeatdb, makeblastdb_type)

    makeblastdb.out
        .collect()
        .set { db_files_ch }

    runBlast(sequences_species_ch, db_files_ch, blast_program, params.blastevalue)

    runBlast.out
        .map  { blast_file, species -> tuple(blast_file, species.trim()) }
        .join (sequences_species_ch, by: 1)
        .map  { species, blast_file, fasta_file -> tuple(blast_file, fasta_file, species) }
        .set  { sf_assign_input_ch }

    sfAssign(sf_assign_input_ch, params.sep)
}
