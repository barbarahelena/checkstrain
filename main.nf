#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    checkstrain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pull a bacterial genome assembly from NCBI, annotate it with Bakta, and
    pull out the osmoadaptation genes that osmotool profiles.
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { CHECKSTRAIN } from './workflows/checkstrain'

def helpMessage() {
    log.info """
    Usage:
      nextflow run main.nf --input samplesheet.csv --outdir results [-profile docker]

    Required params:
      --input   CSV with header 'name,accession' -- one row per genome, e.g.
                    name,accession
                    Escherichia coli K-12,GCF_000005845.2
                'accession' is the NCBI assembly accession (GCF_/GCA_...).
      --outdir  Directory to write results to.

    Optional params:
      --bakta_database     Path to a pre-built Bakta database (skips BAKTA_BAKTADBDOWNLOAD).
      --osmo_db            Path to a pre-built osmo_refdb release (skips OSMOTOOL_DOWNLOAD_DB).
      --osmo_db_release    osmo_refdb release to download when --osmo_db is not set (default: '${params.osmo_db_release}').
      --publish_dir_mode   File publishing mode, e.g. 'copy' or 'symlink' (default: '${params.publish_dir_mode}').
    """.stripIndent()
}

workflow {

    if ( params.help ) {
        helpMessage()
        exit 0
    }
    if ( !params.input ) {
        error "Missing required param --input (samplesheet with columns: name,accession). Run with --help for usage."
    }
    if ( !params.outdir ) {
        error "Missing required param --outdir. Run with --help for usage."
    }

    ch_samplesheet = Channel
        .fromPath( params.input, checkIfExists: true )
        .splitCsv( header: true )
        .map { row ->
            if ( !row.name || !row.accession ) {
                error "Every samplesheet row needs 'name' and 'accession' columns, got: ${row}"
            }
            def id = row.name.trim().replaceAll(/[^A-Za-z0-9]+/, '_').replaceAll(/^_+|_+$/, '').toLowerCase()
            def meta = [ id: id, name: row.name.trim(), accession: row.accession.trim() ]
            tuple( meta, row.accession.trim() )
        }

    CHECKSTRAIN ( ch_samplesheet )
}
