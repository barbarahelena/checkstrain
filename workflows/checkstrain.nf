/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { NCBI_DATASETS_DOWNLOAD              } from '../modules/local/ncbi_datasets_download'
include { BAKTA_BAKTADBDOWNLOAD               } from '../modules/local/bakta/baktadbdownload'
include { BAKTA_BAKTA                         } from '../modules/local/bakta/bakta'
include { OSMOTOOL_DOWNLOAD_DB                } from '../modules/local/osmotool/download_db'
include { OSMOTOOL_ANNOTATE                   } from '../modules/local/osmotool/annotate'
include { CHECKSTRAIN_SELECT_GENES_OF_INTEREST } from '../modules/local/checkstrain/select_genes_of_interest'
include { CHECKSTRAIN_COMPARE_BAKTA_OSMOTOOL   } from '../modules/local/checkstrain/compare_bakta_osmotool'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CHECKSTRAIN {

    take:
    ch_accessions   // channel: [ meta, accession ]  (meta.id/name/accession from --input)

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Pull the assembly (genome + gff) for each accession from NCBI
    //
    NCBI_DATASETS_DOWNLOAD ( ch_accessions )
    ch_versions = ch_versions.mix(NCBI_DATASETS_DOWNLOAD.out.versions)
    ch_fasta    = NCBI_DATASETS_DOWNLOAD.out.fna

    //
    // Bakta database: reuse a pre-built one if given, otherwise download it once
    // and broadcast it to every sample (.collect() turns the singleton output
    // into a value channel so it can pair with N samples, not just 1).
    //
    if ( params.bakta_database ) {
        ch_bakta_db = Channel.value( file(params.bakta_database, checkIfExists: true) )
    } else {
        BAKTA_BAKTADBDOWNLOAD ()
        ch_versions = ch_versions.mix(BAKTA_BAKTADBDOWNLOAD.out.versions)
        ch_bakta_db = BAKTA_BAKTADBDOWNLOAD.out.db.collect()
    }

    //
    // MODULE: Annotate the downloaded genome with Bakta
    //
    BAKTA_BAKTA (
        ch_fasta,
        ch_bakta_db,
        [],
        []
    )
    ch_versions = ch_versions.mix(BAKTA_BAKTA.out.versions)

    //
    // osmotool reference database: same reuse-or-download-once pattern as Bakta's.
    //
    if ( params.osmo_db ) {
        ch_osmo_db = Channel.value( file(params.osmo_db, checkIfExists: true) )
    } else {
        OSMOTOOL_DOWNLOAD_DB ( params.osmo_db_release )
        ch_versions = ch_versions.mix(OSMOTOOL_DOWNLOAD_DB.out.versions)
        ch_osmo_db = OSMOTOOL_DOWNLOAD_DB.out.db.collect()
    }

    //
    // MODULE: Profile the osmoadaptation gene panel on the same genome.
    // This defines the "genes of interest": osmotool's gene_counts.tsv always
    // lists its full family panel (count 0 for unmatched families), so it
    // doubles as the panel definition, not just a per-genome result.
    //
    OSMOTOOL_ANNOTATE (
        ch_fasta,
        ch_osmo_db
    )
    ch_versions = ch_versions.mix(OSMOTOOL_ANNOTATE.out.versions)

    //
    // MODULE: Cross-reference Bakta's annotation against osmotool's panel to
    // pull out just the osmoadaptation genes, keeping Bakta's richer fields
    // (locus tag, product, coordinates) alongside osmotool's per-family counts.
    //
    ch_select_input = BAKTA_BAKTA.out.tsv
        .join( OSMOTOOL_ANNOTATE.out.counts )

    CHECKSTRAIN_SELECT_GENES_OF_INTEREST ( ch_select_input )
    ch_versions = ch_versions.mix(CHECKSTRAIN_SELECT_GENES_OF_INTEREST.out.versions)

    //
    // MODULE: Cross-check Bakta's filtered gene calls against osmotool's own
    // coordinate-level hits for the same genome (matched by gene name, since
    // Bakta renumbers contigs while osmotool keeps the input FASTA's original
    // contig names -- the two aren't coordinate-comparable).
    //
    ch_compare_input = CHECKSTRAIN_SELECT_GENES_OF_INTEREST.out.tsv
        .join( OSMOTOOL_ANNOTATE.out.coordinates )

    CHECKSTRAIN_COMPARE_BAKTA_OSMOTOOL ( ch_compare_input )
    ch_versions = ch_versions.mix(CHECKSTRAIN_COMPARE_BAKTA_OSMOTOOL.out.versions)

    //
    // Collate software versions
    //
    ch_versions
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'checkstrain_software_versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    genes_of_interest = CHECKSTRAIN_SELECT_GENES_OF_INTEREST.out.tsv  // channel: [ meta, genes_of_interest.tsv ]
    comparison        = CHECKSTRAIN_COMPARE_BAKTA_OSMOTOOL.out.tsv    // channel: [ meta, bakta_osmotool_comparison.tsv ]
    bakta_tsv         = BAKTA_BAKTA.out.tsv                           // channel: [ meta, bakta.tsv ]
    osmotool_counts   = OSMOTOOL_ANNOTATE.out.counts                  // channel: [ meta, gene_counts.tsv ]
    versions          = ch_versions                                  // channel: [ path(versions.yml) ]
}
