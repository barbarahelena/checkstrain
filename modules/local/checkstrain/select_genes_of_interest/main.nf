process CHECKSTRAIN_SELECT_GENES_OF_INTEREST {
    tag "$meta.id"
    label 'process_single'

    container "docker://barbarahelena/osmotool:0.5.0"

    input:
    tuple val(meta), path(bakta_tsv), path(gene_counts_tsv)

    output:
    tuple val(meta), path("*.genes_of_interest.tsv"), emit: tsv
    path "versions.yml",                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'select_genes_of_interest.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.genes_of_interest.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
