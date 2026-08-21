process CHECKSTRAIN_COMPARE_BAKTA_OSMOTOOL {
    tag "$meta.id"
    label 'process_single'

    container "docker://barbarahelena/osmotool:0.5.0"

    input:
    tuple val(meta), path(genes_of_interest_tsv), path(gene_coordinates_tsv)

    output:
    tuple val(meta), path("*.bakta_osmotool_comparison.tsv"), emit: tsv
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'compare_bakta_osmotool.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bakta_osmotool_comparison.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
