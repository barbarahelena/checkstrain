process OSMOTOOL_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    container "docker://barbarahelena/osmotool:0.5.0"

    input:
    tuple val(meta), path(fasta)
    path osmo_refdb

    output:
    tuple val(meta), path("*.gene_counts.tsv"),      emit: counts
    tuple val(meta), path("*.aln_stats.tsv"),        emit: stats
    tuple val(meta), path("*.systems.tsv"),          emit: systems
    tuple val(meta), path("*.gene_coordinates.tsv"), emit: coordinates
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    osmotool annotate \\
        ${osmo_refdb} \\
        ${fasta} \\
        --out_prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        osmotool: \$(osmotool --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gene_counts.tsv
    touch ${prefix}.aln_stats.tsv
    touch ${prefix}.systems.tsv
    touch ${prefix}.gene_coordinates.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        osmotool: stub
    END_VERSIONS
    """
}
