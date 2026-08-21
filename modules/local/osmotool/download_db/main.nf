process OSMOTOOL_DOWNLOAD_DB {
    label 'process_low'

    container "docker://barbarahelena/osmotool:0.5.0"

    input:
    val release

    output:
    path "refdb/*", emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p refdb
    osmotool download-db \\
        --release ${release} \\
        --location refdb \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        osmotool: \$(osmotool --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p refdb/${release}
    touch refdb/${release}/osmo_refdb.dmnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        osmotool: stub
    END_VERSIONS
    """
}
