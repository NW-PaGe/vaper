process SUMMARIZE_TAXA {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(sm_meta)
    
    output:
    tuple val(meta), path("*.taxa-summary.json"), emit: sm_summary
    tuple val(meta), path("*.jpg"), emit: plot,  optional: true
    path "versions.yml",                         emit: versions



    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    tool="vaper_metagenome.py"
    """
    # summarize taxa at >= 1X coverage and >= 1% relative abundance
    ${tool} \\
        "${sm_meta}" \\
        ${prefix}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
