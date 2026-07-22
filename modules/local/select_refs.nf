process SELECT_REFS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(assembly), val(include_vals), val(exclude_vals)
    path refs
    
    output:
    tuple val(meta), path("subset.fa.gz"),             emit: fa,    optional: true
    tuple val(meta), path("subset.jsonl.gz"),          emit: jsonl, optional: true
    tuple val(meta), path("genome-fraction.csv"), emit: gf, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:

    prefix = task.ext.prefix ?: "${meta.id}"
    tool = "vaper_refs.py"
    """
    ${tool} \\
        --refs "${refs}" \\
        --query "${assembly}" \\
        --min-query-cov ${params.ref_min_query_cov} \\
        --min-ref-cov ${params.ref_min_ref_cov} \\
        ${include_vals ? "--include '" + include_vals.join(',') + "'" : ''} \\
        ${exclude_vals ? "--exclude '" + exclude_vals.join(',') + "'" : ''}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
