process RASUSA_ALN {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::rasusa=5.1.0"
    container "docker.io/staphb/rasusa:5.1.0"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.unsorted.bam"), emit: bam
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.subsampled"
    """
    # Randomly subsample alignments so that no position exceeds the maximum depth.
    # Input must be coordinate-sorted; output is NOT guaranteed to be sorted.
    rasusa aln \\
        --coverage ${params.cons_max_depth} \\
        --threads ${task.cpus} \\
        ${args} \\
        --output ${prefix}.unsorted.bam \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rasusa: \$(rasusa --version 2>&1 | sed 's/^rasusa //')
    END_VERSIONS
    """
}
