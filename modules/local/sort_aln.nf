process SORT_ALN {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::samtools=1.21"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.subsampled"
    """
    # coordinate-sort the subsampled alignment
    samtools sort ${args} -@ ${task.cpus} -o ${prefix}.bam ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
