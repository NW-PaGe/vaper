process EXPORT_ALN {
    tag "${prefix}"
    label 'process_high'
    stageInMode 'copy'

    input:
    tuple val(meta), val(ref_id), path(ref), path(bam), path(bam_sub)
    
    output:
    tuple val(meta), path("*.bam*", includeInputs: true)
    tuple val(meta), path("*.fa*", includeInputs: true)
    path "versions.yml", emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = "${meta.id}-${ref_id}"
    """
    # index the full and subsampled bam
    samtools index ${bam}
    samtools index ${bam_sub}

    gzip -d "${ref}"
    samtools faidx *.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
