process IVAR_CONSENSUS {
    tag "${prefix}"
    label 'process_high'

    conda "bioconda::ivar"
    container "public.ecr.aws/o8h2f0o1/ivar:1.4.2"

    input:
    tuple val(meta), val(ref_id), path(bam)

    output:
    tuple val(meta), val(ref_id), path('*.fa.gz'), emit: consensus
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = "${meta.id}_${ref_id}"

    """
    # setup for pipe
    set -euxo pipefail

    # index the bam
    samtools index ${bam}

    # create mpilup and call consensus
    # depth is capped upstream by random subsampling (RASUSA_ALN; --cons_max_depth),
    # so the mpileup depth limit is disabled here (-d 0)
    samtools mpileup \\
        -aa -A -Q 0 \
        -d 0 \\
        -r ${ref_id} ${bam} | \\
        ivar consensus \\
        -p ${prefix} \\
        -m ${params.cons_allele_depth} \\
        -t ${params.cons_allele_ratio} \\
        -q ${params.cons_allele_qual} \\
        ${args}    
    sed -i 's/>.*/>${prefix}/g' ${prefix}.fa
    gzip ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools 2>&1 | grep "Version" | cut -f 2 -d ' ')
        ivar: \$(ivar version | head -n 1)
    END_VERSIONS
    """
}
