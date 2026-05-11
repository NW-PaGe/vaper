process FINALIZE_REFS {
    tag "${meta.id}"
    label 'process_low'
    stageInMode 'copy'

    input:
    tuple val(meta), path(ref, stageAs: "input/*")

    output:
    tuple val(meta), path("refs.fa.gz"), emit: refs

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gzip input/* || true
    cat input/* > refs.fa.gz
    """
}
