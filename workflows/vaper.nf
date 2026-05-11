/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_vaper_pipeline'

include { PREPARE                } from '../subworkflows/local/prepare/main'
include { CLASSIFY               } from '../subworkflows/local/classify/main'
include { ASSEMBLE               } from '../subworkflows/local/assemble/main'

include { SUMMARYLINE            } from '../modules/local/summaryline'
include { COMBINE_SUMMARYLINES   } from '../modules/local/combine-summary'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VAPER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    /* 
    =============================================================================================================================
        PREPARE INPUT
    =============================================================================================================================
    */
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files

    ref_set = params.ref_set ? file(params.ref_set) : []

    PREPARE (
        ch_samplesheet,
        (ref_set instanceof List) ? ref_set : [ ref_set ]
    )
    ch_versions = ch_versions.mix(PREPARE.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE.out.fastqc_zip.collect{it[1]})

    ch_reads   = PREPARE.out.reads
    ch_ref_man = PREPARE.out.ref_man
    ch_ref_set = PREPARE.out.ref_set

    /* 
    =============================================================================================================================
        CLASSIFY VIRUSES
    =============================================================================================================================
    */ 

    //
    // SUBWORKFLOW: Classify viruses
    //
    CLASSIFY (
        ch_reads,
        ch_ref_set,
        ch_ref_man
    )
    ch_versions = ch_versions.mix(CLASSIFY.out.versions)

    /* 
    =============================================================================================================================
        ASSEMBLE VIRUSES
    =============================================================================================================================
    */
    // SUBWORKFLOW: Create consensus assemblies
    ASSEMBLE (
        CLASSIFY.out.ref.combine(ch_reads, by: 0)
    )
    ch_versions = ch_versions.mix(ASSEMBLE.out.versions)

    /* 
    =============================================================================================================================
        SUMMARIZE RESULTS
    =============================================================================================================================
    */

    ch_reads
        .map{ meta, reads -> [ meta ] }
        .join(PREPARE.out.fastp_json, by: 0, remainder: true)
        .join(CLASSIFY.out.sm_summary, by: 0, remainder: true)
        .join(CLASSIFY.out.ref_meta, by: 0, remainder: true)
        .join(ASSEMBLE.out.bam_stats, by: 0, remainder: true)
        .join(ASSEMBLE.out.assembly_stats, by: 0, remainder: true)
        .map{ meta, fastp, sm, refs_meta, bam_stats, asm_stats -> [ 
            meta,
            fastp ? fastp : [],
            sm ? sm : [],
            refs_meta ? refs_meta : [],
            bam_stats ? bam_stats : [],
            asm_stats ? asm_stats : []
        ]}
        .set{ ch_all }

    // MODULE: Create summaryline for each sample 
    SUMMARYLINE (
       ch_all
    )
    ch_versions = ch_versions.mix(SUMMARYLINE.out.versions)

    // MODULE: Combine summarylines
    SUMMARYLINE
        .out
        .json
        .map{ meta, json -> [ json ] }
        .collect()
        .set{ all_summaries }
    
    COMBINE_SUMMARYLINES (
        all_summaries,
        Channel.value(groovy.json.JsonOutput.toJson(params)).collectFile(name: 'params.json')
    )
    ch_versions = ch_versions.mix(COMBINE_SUMMARYLINES.out.versions)

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'vaper_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
