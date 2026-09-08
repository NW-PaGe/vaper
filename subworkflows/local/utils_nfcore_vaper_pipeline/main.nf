//
// Subworkflow with functionality specific to the doh-jdj0303/vaper pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        "",
        "",
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { create_sample_channel(it) }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Custom samplesheet process
//

// Function to prepare the samplesheet
def create_sample_channel(row) {
    def (
        sample, fastq_1, fastq_2, sra, ref_file, 
        ref_taxon, ref_species, ref_segment, ref_name, ref_include, ref_exclude
    ) = row

    // Check that at least one read type was provided
    if( ! (fastq_1 && fastq_2) && ! sra ){ exit 1, "ERROR: Reads must be provided via the 'fastq_1' and 'fastq_2', or 'sra' columns for ${sample}." }
    // Check that only one read type was provided
    if( ((fastq_1 && fastq_2) && sra) ){ exit 1, "ERROR: Multiple read inputs provided - either fastq_1 & fastq_2 or sra" }

    // ===================================================
    //    References
    // ===================================================
    // ---------- reference paths (not in reference set) ----------
    def ref_file_ss  = ref_file        ? ref_file.toString().tokenize(';').collect { file(it, checkIfExists: true) }        : []
    def ref_file_cli = params.ref_file ? params.ref_file.toString().tokenize(' ').collect { file(it, checkIfExists: true) } : []
    def ref_files    = ((ref_file_ss + ref_file_cli) as Set) as List   // dedupe

    // ---------- reference names (in reference set) ----------
    def ref_names_ss  = ref_name        ? ref_name.tokenize(';')        : []
    def ref_names_cli = params.ref_name ? params.ref_name.tokenize(';') : []
    def ref_names_pattern = ((ref_names_ss + ref_names_cli) as Set) as List
    ref_names_pattern = ref_names_pattern ? ref_names_pattern.collect { "name=${it}" } : []

    // ---------- segment ----------
    def ref_segment_ss  = ref_species        ? ref_species.tokenize(';')        : []
    def ref_segment_cli = params.ref_species ? params.ref_species.tokenize(';') : []
    def ref_segment_pattern = ((ref_segment_ss + ref_segment_cli) as Set) as List
    ref_segment_pattern = ref_segment_pattern ? ref_segment_pattern.collect { "species=${it}" } : []

    // ---------- species ----------
    def ref_species_ss  = ref_species        ? ref_species.tokenize(';')        : []
    def ref_species_cli = params.ref_species ? params.ref_species.tokenize(';') : []
    def ref_species_pattern = ((ref_species_ss + ref_species_cli) as Set) as List
    ref_species_pattern = ref_species_pattern ? ref_species_pattern.collect { "species=${it}" } : []

    // ---------- taxon ----------
    def ref_taxon_ss  = ref_taxon        ? ref_taxon.tokenize(';')        : []
    def ref_taxon_cli = params.ref_taxon ? params.ref_taxon.tokenize(';') : []
    def ref_taxon_pattern = ((ref_taxon_ss + ref_taxon_cli) as Set) as List
    ref_taxon_pattern  = ref_taxon_pattern ? ref_taxon_pattern.collect { "taxon=${it}" } : []

    // ---------- include/exclude ----------
    // keep your logic: include is the base, then add ref_* and dedupe via Set
    def include_ss  = ref_include ? ref_include.tokenize(';') : []
    def include_cli = ref_include ? ref_include.tokenize(';') : []
    def include_set = ((include_ss + include_cli + ref_names_pattern + ref_segment_pattern + ref_species_pattern + ref_taxon_pattern) as Set) as List

    // exclude stays as you had it (semicolon-sep → comma string)
    def exclude_ss  = ref_exclude        ? ref_exclude.tokenize(';')        : []
    def exclude_cli = params.ref_exclude ? params.ref_exclude.tokenize(';') : []
    def exclude_set = ((exclude_ss + exclude_cli) as Set) as List
    exclude_set     = exclude_set ? ((ref_taxon_pattern || ref_species_pattern || ref_segment_pattern || ref_names_pattern) ? ['*'] : []) : []
    
    //// Build output
    def out = [ 
        meta: sample + ['single_end': false], 
        fastq_12: [ fastq_1, fastq_2 ], 
        sra: sra, 
        reference: ref_files,
        include: include_set,
        exclude: exclude_set
    ]
    
    return out
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
