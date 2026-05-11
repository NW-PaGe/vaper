//
// Check input samplesheet and get read channels
//
include { SHOVILL             } from '../../../modules/nf-core/shovill/main'
include { SOURMASH_SKETCH     } from '../../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_GATHER     } from '../../../modules/nf-core/sourmash/gather/main'
include { SOURMASH_METAGENOME } from '../../../modules/local/sourmash_metagenome'
include { KRONA_KTIMPORTTEXT  } from '../../../modules/nf-core/krona/ktimporttext/main'
include { SELECT_REFS         } from '../../../modules/local/select_refs'
include { SUMMARIZE_TAXA      } from '../../../modules/local/summarize_taxa'
include { KITCHEN_SINK        } from '../../../modules/local/kitchen_sink'
include { FINALIZE_REFS       } from '../../../modules/local/finalize_refs'

workflow CLASSIFY {
    take:
    ch_reads    // channel: [ val(meta), path(reads) ]
    ch_ref_set  // channel: [ path(references.jsonl.gz) ]
    ch_ref_man  // channel: [ val(meta), path(ref), val(include), val(exclude)  ]

    main:

    ch_versions  = Channel.empty()
    ch_ref_meta  = Channel.empty()
    ch_ref       = Channel.empty()

    /* 
    =============================================================================================================================
        CREATE SOURMASH SKETCHES
    =============================================================================================================================
    */
    // MODULE: Sample sketches
    SOURMASH_SKETCH (
        ch_reads.map{ meta, reads -> [ meta, reads.get(0) ] } // only uses forward read
    )
    ch_versions = ch_versions.mix(SOURMASH_SKETCH.out.versions)

    /* 
    =============================================================================================================================
        DETERMINE DOMINANT VIRAL TAXA
    =============================================================================================================================
    */
    ch_sm_summary = Channel.empty()

    if (params.metagenome || params.ref_mode == 'kitchen-sink'){
        // MODULE: Classify viral species in each sample using Sourmash
        SOURMASH_GATHER (
            SOURMASH_SKETCH.out.signatures,
            params.sm_db,
            false,
            false,
            false,
            false
        )
        ch_versions = ch_versions.mix(SOURMASH_GATHER.out.versions)

        // MODULE: Summarize taxa from sourmash gather
        SOURMASH_METAGENOME (
            SOURMASH_GATHER.out.result,
            file(params.sm_taxa, checkIfExists: true)
        )
        ch_versions = ch_versions.mix(SOURMASH_METAGENOME.out.versions)

        // MODULE: Create Krona plot from Sourmash metagenome output
        KRONA_KTIMPORTTEXT (
            SOURMASH_METAGENOME.out.krona
        )
        ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions)

        // Combine reference mapping, sample taxonomy, and reference taxonomy
        SUMMARIZE_TAXA(
            SOURMASH_METAGENOME.out.result
        )
        ch_versions = ch_versions.mix(SUMMARIZE_TAXA.out.versions)
        // Acount for samples that had no metagenome
        SUMMARIZE_TAXA
            .out
            .sm_summary
            .join(ch_reads.map{ meta, reads -> meta }, by: 0, remainder: true)
            .map{ meta, summary -> if(summary){
                    [ meta, summary ]
                }else{
                    fwork = file(workflow.workDir).resolve("${meta.id}-sm_summary.txt")
                    fwork.text = '100% unclassified'
                    [ meta, fwork ]
                }
            }.set{ ch_sm_summary }
    }

    /* 
    =============================================================================================================================
        REFERENCE SELECTION: Accurate Mode
    =============================================================================================================================
    */
    if (params.ref_mode in ["standard", "sensitive"] ){
        // MODULE: Run Shovill
        SHOVILL (
            ch_reads
        )
        ch_versions = ch_versions.mix(SHOVILL.out.versions)

        // MODULE: Accurate reference selection
        SELECT_REFS(
            SHOVILL
                .out
                .contigs
                .join( ch_ref_man.map{meta, ref, inc, exc -> [meta, inc, exc]} ),
            ch_ref_set
        )

        SELECT_REFS.out.fa.map{meta, ref -> [meta, ref]}.set{ ch_ref }
        SELECT_REFS.out.jsonl.set{ ch_ref_meta }
        
    }

    /* 
    =============================================================================================================================
        REFERENCE SELECTION: Kitchen Sink Mode
        - Downloads NCBI assembly for all Sourmash hits and adds them to the reference list
    =============================================================================================================================
    */
    if (params.ref_mode == "kitchen-sink"){
        KITCHEN_SINK (
            SOURMASH_GATHER.out.result.map{meta, file -> file}.collect()
        )
        ch_versions = ch_versions.mix(SOURMASH_GATHER.out.versions)

        SOURMASH_GATHER
            .out
            .result
            .map{ meta, result -> [meta]}
            .combine(KITCHEN_SINK.out.jsonl)
            .set{ch_ref_meta}

        KITCHEN_SINK
            .out
            .fa
            .flatMap{ ref -> 
                ref.collect{ 
                    it -> [ getRefName(assembly, null), it ] 
                }  
            }
            .set{ ch_ref_map }

        KITCHEN_SINK
            .out
            .csv
            .splitCsv(header: true)
            .map{[it.ref_id, it.sample_id]}
            .combine(ch_ref_map, by: 0)
            .map{ ref_name, sample_id, ref_path -> [ sample_id, ref_name, ref_path ]}
            .set{ ch_ref_map }
        
        ch_reads
            .map{ meta, reads -> [ meta.id, meta ] }
            .combine(ch_ref_map, by: 0)
            .map{ sample_id, meta, ref_name, ref_path -> [ meta, ref_path ] }
            .set{ ch_ref }
    }

    /*
    =============================================================================================================================
        REFERENCE SELECTION: Final Clean Up
        - Add selected reference paths
        - Add supplied reference paths
    =============================================================================================================================
    */
    // Combine any manually supplied reference paths
    ch_ref
        .concat(
            ch_ref_man
                .map{meta, ref, inc, exc -> [meta, ref]}
                .transpose()
                .filter{meta, ref -> ref}
        )
        .groupTuple(by: 0)
        .set{ ch_ref }

    FINALIZE_REFS (
        ch_ref
    )

    ch_ref = FINALIZE_REFS.out.refs

    emit:
    ref        = ch_ref                 // channel: [ val(sample_meta), path(ref.fa.gz) ]
    ref_meta   = ch_ref_meta            // channel: [ path(jsonl.gz) ]
    sm_summary = ch_sm_summary          // channel: [ val(meta), val(result) ]
    versions   = ch_versions            // channel: [ versions.yml ]
}

// Remove sample prefix and FASTA extensions, preserving periods in names.
// sampleName: the ID of the sample (string), used to strip "sample-"
def getRefName(filename, String sampleName) {
    // Convert to string and strip any path
    def name = filename.getName()

    // Strip sampleName prefix: "<sampleName>-"
    if (sampleName && name.startsWith(sampleName + "_")) {
        name = name.substring((sampleName + "_").length())
    }

    // Handle .gz first
    if (name.endsWith(".gz")) {
        // strip .gz
        name = name[0..-4]
    }

    // Now strip FASTA-style extensions
    def fastaExts = [".fa", ".fna", ".fasta"]
    for (ext in fastaExts) {
        if (name.toLowerCase().endsWith(ext)) {
            name = name[0..-(ext.size() + 1)]
            break
        }
    }

    return name
}


