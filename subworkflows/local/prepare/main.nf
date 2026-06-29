//
// Check input samplesheet and get read channels
//

include { FASTERQDUMP        } from '../../../modules/local/fasterqdump'
include { SEQTK_SAMPLE       } from '../../../modules/nf-core/seqtk/sample/main'
include { SRA_HUMAN_SCRUBBER } from '../../../modules/local/sra_human_scrubber/main'
include { FASTQC             } from '../../../modules/nf-core/fastqc/main'
include { FASTP              } from '../../../modules/local/fastp/main'
include { REF_CSV2JSONL      } from '../../../modules/local/ref_csv2jsonl'
include { VALIDATE_REFS      } from '../../../modules/local/validate_refs'

workflow PREPARE {
    take:
    ch_samplesheet // channel: [ val(meta), ... ]
    ref_set        // file: /path/to/ref_set.jsonl.gz
    
    main:
    ch_versions = Channel.empty()

    /* 
    =============================================================================================================================
        REFERENCES
    =============================================================================================================================
    */
    ch_ref_set = Channel.empty()
    if(ref_set){
        ch_ref_set = Channel.fromList(ref_set)

        // Reference supplied as CSV
        ch_ref_set_csv = ch_ref_set.filter { it.name.endsWith('.csv') }
        REF_CSV2JSONL (
            ch_ref_set_csv,
            ch_ref_set_csv
                .splitCsv(header: true)
                .map { file(it.assembly, checkIfExists: true) }
                .collect()
        )
        ch_versions = ch_versions.mix(REF_CSV2JSONL.out.versions)
            
        // Combine all references into a single file and validate
        VALIDATE_REFS (
            ch_ref_set
                .filter { it.name.endsWith('.jsonl') || it.name.endsWith('.jsonl.gz') }
                .concat(REF_CSV2JSONL.out.jsonl)
                .flatten()
                .collect()
        )
        ch_versions = ch_versions.mix(VALIDATE_REFS.out.versions)
        ch_ref_set = VALIDATE_REFS.out.jsonl
    }
    

    /* 
    =============================================================================================================================
        SAMPLESHEET
    =============================================================================================================================
    */
    ch_ref_man = ch_samplesheet.map{ [ it.meta, it.reference, it.include, it.exclude ] }
    ch_reads   = ch_samplesheet.filter { !it.sra }.map{ [it.meta, it.fastq_12] }

    /* 
    =============================================================================================================================
        PROCESS READS
    =============================================================================================================================
    */

    // MODULE: Download reads from SRA
    FASTERQDUMP(
        ch_samplesheet.filter{ it -> it.sra }.map{ it -> [ it.meta, it.sra ] }
    )
    ch_versions = ch_versions.mix(FASTERQDUMP.out.versions)
    // Update reads channel with SRA reads
    ch_reads
        .concat(FASTERQDUMP.out.reads)
        .set{ ch_reads }
    
    // Downsample high-coverage samples
    if (params.max_reads) {
        ch_reads
            .map { meta, reads ->
                def read_count
                try {
                    read_count = reads[0].countFastq() * 2
                } catch (Exception e) {
                    log.warn "${meta.id}: Failed to count reads - will attempt downsampling anyway!"
                    read_count = params.max_reads + 1
                }
                [meta: meta, reads: reads, n: read_count]
            }
            .branch {
                ok:   it.n <= params.max_reads
                high: it.n > params.max_reads
            }
            .set { ch_reads_branched }

        // Downsample R1 and R2 separately
        SEQTK_SAMPLE(
            ch_reads_branched.high
                .map { [it.meta, it.reads[0], params.max_reads] }
                .concat(
                    ch_reads_branched.high
                        .map { [it.meta, it.reads[1], params.max_reads] }
                )
        )
        ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions.first())

        // Recombine downsampled R1/R2 pairs
        ch_reads = SEQTK_SAMPLE
            .out
            .reads
            .groupTuple(by: 0)
            .concat(ch_reads_branched.ok.map { [it.meta, it.reads] })
    }

    // MODULE: NCBI Human Read Scrubber    
    if(params.scrub_reads){
        SRA_HUMAN_SCRUBBER(
            ch_reads
        )
        ch_versions = ch_versions.mix(SRA_HUMAN_SCRUBBER.out.versions)

        SRA_HUMAN_SCRUBBER.out.reads.set{ ch_reads }
    }

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    //
    // MODULE: Run Fastp
    //

    FASTP (
        ch_reads.map{ meta, reads -> [ meta, reads, [] ]},
        false,
        false,
        false
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    FASTP.out.reads.set{ch_reads}

    emit:
    reads      = ch_reads               // channel: [ val(meta), [ path(fastq_1), path(fastq_2) ], path(fastq_l), val(sra), val/path(references) ]
    ref_man    = ch_ref_man             // channel: [ val(meta), path(reference), val(inclusions), val(exclusions) ]
    ref_set    = ch_ref_set             // channel: path(reference.valid.jsonl.gz)
    fastqc_zip = FASTQC.out.zip         // channel: path(.zip)
    fastp_json = FASTP.out.json         // channel: [ val(meta), path(json) ]
    versions   = ch_versions            // channel: [ versions.yml ]
}



