//
// Check input samplesheet and get read channels
//

include { BWA_MEM        } from '../../../modules/local/bwa_mem'
include { RASUSA_ALN     } from '../../../modules/local/rasusa_aln'
include { SORT_ALN       } from '../../../modules/local/sort_aln'
include { IVAR_CONSENSUS } from '../../../modules/local/ivar_consensus'
include { ASSEMBLY_TIDY  } from '../../../modules/local/assembly_tidy'
include { CONDENSE       } from '../../../modules/local/condense'
include { BAM_STATS      } from '../../../modules/local/bam_stats'
include { ASSEMBLY_STATS } from '../../../modules/local/assembly_stats'
include { EXPORT_ALN     } from '../../../modules/local/export_aln'
include { EXPORT_FASTQ   } from '../../../modules/local/export_fastq'


workflow ASSEMBLE {
    take:
    ch_ref_set // channel: [meta, ref_path, reads]

    main:

    ch_versions = Channel.empty()
    ch_bam      = Channel.empty()

    /* 
    =============================================================================================================================
        MAP READS
    =============================================================================================================================
    */
    // MODULE: Run BWA MEM
    BWA_MEM (
        ch_ref_set
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    BWA_MEM.out.bam.set{ ch_bam_full }

    // Module: Get BAM stats (always from the full alignment)
    BAM_STATS (
        ch_bam_full
    )
    ch_versions = ch_versions.mix(BAM_STATS.out.versions)

    /* 
    =============================================================================================================================
        SUBSAMPLE ALIGNMENTS
        - randomly subsample each alignment to a maximum per-position depth (--cons_max_depth)
        - the full alignment is used for stats / read export; the subsampled alignment is used for consensus calling
    =============================================================================================================================
    */
    // MODULE: Randomly subsample alignments with Rasusa
    RASUSA_ALN (
        ch_bam_full
    )
    ch_versions = ch_versions.mix(RASUSA_ALN.out.versions)

    // MODULE: Coordinate-sort the subsampled alignments
    SORT_ALN (
        RASUSA_ALN.out.bam
    )
    ch_versions = ch_versions.mix(SORT_ALN.out.versions)
    SORT_ALN.out.bam.set{ ch_bam_sub }

    // Assemblies are created separately for each reference contig
    BAM_STATS
        .out
        .coverage
        .splitCsv(header: true, sep: "\t")
        .map{meta, data -> [meta, data["#rname"]]}
        .set{ ch_contigs }

    // channel: [ meta, ref_id, subsampled_bam ]
    ch_contigs
        .combine(ch_bam_sub, by: 0)
        .set{ ch_bam }

    // MODULE: Run Ivar on the subsampled alignment
    IVAR_CONSENSUS (
        ch_bam
    )
    ch_versions = ch_versions.mix(IVAR_CONSENSUS.out.versions)
    IVAR_CONSENSUS.out.consensus.set{ ch_consensus }


    /* 
    =============================================================================================================================
        ASSEMBLY TIDY
        - generate a tidy version of each assembly
        - can include replacing all mixed sites with N or pruning all terminal N
    =============================================================================================================================
    */
    if (params.cons_no_mixed_sites || params.cons_prune_termini){
        ASSEMBLY_TIDY (
            ch_consensus
        )
        ch_versions = ch_versions.mix(ASSEMBLY_TIDY.out.versions)
    }

    /* 
    =============================================================================================================================
        CONDENSE DUPLICATE ASSEMBLIES
    =============================================================================================================================
    */
    CONDENSE(
        ch_consensus
            .groupTuple(by: 0)
            .join(BAM_STATS.out.coverage, by: 0)
            .map{ meta, ref_names, assemblies, stats -> [ meta, assemblies, stats ] }
    )
    
    ch_versions = ch_versions.mix(CONDENSE.out.versions)

    CONDENSE
        .out
        .assembly
        .transpose()
        .map{ meta, assembly -> [ meta, Utils.getRefName(assembly, meta.id), assembly ] }
        .set{ ch_consensus }
    
    ch_consensus.map{ meta, ref_name, assembly -> [meta, ref_name] }.set{ ch_remaining }

    // Publish alignments (full + subsampled) of condensed assemblies
    // channel: [ meta, ref_id, ref_path, full_bam, subsampled_bam ]
    EXPORT_ALN (
        ch_remaining
          .combine( ch_ref_set.map{ meta, ref_path, reads -> [ meta, ref_path ] }, by: 0 )
          .join( ch_contigs.combine( ch_bam_full.join( ch_bam_sub, by: 0 ), by: 0 ), by: [ 0, 1 ] )
    )

    // Export mapped reads for condensed assemblies
    EXPORT_FASTQ (
        BAM_STATS.out.read_list.join(ch_ref_set.map{ meta, ref_path, reads -> [ meta, reads ] }, by: 0)
    )
    ch_versions = ch_versions.mix(EXPORT_FASTQ.out.versions)


    /* 
    =============================================================================================================================
        ASSEMBLY QC METRICS
    =============================================================================================================================
    */
    // MODULE: Gather assembly stats
    ASSEMBLY_STATS (
        ch_consensus
            .combine(ch_ref_set.map{ meta, ref, reads -> [ meta, ref ] }, by: 0)
    )
    ch_versions = ch_versions.mix(ASSEMBLY_STATS.out.versions)
    
    ASSEMBLY_STATS
        .out
        .json
        .groupTuple(by: 0)
        .map{meta, ref_name, json -> [meta, json]}
        .set{ch_assembly_stats}


    emit:
    bam_stats      = BAM_STATS.out.stats // channel: [ val(meta), path(stats) ]
    assembly_stats = ch_assembly_stats    // channel: [ val(meta), path(json) ]
    consensus      = ch_consensus         // channel: [ val(meta), val(ref_name), path(assembly) ]
    versions       = ch_versions
}
