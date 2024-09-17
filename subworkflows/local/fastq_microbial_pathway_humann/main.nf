//
// SUBWORKFLOW: Identify gene families and pathways associated with reads using HUMAnN 3
//

include { HUMANN_DOWNLOADCHOCOPHLANDB      } from '../../../modules/local/humann/downloadchocophlandb/main'
include { HUMANN_DOWNLOADUNIREFDB          } from '../../../modules/local/humann/downloadunirefdb/main'
include { HUMANN_HUMANN                    } from '../../../modules/local/humann/humann/main'
include { HUMANN_JOIN                      } from '../../../modules/local/humann/join/main'
include { HUMANN_REGROUP                   } from '../../../modules/local/humann/regroup/main'
include { HUMANN_RENAME                    } from '../../../modules/local/humann/rename/main'
include { HUMANN_RENORM                    } from '../../../modules/local/humann/renorm/main'

workflow FASTQ_MICROBIAL_PATHWAY_HUMANN {

    take:
    processed_reads_fastq_gz      // channel: [ val(meta), [ processed_reads_1.fastq.gz, processed_reads_2.fastq.gz ] ] (MANDATORY)
    metaphlan_profile             // channel: [ val(meta2), metaphlan_profile.tsv ] (MANDATORY)
    chocophlan_db                 // channel: [ chocophlan_db ] (OPTIONAL)
    chochophlan_db_version        // value: '' (OPTIONAL)
    uniref_db                     // channel: [ uniref_db ] (OPTIONAL)
    uniref_db_version             // value: '' (OPTIONAL)

    main:

    ch_versions = Channel.empty()

    // if chocophlan_db exists, skip HUMANN_DOWNLOADCHOCOPHLANDB
    if ( chocophlan_db ){
        ch_chocophlan_db = chochophlan_db
    } else {
        //
        // MODULE: Download ChocoPhlAn database
        //
        ch_chocophlan_db = HUMANN_DOWNLOADCHOCOPHLANDB ( chochophlan_db_version ).chochophlan_db
        ch_versions = ch_versions.mix(HUMANN_DOWNLOADCHOCOPHLANDB.out.versions)
    }

    // if uniref_db exists, skip HUMANN_DOWNLOADUNIREFDB
    if ( uniref_db ){
        ch_uniref_db = uniref_db
    } else {
        //
        // MODULE: Download UniRef database
        //
        ch_uniref_db = HUMANN_DOWNLOADUNIREFDB ( uniref_db_version ).uniref_db
        ch_versions = ch_versions.mix(HUMANN_DOWNLOADUNIREFDB.out.versions)
    }

    //
    // MODULE: Run HUMAnN 3 for raw outputs
    //
    ch_humann_genefamilies_raw = HUMANN_HUMANN ( processed_reads_fastq_gz, metaphlan_profile, ch_chochophlan_db, ch_uniref_db ).genefamilies
    ch_humann_pathabundance_raw = HUMANN_HUMANN.out.pathabundance
    ch_humann_pathcoverage_raw = HUMANN_HUMANN.out.pathcoverage // TODO is this still right? looking at humann docs, might not get this file any longer?
    ch_humann_logs = HUMANN_HUMANN.out.log
    ch_versions = ch_versions.mix(HUMANN_HUMANN.out.versions)

    // collect log files and store in a directory
    ch_combined_humann_logs = ch_humann_logs
        .map { [ [ id:'all_samples' ], it[1] ] }
        .groupTuple( sort: 'deep' )

    //
    // MODULE: renormalize raw gene families from HUMAnN outputs to cpm
    //
    ch_humann_genefamilies_cpm = HUMANN_RENORM ( ch_humann_genefamilies_raw, 'cpm' ).renorm
    ch_versions = ch_versions.mix(HUMANN_RENORM.out.versions)

    //
    // MODULE: regroup cpm gene families to EC numbers
    //
    ch_humann_ec = HUMANN_REGROUP( ch_humann_genefamilies_cpm, 'ec').regroup // TODO make sure 'ec' is still valid arg
    ch_versions = ch_versions.mix(HUMANN_REGROUP.out.versions)

    //
    // MODULE: rename ec number outputs to include descriptors
    //
    ch_humann_ec_renamed = HUMANN_RENAME (ch_humann_ec, 'ec').rename // TODO make sure 'ec' is valid arg
    ch_versions = ch_versions.mix(HUMANN_RENAME.out.versions)

    // TODO join all outputs as necessary, then update emit below
    // TODO need to modify modules to return output dirs i suppose first, so they can be passed to join module

    emit:
    humann_genefamilies_cpm   = ch_humann_genefamilies_cpm        // channel: [ val(meta), [ reads_1.fastq.gz, reads_2.fastq.gz  ] ]
    humann_ec                 = ch_humann_ec_renamed              // channel: [ val(meta), read_counts.tsv ]
    humann_pathabundance      = ch_humann_pathabundance_raw       // channel: [ val(meta), pathabundance.tsv ]
    humann_pathcoverage       = ch_humann_pathcoverage_raw        // channel: [ val(meta), pathcoverage.tsv ]
    versions                  = ch_versions                       // channel: [ versions.yml ]
}
