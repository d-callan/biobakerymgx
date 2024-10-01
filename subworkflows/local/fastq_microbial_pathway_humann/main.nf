//
// SUBWORKFLOW: Identify gene families and pathways associated with reads using HUMAnN 3
//

include { HUMANN_DOWNLOADCHOCOPHLANDB      } from '../../../modules/local/humann/downloadchocophlandb/main'
include { HUMANN_DOWNLOADUNIREFDB          } from '../../../modules/local/humann/downloadunirefdb/main'
include { HUMANN_HUMANN                    } from '../../../modules/local/humann/humann/main'
include { HUMANN_JOIN as JOIN_GENES        } from '../../../modules/local/humann/join/main'
include { HUMANN_JOIN as JOIN_PATHABUND    } from '../../../modules/local/humann/join/main'
include { HUMANN_JOIN as JOIN_PATHCOV      } from '../../../modules/local/humann/join/main'
include { HUMANN_JOIN as JOIN_EC           } from '../../../modules/local/humann/join/main'
include { HUMANN_REGROUP                   } from '../../../modules/local/humann/regroup/main'
include { HUMANN_RENAME                    } from '../../../modules/local/humann/rename/main'
include { HUMANN_RENORM                    } from '../../../modules/local/humann/renorm/main'

workflow FASTQ_MICROBIAL_PATHWAY_HUMANN {

    take:
    processed_reads_fastq_gz      // channel: [ val(meta), [ processed_reads_1.fastq.gz, processed_reads_2.fastq.gz ] ] (MANDATORY)
    metaphlan_profile             // channel: [ val(meta2), metaphlan_profile.tsv ] (MANDATORY)
    chocophlan_db                 // channel: [ chocophlan_db ] (OPTIONAL)
    chocophlan_db_version        // value: '' (OPTIONAL)
    uniref_db                     // channel: [ uniref_db ] (OPTIONAL)
    uniref_db_version             // value: '' (OPTIONAL)

    main:

    ch_versions = Channel.empty()

    // if chocophlan_db exists, skip HUMANN_DOWNLOADCHOCOPHLANDB
    if ( chocophlan_db ){
        ch_chocophlan_db = chocophlan_db
    } else {
        //
        // MODULE: Download ChocoPhlAn database
        //
        ch_chocophlan_db = HUMANN_DOWNLOADCHOCOPHLANDB ( chocophlan_db_version ).chocophlan_db
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
    ch_humann_genefamilies_raw = HUMANN_HUMANN ( processed_reads_fastq_gz, metaphlan_profile, ch_chocophlan_db, ch_uniref_db ).genefamilies
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
    ch_humann_ec = HUMANN_REGROUP(ch_humann_genefamilies_cpm, 'uniref90_rxn').regroup
    ch_versions = ch_versions.mix(HUMANN_REGROUP.out.versions)

    //
    // MODULE: rename ec number outputs to include descriptors
    //
    ch_humann_ec_renamed = HUMANN_RENAME(ch_humann_ec, 'ec').renamed // TODO make sure 'ec' is valid arg
    ch_versions = ch_versions.mix(HUMANN_RENAME.out.versions)

    //
    // MODULE: join gene abundances across all samples into one file
    //
    ch_humann_genefamilies_joined = JOIN_GENES(ch_humann_genefamilies_cpm, 'genefamilies').joined

    //
    // MODULE: join ec abundances across all samples into one file
    //
    ch_humann_ec_joined = JOIN_EC(ch_humann_ec_renamed, 'ec').joined // TODO check the file name pattern

    //
    // MODULE: join pathway abundances across all samples into one file
    //
    ch_humann_pathabundance_joined = JOIN_PATHABUND(ch_humann_pathabundance_raw, 'pathabundance').joined

    //
    // MODULE: join pathway coverage across all samples into one file
    //
    ch_humann_pathcoverage_joined = JOIN_PATHCOV(ch_humann_pathcoverage_raw, 'pathcoverage').joined

    emit:
    humann_genefamilies       = ch_humann_genefamilies_joined     // channel: [ val(meta), genefamilies.tsv ]
    humann_ec                 = ch_humann_ec_joined               // channel: [ val(meta), read_counts.tsv ]
    humann_pathabundance      = ch_humann_pathabundance_joined    // channel: [ val(meta), pathabundance.tsv ]
    humann_pathcoverage       = ch_humann_pathcoverage_joined     // channel: [ val(meta), pathcoverage.tsv ]
    versions                  = ch_versions                       // channel: [ versions.yml ]
}
