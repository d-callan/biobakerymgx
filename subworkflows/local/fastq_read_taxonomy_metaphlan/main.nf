//
// SUBWORKFLOW: Taxonomically classify reads using MetaPhlAn
//

include { METAPHLAN_MAKEDB                  } from '../../../modules/nf-core/metaphlan/makedb/main'
include { METAPHLAN_METAPHLAN               } from '../../../modules/nf-core/metaphlan/metaphlan/main'
include { METAPHLAN_SGBTOGTDBPROFILE        } from '../../../modules/local/metaphlan/sgbtogtdbprofile/main'
include { METAPHLAN_MERGEMETAPHLANTABLES    } from '../../../modules/nf-core/metaphlan/mergemetaphlantables/main'

workflow FASTQ_READ_TAXONOMY_METAPHLAN {

    take:
    preprocessed_reads_fastq_gz // channel: [ val(meta), [ reads_1.fastq.gz, reads_2.fastq.gz ] ] (MANDATORY)
    metaphlan_sgb2gtbd_file     // channel: [ metaphlan_sgb2gtbd_file.tsv ] (MANDATORY)
    metaphlan_db                // channel: [ metaphlan_db ] (OPTIONAL)
    metaphlan_db_version        // value: mpa_vJan21_TOY_CHOCOPhlAnSGB_202103 (OPTIONAL)

    main:

    ch_versions = Channel.empty()

    // if metaphlan_db exists, skip METAPHLAN_MAKEDB
    if ( metaphlan_db ){
        ch_metaphlan_db = metaphlan_db
    } else {
        //
        // MODULE: Download KneadData database
        //
        ch_metaphlan_db = METAPHLAN_MAKEDB ( metaphlan_db_version ).db
        ch_versions = ch_versions.mix(METAPHLAN_MAKEDB.out.versions)
    }

    //
    // MODULE: Trim and remove human reads
    //
    ch_metaphlan_profile_txt = METAPHLAN_METAPHLAN ( preprocessed_reads_fastq_gz, ch_metaphlan_db ).profile
    ch_versions = ch_versions.mix(METAPHLAN_METAPHLAN.out.versions)

    //
    // MODULE: Convert SGB to GTDB taxonmy
    //
    ch_gtdb_profile_tsv = METAPHLAN_SGBTOGTDBPROFILE ( ch_metaphlan_profile_txt, metaphlan_sgb2gtbd_file ).gtdb_profile
    ch_versions = ch_versions.mix(METAPHLAN_SGBTOGTDBPROFILE.out.versions)

    // combine metaphlan profiles for merging
    ch_metaphlan_profiles_combined_tsv = ch_metaphlan_profile_txt
    .map { [ [ id:'all_samples' ], it[1] ] }
    .groupTuple( sort: 'deep' )

    //
    // MODULE: Merge Metaphlan profiles
    //
    ch_metaphlan_profiles_merged_tsv = METAPHLAN_MERGEMETAPHLANTABLES ( ch_metaphlan_profiles_combined_tsv ).txt
    ch_versions = ch_versions.mix(METAPHLAN_MERGEMETAPHLANTABLES.out.versions)

    emit:
    metaphlan_profiles_merged_tsv = ch_metaphlan_profiles_merged_tsv    // channel: [ val(meta), [ merged_metaphlan_profiles.tsv ] ]
    versions                        = ch_versions                       // channel: [ versions.yml ]
}
