//
// SUBWORKFLOW: Trim, remove human reads, and calculate read counts with KneadData
//

include { KNEADDATA_DATABASE            } from '../../../modules/local/kneaddata/database/main'
include { KNEADDATA_KNEADDATA           } from '../../../modules/local/kneaddata/kneaddata/main'
include { KNEADDATA_READCOUNTS          } from '../../../modules/local/kneaddata/readcounts/main'

workflow FASTQ_READ_PREPROCESSING_KNEADDATA {

    take:
    raw_reads_fastq_gz      // channel: [ val(meta), [ reads_1.fastq.gz, reads_2.fastq.gz ] ] (MANDATORY)
    kneaddata_db            // channel: [ kneaddata_db ] (OPTIONAL)
    kneaddata_db_version    // value: 'human_genome' (OPTIONAL)

    main:

    ch_versions = Channel.empty()

    // if kneaddata_db exists, skip KNEADDATA_DATABASE
    if ( kneaddata_db ){
        ch_kneaddata_db = kneaddata_db
    } else {
        //
        // MODULE: Download KneadData database
        //
        ch_kneaddata_db = KNEADDATA_DATABASE ( kneaddata_db_version ).kneaddata_db
        ch_versions = ch_versions.mix(KNEADDATA_DATABASE.out.versions)
    }

    //
    // MODULE: Trim and remove human reads
    //
    ch_preprocessed_reads_fastq_gz = KNEADDATA_KNEADDATA ( raw_reads_fastq_gz, ch_kneaddata_db ).preprocessed_reads
    ch_kneaddata_logs = KNEADDATA_KNEADDATA.out.kneaddata_log
    ch_versions = ch_versions.mix(KNEADDATA_KNEADDATA.out.versions)

    // collect log files and store in a directory
    ch_combined_kneaddata_logs = ch_kneaddata_logs
        .map { [ [ id:'all_samples' ], it[1] ] }
        .groupTuple( sort: 'deep' )

    //
    // MODULE: Calculate read counts during preprocessing
    //
    kneaddata_read_counts_tsv = KNEADDATA_READCOUNTS ( ch_combined_kneaddata_logs ).kneaddata_read_counts
    ch_versions = ch_versions.mix(KNEADDATA_READCOUNTS.out.versions)

    emit:
    preprocessed_fastq_gz   = ch_preprocessed_reads_fastq_gz    // channel: [ val(meta), [ reads_1.fastq.gz, reads_2.fastq.gz  ] ]
    read_counts_tsv         = kneaddata_read_counts_tsv         // channel: [ val(meta), read_counts.tsv ]
    versions                = ch_versions                       // channel: [ versions.yml ]
}
