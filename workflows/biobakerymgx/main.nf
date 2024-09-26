/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOWS: Local subworkflows
//
include { FASTQ_READ_PREPROCESSING_KNEADDATA    } from '../../subworkflows/local/fastq_read_preprocessing_kneaddata/main'
include { FASTQ_READ_TAXONOMY_METAPHLAN         } from '../../subworkflows/local/fastq_read_taxonomy_metaphlan/main'
include { FASTQ_MICROBIAL_PATHWAY_HUMANN        } from '../../subworkflows/local/fastq_microbial_pathway_humann/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// PLUGIN
//
include { paramsSummaryMap       } from 'plugin/nf-validation'

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                 } from '../../modules/nf-core/fastqc/main'
include { CAT_FASTQ              } from '../../modules/nf-core/cat/fastq/main'
include { MULTIQC                } from '../../modules/nf-core/multiqc/main'

//
// SUBWORKFLOWS: Installed directory from nf-core/subworkflows
//
include { paramsSummaryMultiqc   } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../subworkflows/local/utils_nfcore_biobakerymgx_pipeline'




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BIOBAKERYMGX {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    /*-----------------------------------------------------------------------------------
        Prepare reads for biobakerymgx
    -----------------------------------------------------------------------------------*/
    // group reads from the same sample
    if ( params.run_runmerging ) {
        ch_reads_for_cat_branch = ch_samplesheet
            .map {
                meta, reads ->
                    def meta_new = meta - meta.subMap('replicate')
                    [ meta_new, reads ]
            }
            .groupTuple()
            .map {
                meta, reads ->
                    [ meta, reads.flatten() ]
            }
            .branch {
                meta, reads ->
                // we can't concatenate files if there is not a second run, we branch
                // here to separate them out, and mix back in after for efficiency
                cat: reads.size() > 2
                skip: true
            }

        //
        // MODULE: Concatenate fastq files within samples
        //
        ch_reads_runmerged = CAT_FASTQ ( ch_reads_for_cat_branch.cat ).reads
            .mix( ch_reads_for_cat_branch.skip )
            .map {
                meta, reads ->
                [ meta, [ reads ].flatten() ]
            }
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)
    } else {
        ch_reads_runmerged = ch_samplesheet
    }

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())


    /*-----------------------------------------------------------------------------------
        Read-preprocessing: KneadData
    -----------------------------------------------------------------------------------*/
    if ( params.run_kneaddata ) {
        // create channel from params.kneaddata_db
        if ( !params.kneaddata_db ){
            ch_kneaddata_db = null
        } else {
            ch_kneaddata_db = Channel.value( file( params.kneaddata_db, checkIfExists:true ) )
        }

        //
        // SUBWORKFLOW: KneadData
        //
        ch_preprocessed_fastq_gz = FASTQ_READ_PREPROCESSING_KNEADDATA ( ch_reads_runmerged, ch_kneaddata_db, params.kneaddata_db_version ).preprocessed_fastq_gz
        ch_preprocessed_read_counts_tsv = FASTQ_READ_PREPROCESSING_KNEADDATA.out.read_counts_tsv
        ch_versions = ch_versions.mix(FASTQ_READ_PREPROCESSING_KNEADDATA.out.versions)
    } else {
        ch_preprocessed_fastq_gz = ch_reads_runmerged
        ch_preprocessed_read_counts_tsv = Channel.empty()
    }


    /*-----------------------------------------------------------------------------------
        Taxonomic classification: MetaPhlAn
    -----------------------------------------------------------------------------------*/
    if ( params.run_metaphlan ) {
        // create channel from params.metaphlan_db
        if ( !params.metaphlan_db ){
            ch_metaphlan_db = null
        } else {
            ch_metaphlan_db = Channel.value( file( params.metaphlan_db, checkIfExists:true ) )
        }
        // create channel from params.metaphlan_sgb2gtbd_file
        if ( !params.metaphlan_sgb2gtbd_file ){
            ch_metaphlan_sgb2gtbd_file = null
        } else {
            ch_metaphlan_sgb2gtbd_file = Channel.value( file( params.metaphlan_sgb2gtbd_file, checkIfExists:true ) )
        }

        //
        // SUBWORKFLOW: MetaPhlAn
        //
        ch_read_taxonomy_tsv = FASTQ_READ_TAXONOMY_METAPHLAN (
            ch_preprocessed_fastq_gz,
            ch_metaphlan_sgb2gtbd_file,
            ch_metaphlan_db,
            params.metaphlan_db_version
        ).metaphlan_profiles_merged_tsv
        ch_versions = ch_versions.mix(FASTQ_READ_TAXONOMY_METAPHLAN.out.versions)
    } else {
        ch_read_taxonomy_tsv = Channel.empty()
    }


    /*-----------------------------------------------------------------------------------
        Functional classification: HUMAnN
    -----------------------------------------------------------------------------------*/
    if ( params.run_humann ) {
        // create channel from params.chocophlan_db
        if ( !params.chocophlan_db ) {
            ch_chocophlan_db = null
        } else {
            ch_chocophlan_db = Channel.value( file( params.chocophlan_db, checkIfExists: true ) )
        }

        // create channel from params.uniref_db
        if ( !params.uniref_db ) {
            ch_uniref_db = null
        } else {
            ch_uniref_db = Channel.value( file( params.uniref_db, checkIfExists: true ) )
        }

        // theres probably a better way to handle this. but good enough for me for now..
        if ( !params.run_metaphlan ) {
            error "Error: run_humann is true but run_metaphlan is false. Cannot run HUMAnN without MetaPhlAn."
        }

        //
        // SUBWORKFLOW: HUMAnN
        //
        // TODO double check the metaphlan output channel. not sure its the format i was expecting in the module
        ch_genefamilies_tsv = FASTQ_MICROBIAL_PATHWAY_HUMANN(
            ch_preprocessed_fastq_gz,
            ch_read_taxonomy_tsv,
            ch_chocophlan_db,
            params.chocophlan_db_version,
            ch_uniref_db,
            params.uniref_db_version).humann_genefamilies
        ch_ec_tsv = FASTQ_MICROBIAL_PATHWAY_HUMANN.out.humann_ec
        ch_pathabundance_tsv = FASTQ_MICROBIAL_PATHWAY_HUMANN.out.humann_pathabundance
        ch_pathcoverage_tsv = FASTQ_MICROBIAL_PATHWAY_HUMANN.out.humann_pathcoverage
        ch_versions = ch_versions.mix(FASTQ_MICROBIAL_PATHWAY_HUMANN.out.versions)
    } else {
        ch_genefamilies_tsv = Channel.empty()
        ch_ec_tsv = Channel.empty()
        ch_pathabundance_tsv = Channel.empty()
        ch_pathcoverage_tsv = Channel.empty()
    }


    /*-----------------------------------------------------------------------------------
        Pipeline report utilities
    -----------------------------------------------------------------------------------*/
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    //
    // MODULE: Run MULTIQC to create the pipeline report
    //
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    preprocessed_fastq_gz           = ch_preprocessed_fastq_gz
    preprocessed_read_counts_tsv    = ch_preprocessed_read_counts_tsv
    read_taxonomy_tsv               = ch_read_taxonomy_tsv
    genefamilies_tsv                = ch_genefamilies_tsv
    ec_tsv                          = ch_ec_tsv
    pathabundance_tsv               = ch_pathabundance_tsv
    pathcoverage_tsv                = ch_pathcoverage_tsv
    multiqc_report                  = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions                        = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
