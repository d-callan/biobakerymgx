process HUMANN_JOIN {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann:3.8--pyh7cba7a3_0':
        'biocontainers/humann:3.8--pyh7cba7a3_0' }"

    input:
    path(input_dir)
    val file_name_pattern

    output:
    path("*_joined.tsv.gz"), emit: joined
    path "versions.yml"    , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    if compgen -G "$input_dir/*$file_name_pattern*.gz" > /dev/null; then
        find $input_dir \( -name '*$file_name_pattern*' \) -exec gunzip --verbose {} \;
    fi
    humann_join_table \\
        --input $input_dir \\
        --output ${file_name_pattern}_joined.tsv \\
        --file_name $file_name_pattern \\
        ${args}
    gzip -n ${file_name_pattern}_joined.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: \$(echo \$(humann --version 2>&1 | sed 's/^.*humann //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    """
    touch ${file_name_pattern}_joined.tsv.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: \$(echo \$(humann --version 2>&1 | sed 's/^.*humann //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
