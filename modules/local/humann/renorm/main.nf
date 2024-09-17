process HUMANN_RENORM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann:3.8--pyh7cba7a3_0':
        'biocontainers/humann:3.8--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(input)
    val units

    output:
    tuple val(meta), path("*_renorm.tsv.gz"), emit: renorm
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [[ $input == *.gz ]]; then
        gunzip -c $input > input.tsv
    else
        mv $input input.tsv
    fi
    humann_renorm_table \\
        --input input.tsv \\
        --output ${prefix}_${units}_renorm.tsv \\
        --units $units \\
        --update-snames \\
        ${args}
    gzip -n ${prefix}_${units}_renorm.tsv

    stub:
    def args = task.ext.args ?: ''
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${units}_renorm.tsv.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: \$(echo \$(humann --version 2>&1 | sed 's/^.*humann //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
