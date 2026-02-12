/*
 * MMseqs2 homology search — parameterized for any target database
 *
 * Usage with DSL2 aliasing:
 *   include { MMSEQS2_SEARCH as MMSEQS2_SEARCH_SWISSPROT } from './modules/mmseqs2_search'
 *   include { MMSEQS2_SEARCH as MMSEQS2_SEARCH_PFAM      } from './modules/mmseqs2_search'
 *
 * The database path is passed as val() (not path()) to avoid staging
 * multi-GB pre-built databases to scratch — they are read directly from
 * the shared filesystem.
 */

process MMSEQS2_SEARCH {
    tag "${tag_name}"

    input:
    path(query_pep)
    val(db_path)
    val(tag_name)

    output:
    path("${tag_name}_alnRes.m8"), emit: m8

    script:
    def args = task.ext.args ?: ''
    """
    mmseqs easy-search \\
        ${query_pep} \\
        ${db_path} \\
        ${tag_name}_alnRes.m8 \\
        tmp_${tag_name} \\
        -s ${params.mmseqs2_search_sens} \\
        ${args} \\
        --threads ${task.cpus}
    """
}
