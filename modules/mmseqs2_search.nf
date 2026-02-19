/*
 * MMseqs2 homology search (single-process, for DBs like SwissProt/Pfam)
 * DB path is passed as val to avoid staging multi-GB databases.
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
    def mem_gb = (task.memory.toGiga() * 0.85).intValue()
    """
    mmseqs easy-search \\
        ${query_pep} \\
        ${db_path} \\
        ${tag_name}_alnRes.m8 \\
        tmp_${tag_name} \\
        --split-memory-limit ${mem_gb}G \\
        ${args} \\
        --threads ${task.cpus}
    """
}
