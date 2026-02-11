/*
 * MMseqs2 homology search — reusable for SwissProt and Pfam
 */

process MMSEQS2_SEARCH_SWISSPROT {
    label 'process_high'
    tag "swissprot"

    input:
    path(query_pep)
    path(db, stageAs: 'swissprot_db/*')

    output:
    path("swissprot_alnRes.m8"), emit: m8

    script:
    // The db path points to the directory; the actual DB prefix is the basename
    def db_prefix = db.find { it.name.endsWith('.dbtype') || it.name == params.mmseqs2_swissprot.split('/')[-1] }
    """
    mmseqs easy-search \\
        ${query_pep} \\
        ${params.mmseqs2_swissprot} \\
        swissprot_alnRes.m8 \\
        tmp_sp \\
        -s ${params.mmseqs2_search_sens} \\
        --threads ${task.cpus}
    """
}

process MMSEQS2_SEARCH_PFAM {
    label 'process_high'
    tag "pfam"

    input:
    path(query_pep)
    path(db, stageAs: 'pfam_db/*')

    output:
    path("pfam_alnRes.m8"), emit: m8

    script:
    """
    mmseqs easy-search \\
        ${query_pep} \\
        ${params.mmseqs2_pfam} \\
        pfam_alnRes.m8 \\
        tmp_pfam \\
        -s ${params.mmseqs2_search_sens} \\
        --threads ${task.cpus}
    """
}
