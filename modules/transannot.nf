/*
 * TransAnnot — functional annotation against SwissProt + Pfam + eggNOG
 *
 * TransAnnot 4.0.0 hardcodes a curl download of the eggNOG 5.0 annotation
 * file (e5.og_annotations.tsv) for OG→description mapping. To support
 * eggNOG 7, we inject a pre-converted annotation file via a curl wrapper
 * that intercepts the eggnog5 URL and returns the local file instead.
 */

process TRANSANNOT {
    label 'process_high'
    tag "${species_label}"

    publishDir "${params.outdir}/transannot", mode: 'copy'

    input:
    path(faa)
    path(eggnog_annot)
    val(species_label)

    output:
    path("${species_label}_transannot.tsv"), emit: annotation

    script:
    """
    # Create a curl wrapper that intercepts the eggNOG annotation download
    # and returns our pre-converted eggNOG7 file instead
    mkdir -p fake_bin
    cat > fake_bin/curl << 'WRAPPER'
#!/bin/bash
# Intercept eggNOG annotation download, pass everything else through
for arg in "\$@"; do
    if echo "\$arg" | grep -q "eggnog.*og_annotations"; then
        # Find the -o output path from the arguments
        outfile=""
        next_is_out=false
        for a in "\$@"; do
            if \$next_is_out; then outfile="\$a"; break; fi
            if [ "\$a" = "-o" ]; then next_is_out=true; fi
        done
        if [ -n "\$outfile" ]; then
            cp "${eggnog_annot}" "\$outfile"
            exit 0
        fi
    fi
done
# Fall through to real curl for anything else
exec /usr/bin/curl "\$@"
WRAPPER
    chmod +x fake_bin/curl
    export PATH="\$(pwd)/fake_bin:\$PATH"

    # Create MMseqs2 query database from the final .faa
    transannot createquerydb ${faa} queryDB tmp_createdb

    # Annotate against all 3 databases
    transannot annotate \\
        queryDB \\
        ${params.mmseqs2_pfam} \\
        ${params.mmseqs2_eggnog} \\
        ${params.mmseqs2_swissprot} \\
        ${species_label}_transannot.tsv \\
        tmp_annotate \\
        --no-run-clust \\
        --threads ${task.cpus}
    """
}
