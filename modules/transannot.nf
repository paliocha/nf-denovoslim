/*
 * TransAnnot — functional annotation against SwissProt + Pfam + eggNOG
 *
 * TransAnnot 4.0.0 embeds an annotate.sh that downloads the eggNOG 5.0
 * annotation file and uses AWK logic that assumes e5 naming (no pipes in
 * OG family names). The eggNOG 7 profile DB headers use pipes in family
 * names (e.g. "14|3|3@171637|ApB-57.faa.gz") while the annotation file
 * uses hyphens ("14-3-3@171637|ApB-57.faa.gz"), causing exact matching
 * to fail for ~98% of eggNOG hits.
 *
 * Fix (matching PR soedinglab/transannot#7):
 *   1. Curl wrapper intercepts the e5 download, providing our e7 file
 *   2. After transannot finishes (Pfam + SwissProt correct, eggNOG broken),
 *      we re-run convertalis on the preserved profile search results
 *   3. Apply corrected AWK with pipe→hyphen normalization + .faa.gz stripping
 *   4. Merge corrected eggNOG rows into the final output
 */

process TRANSANNOT {
    label 'process_high'
    tag "${species_label}"

    input:
    path(faa)
    path(eggnog_annot)
    val(species_label)
    val(pfam_db)
    val(eggnog_db)
    val(swissprot_db)

    output:
    path("${species_label}_transannot.tsv"), emit: annotation

    script:
    """
    # ── Curl wrapper: intercept eggNOG5 annotation download ──
    mkdir -p fake_bin
    cat > fake_bin/curl << 'WRAPPER'
#!/bin/bash
for arg in "\$@"; do
    if echo "\$arg" | grep -q "eggnog.*og_annotations"; then
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
exec /usr/bin/curl "\$@"
WRAPPER
    chmod +x fake_bin/curl
    export PATH="\$(pwd)/fake_bin:\$PATH"

    # ── Run transannot (Pfam + SwissProt correct; eggNOG broken) ──
    transannot createquerydb ${faa} queryDB tmp_createdb

    transannot annotate \\
        queryDB \\
        ${pfam_db} \\
        ${eggnog_db} \\
        ${swissprot_db} \\
        ${species_label}_transannot_raw.tsv \\
        tmp_annotate \\
        --no-run-clust \\
        --threads ${task.cpus}

    # ── Fix eggNOG v7 annotation ──
    # The embedded annotate.sh's AWK does exact matching on column 5
    # (theader from profile search) which has pipes in family names.
    # Our annotation file has hyphens → no match. We redo this step
    # with proper normalization: strip .faa.gz and convert pipes→hyphens
    # on both the annotation keys and the search result headers.

    # Locate preserved MMseqs2 search results
    ANNOT_TMP=\$(find tmp_annotate -name "prof2_searchDB.dbtype" -exec dirname {} \\; | head -1)

    if [ -n "\$ANNOT_TMP" ] && [ -f "\$ANNOT_TMP/prof2_searchDB.dbtype" ]; then
        echo "Re-running eggNOG annotation with corrected pipe→hyphen AWK..."

        # Re-generate profile search TSV (was deleted by annotate.sh)
        transannot convertalis \\
            "\$ANNOT_TMP/clu_rep" \\
            ${eggnog_db} \\
            "\$ANNOT_TMP/prof2_searchDB" \\
            prof2_regen.tsv \\
            --format-output "query,target,qstart,qend,theader,evalue,pident,bits" \\
            --format-mode 4

        # Build OG→family mapping with normalization:
        #   annotation col 2: "14-3-3@171637|ApB-57.faa.gz"
        #   → strip .faa.gz → all pipes→hyphens → key: "14-3-3@171637-ApB-57"
        awk -F'\\t' -v OFS='\\t' '{
            key = \$2
            sub(/\\.faa\\.gz\$/, "", key)
            gsub(/\\|/, "-", key)
            print key, \$4
        }' ${eggnog_annot} > eggnog_mapping.tsv

        # Match profile search results with same normalization:
        #   theader col 5: "14|3|3@171637|ApB-57.faa.gz"
        #   → strip .faa.gz → all pipes→hyphens → key: "14-3-3@171637-ApB-57"
        awk -F'\\t' -v OFS='\\t' '
            NR==FNR { clr[\$1] = \$2; next }
            {
                key = \$5
                sub(/\\.faa\\.gz\$/, "", key)
                gsub(/\\|/, "-", key)
                if (key in clr) { \$5 = clr[key]; print }
            }
        ' eggnog_mapping.tsv prof2_regen.tsv | \\
            LC_ALL=C sort -s -k1b,1 | \\
            awk -F'\\t' -v OFS='\\t' '{ \$(NF+1) = "seq-prof search"; print }' | \\
            awk -F'\\t' -v OFS='\\t' '{ \$(NF+1) = "eggNOG"; print }' > eggnog_corrected.tsv

        EGGNOG_HITS=\$(wc -l < eggnog_corrected.tsv)
        echo "eggNOG v7: \$EGGNOG_HITS hits with corrected mapping"

        # Merge: header + Pfam/SwissProt (remove any broken eggNOG) + corrected eggNOG
        head -1 ${species_label}_transannot_raw.tsv > ${species_label}_transannot.tsv
        tail -n +2 ${species_label}_transannot_raw.tsv | \\
            awk -F'\\t' '\$NF != "eggNOG"' >> ${species_label}_transannot.tsv
        cat eggnog_corrected.tsv >> ${species_label}_transannot.tsv

        echo "Final annotation: \$(tail -n +2 ${species_label}_transannot.tsv | wc -l) total hits"
    else
        echo "WARNING: prof2_searchDB not found — using raw output (eggNOG may be incomplete)"
        cp ${species_label}_transannot_raw.tsv ${species_label}_transannot.tsv
    fi
    """
}
