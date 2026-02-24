/**
 * Shared helper functions for nf-denovoslim.
 *
 * Placed in lib/ so the strict syntax parser (v2) doesn't need
 * top-level function definitions in .nf scripts.
 */
class Utils {

    /**
     * Validate that all required pipeline parameters are set (non-null).
     * Call at the top of the workflow block.
     */
    static void validateParams(params) {
        def required = [
            'trinity_fasta',
            'samplesheet',
            'diamond_db',
            'busco_lineage',
            'mmseqs2_swissprot',
            'mmseqs2_pfam',
            'mmseqs2_eggnog',
            'mmseqs2_taxonomy_db'
        ]
        def missing = required.findAll { params[it] == null }
        if (missing) {
            def flags = missing.collect { "--${it}" }.join(', ')
            throw new IllegalArgumentException(
                "Missing required parameter(s): ${flags}"
            )
        }
    }

    /**
     * Extract condition from a sample name.
     * e.g. BMAX56_T4_L  →  T4_L
     */
    static String extractCondition(String sample_name) {
        def parts = sample_name.split('_')
        if (parts.size() < 3) {
            return sample_name
        }
        return "${parts[-2]}_${parts[-1]}"
    }
}
