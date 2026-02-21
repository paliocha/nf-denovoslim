/**
 * Shared helper functions for nf-denovoslim.
 *
 * Placed in lib/ so the strict syntax parser (v2) doesn't need
 * top-level function definitions in .nf scripts.
 */
class Utils {

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
