#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


/*
 * DISS analysis workflow
 *
 * Stage 16:
 *   1. DISS_FINAL_GUIDE_ASSIGNMENTS
 *   2. DISS_GUIDE_INTEGRATION
 *   3. DISS_GUIDE_RESCUE_AUDIT
 *   4. DISS_GUIDE_FAMILY_CLASSIFICATION
 *   5. DISS_FINAL_EXPORT
 *
 * Inputs:
 *   Stage05k frozen guide evidence
 *   Stage05l P0/P3 guide evidence and QC
 *   Stage05l call identity and fixed cell universe
 *   Stage05i interoperability exports
 *
 * Scientific policies:
 *   - P3_SPECIFIC is the primary guide calling policy.
 *   - Primary guide calls require BH-adjusted q < 0.05.
 *   - FROZEN_ONLY assignments are provenance only.
 *   - Join cells using (well, cell_id).
 *   - Preserve all DISS transcripts and cells.
 *   - Do not annotate unassociated transcripts with guides.
 *   - Do not infer editing from guide assignment.
 */


// ---------------------------------------------------------------------
// Module imports
// ---------------------------------------------------------------------

include {
    DISS_FINAL_GUIDE_ASSIGNMENTS
} from './modules/local/diss_final_guide_assignments/main'

include {
    DISS_GUIDE_INTEGRATION
} from './modules/local/diss_guide_integration/main'

include {
    DISS_GUIDE_RESCUE_AUDIT
} from './modules/local/diss_guide_rescue_evidence_audit/main'

include {
    DISS_GUIDE_FAMILY_CLASSIFICATION
} from './modules/local/diss_guide_family_classification/main'

include {
    DISS_FINAL_EXPORT
} from './modules/local/diss_final_export/main'


// ---------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------

params.diss_analysis_root = null

params.diss_pilot_regression = false

params.diss_outdir = './results/diss_nextflow'


// ---------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------

workflow {

    // -------------------------------------------------------------
    // Validate parameters
    // -------------------------------------------------------------

    if (!params.diss_analysis_root) {
        error "Required parameter: --diss_analysis_root"
    }

    def root = params.diss_analysis_root

    def meta = [
        id: 'DISS'
    ]


    // -------------------------------------------------------------
    // Input directories
    // -------------------------------------------------------------

    def frozen_dir =
        "${root}/stage05k_final_guide_calls"

    def background_dir =
        "${root}/stage05l_background_adjusted_guide_calling"

    def qc_dir =
        "${root}/stage05l_call_quality_assessment"

    def stability_dir =
        "${root}/stage05l_call_stability_and_calibration"

    def interoperability =
        "${root}/stage05i_interoperability_exports"


    // -------------------------------------------------------------
    // Stage05k frozen guide evidence
    // -------------------------------------------------------------

    def frozen_stage05k = file(
        "${frozen_dir}/cell_guide_evidence.tsv",
        checkIfExists: true
    )


    // -------------------------------------------------------------
    // Stage05l guide calling evidence
    // -------------------------------------------------------------

    def p0_evidence = file(
        "${background_dir}/P0_EXACT_cell_construct_evidence.parquet",
        checkIfExists: true
    )

    def p3_evidence = file(
        "${background_dir}/P3_SPECIFIC_cell_construct_evidence.parquet",
        checkIfExists: true
    )


    // -------------------------------------------------------------
    // Stage05l guide QC
    // -------------------------------------------------------------

    def qc_evidence = file(
        "${qc_dir}/called_and_lost_call_evidence.parquet",
        checkIfExists: true
    )


    // -------------------------------------------------------------
    // Stage05l stability and calibration
    // -------------------------------------------------------------

    def identity = file(
        "${stability_dir}/call_identity_by_policy.tsv",
        checkIfExists: true
    )

    def fixed_cell_universe = file(
        "${stability_dir}/fixed_universe_cell_status.tsv",
        checkIfExists: true
    )


    // -------------------------------------------------------------
    // Stage05i DISS inputs
    // -------------------------------------------------------------

    def diss_transcripts = file(
        "${interoperability}/DISS_transcripts_final.parquet",
        checkIfExists: true
    )

    def diss_cell_metadata = file(
        "${interoperability}/DISS_cell_metadata.tsv.gz",
        checkIfExists: true
    )

    def diss_feature_metadata = file(
        "${interoperability}/DISS_feature_metadata.tsv",
        checkIfExists: true
    )

    def diss_target_counts = file(
        "${interoperability}/DISS_cell_target_counts.tsv.gz",
        checkIfExists: true
    )

    def diss_allele_counts = file(
        "${interoperability}/DISS_cell_allele_counts.tsv.gz",
        checkIfExists: true
    )


    // -------------------------------------------------------------
    // Python scripts
    // -------------------------------------------------------------

    def guide_export_script = file(
        "${projectDir}/bin/diss_final_guide_assignment_export.py",
        checkIfExists: true
    )

    def integration_script = file(
        "${projectDir}/bin/diss_guide_integration.py",
        checkIfExists: true
    )

    def audit_script = file(
        "${projectDir}/bin/diss_guide_rescue_evidence_audit.py",
        checkIfExists: true
    )

    def classifier_script = file(
        "${projectDir}/bin/diss_guide_family_classification.py",
        checkIfExists: true
    )

    def export_script = file(
        "${projectDir}/bin/diss_final_interoperability_exports.py",
        checkIfExists: true
    )


    // =============================================================
    // Process 1: Final guide assignment export
    // =============================================================

    DISS_FINAL_GUIDE_ASSIGNMENTS(

        Channel.value([
            meta,
            frozen_stage05k,
            p0_evidence,
            p3_evidence,
            qc_evidence,
            identity,
            fixed_cell_universe
        ]),

        guide_export_script,

        params.diss_pilot_regression
    )


    // -------------------------------------------------------------
    // Reuse generated guide outputs
    // -------------------------------------------------------------

    /*
     * The guide_exports channel contains:
     *
     * meta
     * final_cell_guide_assignments.parquet
     * final_cell_guide_matrix.parquet
     * final_cell_guide_status.tsv
     * final_export_manifest.json
     *
     * Nextflow dataflow channels allow the same output to be consumed
     * independently by integration and audit processes.
     */


    // =============================================================
    // Process 2: DISS/guide integration
    // =============================================================

    def integration_inputs =
        DISS_FINAL_GUIDE_ASSIGNMENTS.out.guide_exports.map {
            guide_meta,
            guide_assignments,
            guide_matrix,
            guide_status,
            guide_manifest ->

            [
                guide_meta,
                diss_transcripts,
                diss_cell_metadata,
                guide_assignments,
                guide_status,
                guide_manifest
            ]
        }


    DISS_GUIDE_INTEGRATION(

        integration_inputs,

        integration_script,

        params.diss_pilot_regression
    )


    // =============================================================
    // Process 3: Guide rescue evidence audit
    // =============================================================

    def audit_inputs =
        DISS_FINAL_GUIDE_ASSIGNMENTS.out.guide_exports.map {
            guide_meta,
            guide_assignments,
            guide_matrix,
            guide_status,
            guide_manifest ->

            [
                guide_meta,
                guide_assignments,
                guide_matrix
            ]
        }


    DISS_GUIDE_RESCUE_AUDIT(

        audit_inputs,

        audit_script,

        params.diss_pilot_regression
    )


    // =============================================================
    // Process 4: Guide family classification
    // =============================================================

    def classification_inputs =
        DISS_GUIDE_INTEGRATION.out.integrated.map {
            integration_meta,
            integrated_transcripts,
            integrated_cells,
            integrated_assignments,
            integrated_guide_cells ->

            [
                integration_meta,
                integrated_cells,
                integrated_assignments
            ]
        }


    DISS_GUIDE_FAMILY_CLASSIFICATION(

        classification_inputs,

        classifier_script
    )


    // =============================================================
    // Process 5: Final interoperability export
    // =============================================================

    // -------------------------------------------------------------
    // Key integration results
    // -------------------------------------------------------------

    def integration_for_export =
        DISS_GUIDE_INTEGRATION.out.integrated.map {
            integration_meta,
            integrated_transcripts,
            integrated_cells,
            integrated_assignments,
            integrated_guide_cells ->

            [
                integration_meta.id,
                integration_meta,
                integrated_transcripts,
                integrated_cells,
                integrated_assignments,
                integrated_guide_cells
            ]
        }


    // -------------------------------------------------------------
    // Key guide family classification results
    // -------------------------------------------------------------

    def classification_for_export =
        DISS_GUIDE_FAMILY_CLASSIFICATION.out.classified.map {
            classification_meta,
            family_file,
            family_manifest ->

            [
                classification_meta.id,
                family_file,
                family_manifest
            ]
        }


    // -------------------------------------------------------------
    // Key guide rescue audit results
    // -------------------------------------------------------------

    def audit_for_export =
        DISS_GUIDE_RESCUE_AUDIT.out.evidence.map {
            audit_meta,
            audit_file ->

            [
                audit_meta.id,
                audit_file
            ]
        }


    // -------------------------------------------------------------
    // Join integration and classification
    // -------------------------------------------------------------

    def integrated_classified =
        integration_for_export.join(
            classification_for_export
        )


    // -------------------------------------------------------------
    // Join rescue audit evidence
    // -------------------------------------------------------------

    def complete_export_inputs =
        integrated_classified.join(
            audit_for_export
        )


    // -------------------------------------------------------------
    // Construct final exporter input tuple
    // -------------------------------------------------------------

    def export_inputs =
        complete_export_inputs.map {
            sample_id,
            export_meta,
            integrated_transcripts,
            integrated_cells,
            integrated_assignments,
            integrated_guide_cells,
            family_file,
            family_manifest,
            audit_file ->

            [
                export_meta,
                integrated_transcripts,
                integrated_cells,
                integrated_assignments,
                integrated_guide_cells,
                audit_file,
                family_file,
                diss_feature_metadata,
                diss_cell_metadata,
                diss_target_counts,
                diss_allele_counts
            ]
        }


    // -------------------------------------------------------------
    // Execute final exporter
    // -------------------------------------------------------------

    DISS_FINAL_EXPORT(

        export_inputs,

        export_script,

        params.diss_pilot_regression
    )

}
