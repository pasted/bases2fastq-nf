#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


/*
 * DISS analysis workflow
 *
 * Processes:
 *   1. DISS_BACKGROUND_ADJUSTED_GUIDE_CALLING
 *   2. DISS_FINAL_GUIDE_ASSIGNMENTS
 *   3. DISS_GUIDE_INTEGRATION
 *   4. DISS_GUIDE_RESCUE_AUDIT
 *   5. DISS_GUIDE_FAMILY_CLASSIFICATION
 *   6. DISS_FINAL_EXPORT
 *
 * Scientific policies:
 *   - P3_SPECIFIC is the primary guide-calling policy.
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
    DISS_BACKGROUND_ADJUSTED_GUIDE_CALLING
} from './modules/local/diss_background_adjusted_guide_calling/main'

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

    def exact_dir =
        "${root}/stage05k_cell_target_matrix"

    def frozen_dir =
        "${root}/stage05k_final_guide_calls"

    def rescue_dir =
        "${root}/stage05l_rescue_validation"

    def specificity_dir =
        "${root}/stage05l_cross_guide_specificity"

    def qc_dir =
        "${root}/stage05l_call_quality_assessment"

    def stability_dir =
        "${root}/stage05l_call_stability_and_calibration"

    def interoperability =
        "${root}/stage05i_interoperability_exports"


    // -------------------------------------------------------------
    // Stage05k exact annotated reads
    // -------------------------------------------------------------

    def exact_reads = file(
        "${exact_dir}/B02_annotated_reads_with_cell_assignment.parquet",
        checkIfExists: true
    )


    // -------------------------------------------------------------
    // Stage05k frozen guide evidence
    // -------------------------------------------------------------

    def frozen_stage05k = file(
        "${frozen_dir}/cell_guide_evidence.tsv",
        checkIfExists: true
    )


    // -------------------------------------------------------------
    // Stage05l rescued reads
    // -------------------------------------------------------------

    def rescued_reads = file(
        "${rescue_dir}/rescue_reads_spatially_annotated.parquet",
        checkIfExists: true
    )


    // -------------------------------------------------------------
    // Stage05l cross-guide specificity
    // -------------------------------------------------------------

    def specificity = file(
        "${specificity_dir}/rescue_all_guide_distances.parquet",
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
    // Process 1: Background-adjusted guide calling
    // =============================================================

    DISS_BACKGROUND_ADJUSTED_GUIDE_CALLING(

        Channel.value([
            meta,
            exact_reads,
            frozen_stage05k,
            rescued_reads,
            specificity
        ])
    )


    // -------------------------------------------------------------
    // Regenerated P0 and P3 evidence
    // -------------------------------------------------------------

    def final_guide_inputs =
        DISS_BACKGROUND_ADJUSTED_GUIDE_CALLING.out.primary_evidence.map {
            meta_out, p0_evidence, p3_evidence ->

            [
                meta_out,
                frozen_stage05k,
                p0_evidence,
                p3_evidence,
                qc_evidence,
                identity,
                fixed_cell_universe
            ]
        }


    // =============================================================
    // Process 2: Final guide assignment export
    // =============================================================

    DISS_FINAL_GUIDE_ASSIGNMENTS(
        final_guide_inputs,
        guide_export_script,
        params.diss_pilot_regression
    )


    // -------------------------------------------------------------
    // Shared final guide assignment outputs
    // -------------------------------------------------------------

    def guide_exports =
        DISS_FINAL_GUIDE_ASSIGNMENTS.out.guide_exports


    // =============================================================
    // Process 3: Guide integration
    // =============================================================

    def integration_inputs = guide_exports.map {
        meta_out, assignments, matrix, status, manifest ->

        [
            meta_out,
            diss_transcripts,
            diss_cell_metadata,
            assignments,
            status,
            manifest
        ]
    }

    DISS_GUIDE_INTEGRATION(
        integration_inputs,
        integration_script,
        params.diss_pilot_regression
    )


    // =============================================================
    // Process 4: Rescue evidence audit
    // =============================================================

    def audit_inputs = guide_exports.map {
        meta_out, assignments, matrix, status, manifest ->

        [
            meta_out,
            assignments,
            matrix
        ]
    }

    DISS_GUIDE_RESCUE_AUDIT(
        audit_inputs,
        audit_script,
        params.diss_pilot_regression
    )


    // =============================================================
    // Process 5: Guide family classification
    // =============================================================

    def integrated =
        DISS_GUIDE_INTEGRATION.out.integrated

    def classification_inputs = integrated.map {
        meta_out, transcripts, cells, assignments, guide_cells ->

        [
            meta_out,
            cells,
            assignments
        ]
    }

    DISS_GUIDE_FAMILY_CLASSIFICATION(
        classification_inputs,
        classifier_script
    )


    // =============================================================
    // Process 6: Final interoperability exports
    // =============================================================

    // Convert metadata maps into scalar join keys.

    def integrated_keyed = integrated.map {
        meta_out, transcripts, cells, assignments, guide_cells ->

        [
            meta_out.id,
            meta_out,
            transcripts,
            cells,
            assignments,
            guide_cells
        ]
    }

    def audit_keyed =
        DISS_GUIDE_RESCUE_AUDIT.out.evidence.map {
            meta_out, evidence ->

            [
                meta_out.id,
                evidence
            ]
        }

    def classified_keyed =
        DISS_GUIDE_FAMILY_CLASSIFICATION.out.classified.map {
            meta_out, family, family_manifest ->

            [
                meta_out.id,
                family
            ]
        }


    // -------------------------------------------------------------
    // Join integrated, audit and classification outputs
    // -------------------------------------------------------------

    def final_export_inputs = integrated_keyed
        .join(audit_keyed)
        .join(classified_keyed)
        .map {
            sample_id,
            meta_out,
            transcripts,
            cells,
            assignments,
            guide_cells,
            evidence,
            family ->

            [
                meta_out,
                transcripts,
                cells,
                assignments,
                guide_cells,
                evidence,
                family,
                diss_feature_metadata,
                diss_cell_metadata,
                diss_target_counts,
                diss_allele_counts
            ]
        }


    DISS_FINAL_EXPORT(
        final_export_inputs,
        export_script,
        params.diss_pilot_regression
    )

}
