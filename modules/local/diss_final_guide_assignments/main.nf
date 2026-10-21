process DISS_FINAL_GUIDE_ASSIGNMENTS {

    tag "${meta.id}"

    input:
    tuple val(meta),
        path(frozen_stage05k),
        path(p0_evidence),
        path(p3_evidence),
        path(qc_evidence),
        path(identity),
        path(fixed_cell_universe)

    path guide_export_script

    val pilot_regression

    output:
    tuple val(meta),
        path("guide_assignments/final_cell_guide_assignments.parquet"),
        path("guide_assignments/final_cell_guide_matrix.parquet"),
        path("guide_assignments/final_cell_guide_status.tsv"),
        path("guide_assignments/final_export_manifest.json"),
        emit: guide_exports

    path "guide_assignments/final_cell_guide_assignments.tsv",
        emit: assignments_tsv

    path "guide_assignments/final_guide_assignment_summary.tsv",
        emit: summary

    path "guide_assignments/final_guide_assignment_qc.tsv",
        emit: qc

    script:
    def pilot_arg = pilot_regression ? '--pilot-regression' : ''

    """
    python ${guide_export_script} \\
        --frozen-stage05k ${frozen_stage05k} \\
        --p0-evidence ${p0_evidence} \\
        --p3-evidence ${p3_evidence} \\
        --qc-evidence ${qc_evidence} \\
        --identity ${identity} \\
        --fixed-cell-universe ${fixed_cell_universe} \\
        --outdir guide_assignments \\
        ${pilot_arg}
    """
}
