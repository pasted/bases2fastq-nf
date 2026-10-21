process DISS_GUIDE_RESCUE_AUDIT {

    tag "${meta.id}"

    input:
    tuple val(meta), path(assignments), path(matrix)
    path audit_script
    val pilot_regression

    output:
    tuple val(meta),
        path("audit/guide_rescue_assignment_evidence.parquet"),
        emit: evidence

    path "audit/guide_rescue_evidence_summary.tsv",
        emit: summary

    path "audit/guide_rescue_evidence_flags.tsv",
        emit: flags

    path "audit/rescued_guide_evidence_strata.tsv",
        emit: strata

    path "audit/guide_rescue_evidence_audit_manifest.json",
        emit: manifest

    script:
    def pilot_arg = pilot_regression ? '--pilot-regression' : ''

    """
    python ${audit_script} \\
        --assignments ${assignments} \\
        --matrix ${matrix} \\
        --outdir audit \\
        ${pilot_arg}
    """
}
