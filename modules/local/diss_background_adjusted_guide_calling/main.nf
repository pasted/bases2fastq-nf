process DISS_BACKGROUND_ADJUSTED_GUIDE_CALLING {

    tag "${meta.id}"

    publishDir "${params.diss_outdir}/stage05l_background_adjusted_guide_calling",
        mode: 'copy'

    input:
    tuple val(meta),
          path(exact),
          path(frozen),
          path(rescue),
          path(specificity)

    output:
    tuple val(meta),
          path("P0_EXACT_cell_construct_evidence.parquet"),
          path("P3_SPECIFIC_cell_construct_evidence.parquet"),
          emit: primary_evidence

    path "P1_EXACT_D1_cell_construct_evidence.parquet",
        emit: p1_evidence

    path "P2_EXACT_D1_D2_cell_construct_evidence.parquet",
        emit: p2_evidence

    path "policy_call_comparison.tsv",
        emit: policy_comparison

    path "newly_called_cell_constructs.tsv",
        emit: new_calls

    path "lost_cell_construct_calls.tsv",
        emit: lost_calls

    path "policy_background_fractions.tsv",
        emit: background_fractions

    path "stage05l_calling_summary.tsv",
        emit: calling_summary

    script:
    """
    python ${projectDir}/bin/diss_background_adjusted_guide_calling.py \\
        --exact ${exact} \\
        --frozen ${frozen} \\
        --rescue ${rescue} \\
        --specificity ${specificity} \\
        --outdir .
    """
}
