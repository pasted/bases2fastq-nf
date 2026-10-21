process DISS_FINAL_EXPORT {

    tag "${meta.id}"

    cpus 2
    memory '8 GB'
    time '2h'

    input:
    tuple val(meta),
          path(transcripts),
          path(cells),
          path(assignments),
          path(guide_cells),
          path(evidence),
          path(family),
          path(features),
          path(cell_metadata),
          path(target_counts),
          path(allele_counts)

    path exporter_script

    val pilot_regression

    output:
    tuple val(meta), path("final_exports/*"), emit: exports

    script:
    def regression_arg = pilot_regression ? '--pilot-regression' : ''

    """
    mkdir -p analysis/stage05l_diss_guide_integration
    mkdir -p analysis/stage05l_guide_rescue_evidence_audit
    mkdir -p analysis/stage05l_guide_specific_cutsite_detection
    mkdir -p analysis/stage05i_interoperability_exports
    mkdir -p final_exports

    ln -s "\$(realpath ${transcripts})" \
        analysis/stage05l_diss_guide_integration/DISS_transcripts_with_guides.parquet

    ln -s "\$(realpath ${cells})" \
        analysis/stage05l_diss_guide_integration/DISS_cells_with_guides.parquet

    ln -s "\$(realpath ${assignments})" \
        analysis/stage05l_diss_guide_integration/DISS_cell_guide_assignments.parquet

    ln -s "\$(realpath ${guide_cells})" \
        analysis/stage05l_diss_guide_integration/guide_cells_with_DISS.parquet

    ln -s "\$(realpath ${evidence})" \
        analysis/stage05l_guide_rescue_evidence_audit/guide_rescue_assignment_evidence.parquet

    ln -s "\$(realpath ${family})" \
        analysis/stage05l_guide_specific_cutsite_detection/diss_cells_guide_family_classes.parquet

    ln -s "\$(realpath ${features})" \
        analysis/stage05i_interoperability_exports/DISS_feature_metadata.tsv

    ln -s "\$(realpath ${cell_metadata})" \
        analysis/stage05i_interoperability_exports/DISS_cell_metadata.tsv.gz

    ln -s "\$(realpath ${target_counts})" \
        analysis/stage05i_interoperability_exports/DISS_cell_target_counts.tsv.gz

    ln -s "\$(realpath ${allele_counts})" \
        analysis/stage05i_interoperability_exports/DISS_cell_allele_counts.tsv.gz

    python ${exporter_script} \
        --analysis-root analysis \
        --outdir final_exports \
        ${regression_arg}

    test -s final_exports/DISS_export_manifest.json
    test -s final_exports/DISS_export_validation.tsv
    """
}
