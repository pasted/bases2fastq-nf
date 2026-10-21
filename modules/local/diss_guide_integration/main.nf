process DISS_GUIDE_INTEGRATION {

    tag "${meta.id}"

    input:
    tuple val(meta),
        path(diss_transcripts),
        path(diss_cell_metadata),
        path(guide_assignments),
        path(guide_cell_status),
        path(guide_export_manifest)

    path integration_script
    val pilot_regression

    output:
    tuple val(meta),
        path("integration/DISS_transcripts_with_guides.parquet"),
        path("integration/DISS_cells_with_guides.parquet"),
        path("integration/DISS_cell_guide_assignments.parquet"),
        path("integration/guide_cells_with_DISS.parquet"),
        emit: integrated

    path "integration/DISS_transcripts_with_guides.tsv.gz",
        emit: transcripts_tsv

    path "integration/DISS_cells_with_guides.tsv.gz",
        emit: cells_tsv

    path "integration/DISS_guide_integration_by_well.tsv",
        emit: by_well

    path "integration/DISS_guide_integration_manifest.json",
        emit: manifest

    script:
    def pilot_arg = pilot_regression ? '--pilot-regression' : ''

    """
    python ${integration_script} \\
        --diss-transcripts ${diss_transcripts} \\
        --diss-cell-metadata ${diss_cell_metadata} \\
        --guide-assignments ${guide_assignments} \\
        --guide-cell-status ${guide_cell_status} \\
        --guide-export-manifest ${guide_export_manifest} \\
        --outdir integration \\
        ${pilot_arg}
    """
}

