process DISS_GUIDE_FAMILY_CLASSIFICATION {

    tag "${meta.id}"

    input:
    tuple val(meta), path(cells), path(assignments)
    path classifier_script

    output:
    tuple val(meta),
        path("guide_family/diss_cells_guide_family_classes.parquet"),
        path("guide_family/guide_family_classification_manifest.json"),
        emit: classified

    script:
    """
    python ${classifier_script} \\
        --cells ${cells} \\
        --assignments ${assignments} \\
        --outdir guide_family
    """
}
