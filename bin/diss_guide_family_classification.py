#!/usr/bin/env python3

"""Classify DISS cells by their primary guide gene families.

This stage performs guide-family classification only.

It does not infer editing, assign individual cut sites to guides,
or perform statistical modelling.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


KEY = ["well", "cell_id"]


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)

    return parser.parse_args()


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def sha256(path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def main():
    args = parse_args()

    paths = {
        "cells": args.cells.resolve(),
        "assignments": args.assignments.resolve(),
    }

    for name, path in paths.items():
        require(path.is_file(), f"Missing {name}: {path}")

    input_hashes = {
        name: sha256(path)
        for name, path in paths.items()
    }

    cells = pd.read_parquet(paths["cells"])
    assignments = pd.read_parquet(paths["assignments"])

    require(
        not cells.duplicated(KEY).any(),
        "Duplicate DISS cell keys",
    )

    require(
        not assignments.duplicated(KEY + ["guide"]).any(),
        "Duplicate guide assignment keys",
    )

    require(
        {"guide_annotation_state", "primary_n_called"}.issubset(cells.columns),
        "Missing guide annotation columns",
    )

    require(
        {"guide", "stage05l_p3_called"}.issubset(assignments.columns),
        "Missing guide assignment columns",
    )

    # Only primary P3 calls contribute to guide-family classification.
    primary = assignments.loc[
        assignments["stage05l_p3_called"].astype(bool),
        KEY + ["guide"],
    ].copy()

    primary["guide_gene"] = (
        primary["guide"]
        .astype(str)
        .str.split("_")
        .str[0]
    )

    # Pilot implementation: explicitly preserve the original gene families.
    genes = ["PLCG2", "PTK2B"]

    require(
        primary["guide_gene"].isin(genes).all(),
        "Unexpected guide gene family",
    )

    # Exclude guide-positive cells without DISS observations.
    primary = primary.merge(
        cells[KEY],
        on=KEY,
        how="inner",
        validate="many_to_one",
    )

    counts = (
        primary.groupby(KEY, observed=True)
        .agg(
            n_p3=("guide", "size"),
            n_plcg2=("guide_gene", lambda x: int(x.eq("PLCG2").sum())),
            n_ptk2b=("guide_gene", lambda x: int(x.eq("PTK2B").sum())),
            guide_list=(
                "guide",
                lambda x: ";".join(sorted(x.astype(str))),
            ),
        )
        .reset_index()
    )

    result = cells.merge(
        counts,
        on=KEY,
        how="left",
        validate="one_to_one",
    )

    for column in ["n_p3", "n_plcg2", "n_ptk2b"]:
        result[column] = result[column].fillna(0).astype(int)

    result["guide_list"] = result["guide_list"].fillna("")

    require(
        result["n_p3"].eq(
            result["primary_n_called"].fillna(0).astype(int)
        ).all(),
        "Primary guide multiplicity mismatch",
    )

    result["guide_family_class"] = np.select(
        [
            result.n_p3.eq(0)
            & result.guide_annotation_state.eq("GUIDE_EVIDENCE_UNCALLED"),

            result.n_p3.eq(0)
            & result.guide_annotation_state.eq("NO_GUIDE_EVIDENCE"),

            result.n_p3.eq(1) & result.n_plcg2.eq(1),

            result.n_p3.eq(1) & result.n_ptk2b.eq(1),

            result.n_p3.gt(1)
            & result.n_plcg2.gt(0)
            & result.n_ptk2b.gt(0),

            result.n_p3.gt(1) & result.n_plcg2.gt(0),

            result.n_p3.gt(1) & result.n_ptk2b.gt(0),
        ],
        [
            "GUIDE_EVIDENCE_UNCALLED",
            "NO_GUIDE_EVIDENCE",
            "SINGLE_PLCG2",
            "SINGLE_PTK2B",
            "MULTIPLE_MIXED",
            "MULTIPLE_PLCG2_ONLY",
            "MULTIPLE_PTK2B_ONLY",
        ],
        default="INVALID",
    )

    require(
        not result.guide_family_class.eq("INVALID").any(),
        "Unclassified DISS cells",
    )

    require(
        len(result) == len(cells),
        "DISS cell count changed",
    )

    require(
        result[KEY].equals(cells[KEY]),
        "DISS cell ordering changed",
    )

    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    output = outdir / "diss_cells_guide_family_classes.parquet"

    result.to_parquet(output, index=False)

    for name, path in paths.items():
        require(
            sha256(path) == input_hashes[name],
            f"Input changed: {name}",
        )

    manifest = {
        "stage": "diss_guide_family_classification",
        "inputs": {
            name: {
                "path": str(path),
                "sha256": input_hashes[name],
            }
            for name, path in paths.items()
        },
        "outputs": {
            output.name: {
                "rows": len(result),
                "sha256": sha256(output),
            }
        },
        "guide_families": genes,
        "primary_policy": "P3_SPECIFIC",
        "classification_counts": {
            str(key): int(value)
            for key, value in result.guide_family_class.value_counts().items()
        },
        "editing_status": "NOT_ASSESSED",
    }

    manifest_path = outdir / "guide_family_classification_manifest.json"

    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n"
    )

    print(result.guide_family_class.value_counts().to_string())
    print(f"\nPASS: {len(result):,} DISS cells classified")
    print(f"Output: {output}")


if __name__ == "__main__":
    main()
