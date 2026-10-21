#!/usr/bin/env python3

"""
Stage16g: Read-only evidence audit of P3 rescued guide assignments.

Refactored from:
    20261015_14_stage05l_guide_rescue_evidence_audit.py

Purpose
-------
Generate descriptive QC evidence for final P3 guide assignments.

This is NOT a new guide-calling policy.

The audit does not:
    - Modify existing guide assignments.
    - Change the frozen P3 q < 0.05 calling threshold.
    - Reclassify rescued assignments.
    - Interpret DISS cut-site mismatches as guide-rescue evidence.

Inputs
------
1. final_cell_guide_assignments.parquet
2. final_cell_guide_matrix.parquet

Outputs
-------
1. guide_rescue_assignment_evidence.parquet
2. guide_rescue_evidence_summary.tsv
3. guide_rescue_evidence_flags.tsv
4. rescued_guide_evidence_strata.tsv
5. guide_rescue_evidence_audit_manifest.json

Usage
-----
python scripts/20261020_16g_guide_rescue_evidence_audit.py \
    --assignments path/to/final_cell_guide_assignments.parquet \
    --matrix path/to/final_cell_guide_matrix.parquet \
    --outdir path/to/output_directory \
    --pilot-regression

The --pilot-regression option enforces the original Teton Atlas
pilot dimensions and assignment-status counts.
"""

import argparse
import hashlib
import json

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# Constants and configuration
# ---------------------------------------------------------------------

KEY = ["well_id", "cell_id", "guide"]

CELL_KEY = ["well_id", "cell_id"]

EXPECTED = {
    "CONSENSUS_CALL": 7327,
    "RESCUED_CALL": 4265,
    "FROZEN_ONLY": 176,
}

EXPECTED_ASSIGNMENTS = 11768

EXPECTED_MATRIX_ROWS = 112910

P3_Q_THRESHOLD = 0.05


# ---------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Read-only evidence audit of final P3 guide assignments."
        )
    )

    parser.add_argument(
        "--assignments",
        type=Path,
        required=True,
        help=(
            "Input final_cell_guide_assignments.parquet."
        ),
    )

    parser.add_argument(
        "--matrix",
        type=Path,
        required=True,
        help=(
            "Input final_cell_guide_matrix.parquet."
        ),
    )

    parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Output directory.",
    )

    parser.add_argument(
        "--pilot-regression",
        action="store_true",
        help=(
            "Enforce validated Teton Atlas pilot dimensions "
            "and assignment-status counts."
        ),
    )

    return parser.parse_args()


# ---------------------------------------------------------------------
# General utilities
# ---------------------------------------------------------------------

def require(test, message):

    if not bool(test):
        raise RuntimeError(message)


def sha256(path):

    h = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            h.update(block)

    return h.hexdigest()


def save(df, outdir, name):

    path = outdir / name

    if path.suffix == ".parquet":
        df.to_parquet(path, index=False)
    else:
        df.to_csv(
            path,
            sep="\t",
            index=False,
        )

    print(
        f"{name:55s} "
        f"{len(df):>9,d} rows"
    )

    return path


# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------

def validate_inputs(a, m, pilot_regression):

    required = set(
        KEY
        + [
            "assignment_status",
            "guide_reads",
            "exact_reads",
            "d1_reads",
            "d2_reads",
            "rescue_reads",
            "rescue_fraction",
            "d2_fraction",
            "combined_both_ends",
            "exact_both_ends",
            "is_unique_top",
            "binomial_q",
            "guide_fraction_in_cell",
            "qc_category",
            "stage05l_p3_called",
        ]
    )

    missing = sorted(required - set(a.columns))

    require(
        not missing,
        f"Missing assignment columns: {missing}",
    )

    require(
        set(KEY + ["guide_reads", "guide_called"]).issubset(m.columns),
        "Missing P3 matrix evidence columns",
    )

    require(
        not a.duplicated(KEY).any(),
        "Duplicate assignment keys",
    )

    require(
        not m.duplicated(KEY).any(),
        "Duplicate matrix keys",
    )

    if pilot_regression:

        require(
            len(a) == EXPECTED_ASSIGNMENTS,
            (
                "Unexpected pilot assignment dimensions: "
                f"{len(a)} != {EXPECTED_ASSIGNMENTS}"
            ),
        )

        require(
            len(m) == EXPECTED_MATRIX_ROWS,
            (
                "Unexpected pilot matrix dimensions: "
                f"{len(m)} != {EXPECTED_MATRIX_ROWS}"
            ),
        )

    counts = a.assignment_status.value_counts()

    if pilot_regression:

        for status, expected_n in EXPECTED.items():

            observed_n = int(counts.get(status, 0))

            require(
                observed_n == expected_n,
                (
                    f"Unexpected pilot {status} count: "
                    f"{observed_n} != {expected_n}"
                ),
            )

    require(
        a[KEY]
        .merge(
            m[KEY],
            on=KEY,
            how="left",
            indicator=True,
        )
        ._merge.eq("both")
        .all(),
        "Assignment missing from complete guide matrix",
    )

    require(
        (
            a[["exact_reads", "d1_reads", "d2_reads"]]
            .sum(axis=1)
            == a.guide_reads
        ).all(),
        "Read-origin conservation failed",
    )

    require(
        (
            a.assignment_status.eq("FROZEN_ONLY")
            == ~a.stage05l_p3_called
        ).all(),
        "P3 call/status inconsistency",
    )

    require(
        (
            a.loc[
                a.stage05l_p3_called,
                "binomial_q",
            ]
            < P3_Q_THRESHOLD
        ).all(),
        "Significant P3 assignment has q >= 0.05",
    )

    return counts


# ---------------------------------------------------------------------
# Independent guide-read ranking
# ---------------------------------------------------------------------

def calculate_guide_ranks(a, m):

    """
    Independently calculate within-cell guide-read ranks.

    Ranking uses:
        well_id
        cell_id
        guide_reads descending
        guide ascending

    Read-rank uniqueness is not sequence-mapping uniqueness.
    """

    sorted_m = m.sort_values(
        [
            "well_id",
            "cell_id",
            "guide_reads",
            "guide",
        ],
        ascending=[
            True,
            True,
            False,
            True,
        ],
    ).copy()

    sorted_m["rank_by_reads"] = (
        sorted_m
        .groupby(CELL_KEY)
        .cumcount()
        + 1
    )

    ranked = sorted_m[
        KEY + ["rank_by_reads"]
    ].copy()

    max_reads = (
        m.groupby(CELL_KEY)
        .guide_reads
        .transform("max")
    )

    m["n_tied_top"] = (
        m.guide_reads.eq(max_reads)
        .groupby(
            [
                m.well_id,
                m.cell_id,
            ]
        )
        .transform("sum")
    )

    top = (
        m[
            CELL_KEY + ["n_tied_top"]
        ]
        .drop_duplicates(CELL_KEY)
    )

    a = (
        a.merge(
            ranked,
            on=KEY,
            validate="one_to_one",
        )
        .merge(
            top,
            on=CELL_KEY,
            validate="many_to_one",
        )
    )

    return a


# ---------------------------------------------------------------------
# Assignment-level evidence flags
# ---------------------------------------------------------------------

def calculate_evidence_flags(a):

    a["exact_support"] = (
        a.exact_reads.gt(0)
    )

    a["mismatch_only_support"] = (
        a.exact_reads.eq(0)
        & a.rescue_reads.gt(0)
    )

    a["both_end_support"] = (
        a.combined_both_ends.astype(bool)
    )

    a["unique_read_top"] = (
        a.n_tied_top.eq(1)
        & a.rank_by_reads.eq(1)
    )

    a["near_q_threshold"] = (
        a.binomial_q.ge(0.01)
        & a.binomial_q.lt(P3_Q_THRESHOLD)
    )

    a["evidence_flag"] = np.select(
        [
            ~a.both_end_support,
            ~a.unique_read_top,
            a.mismatch_only_support,
            a.near_q_threshold,
        ],
        [
            "NO_BOTH_END_SUPPORT",
            "NOT_UNIQUE_READ_TOP",
            "MISMATCH_ONLY_READ_SUPPORT",
            "Q_NEAR_THRESHOLD",
        ],
        default="NO_LISTED_FLAG",
    )

    # Descriptive evidence only.
    # These flags do not modify frozen P3 guide calls.

    return a


# ---------------------------------------------------------------------
# Evidence summary tables
# ---------------------------------------------------------------------

def make_evidence_summary(a):

    dimensions = [
        "assignment_status",
        "guide",
        "well_id",
        "qc_category",
    ]

    rows = []

    for dimension in dimensions:

        for (status, group), sub in a.groupby(
            ["assignment_status", dimension],
            dropna=False,
        ):

            rows.append(
                {
                    "dimension": dimension,
                    "assignment_status": status,
                    "group": str(group),
                    "assignments": len(sub),
                    "exact_support_n": int(
                        sub.exact_support.sum()
                    ),
                    "mismatch_only_n": int(
                        sub.mismatch_only_support.sum()
                    ),
                    "both_end_n": int(
                        sub.both_end_support.sum()
                    ),
                    "unique_read_top_n": int(
                        sub.unique_read_top.sum()
                    ),
                    "near_q_threshold_n": int(
                        sub.near_q_threshold.sum()
                    ),
                    "median_guide_reads": float(
                        sub.guide_reads.median()
                    ),
                    "median_q": (
                        float(sub.binomial_q.median())
                        if sub.binomial_q.notna().any()
                        else np.nan
                    ),
                    "median_rescue_fraction": float(
                        sub.rescue_fraction.median()
                    ),
                }
            )

    return pd.DataFrame(rows)


def make_evidence_flag_summary(a):

    flags = (
        a.groupby(
            [
                "assignment_status",
                "evidence_flag",
            ],
            dropna=False,
        )
        .size()
        .rename("assignments")
        .reset_index()
    )

    return flags


# ---------------------------------------------------------------------
# Rescued-assignment evidence strata
# ---------------------------------------------------------------------

def make_rescue_strata(a, pilot_regression):

    rescue = a.loc[
        a.assignment_status.eq("RESCUED_CALL")
    ].copy()

    if pilot_regression:

        require(
            len(rescue) == EXPECTED["RESCUED_CALL"],
            "Pilot rescue count changed",
        )

    # Empirical evidence strata.
    # These are not calibrated error probabilities.

    rescue["read_support_bin"] = pd.cut(
        rescue.guide_reads,
        bins=[
            -1,
            1,
            2,
            4,
            9,
            np.inf,
        ],
        labels=[
            "0-1",
            "2",
            "3-4",
            "5-9",
            "10+",
        ],
    )

    rescue["rescue_fraction_bin"] = pd.cut(
        rescue.rescue_fraction,
        bins=[
            -0.001,
            0.25,
            0.5,
            0.75,
            1.001,
        ],
        labels=[
            "0-25%",
            "25-50%",
            "50-75%",
            "75-100%",
        ],
        include_lowest=True,
    )

    strat = (
        rescue.groupby(
            [
                "read_support_bin",
                "rescue_fraction_bin",
            ],
            observed=True,
            dropna=False,
        )
        .agg(
            assignments=("guide", "size"),
            both_end=("both_end_support", "sum"),
            unique_read_top=("unique_read_top", "sum"),
            mismatch_only=("mismatch_only_support", "sum"),
        )
        .reset_index()
    )

    return strat


# ---------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------

def write_manifest(
    outdir,
    input_hashes,
    assignment_counts,
    files,
    pilot_regression,
):

    outputs = {}

    for name, path in files.items():

        if path.suffix == ".parquet":
            df = pd.read_parquet(path)
        else:
            df = pd.read_csv(
                path,
                sep="\t",
            )

        outputs[name] = {
            "rows": len(df),
            "sha256": sha256(path),
        }

    manifest = {
        "stage": "stage05l_step14_guide_rescue_evidence_audit",
        "created_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "input_sha256": input_hashes,
        "assignment_counts": {
            str(status): int(count)
            for status, count in assignment_counts.items()
        },
        "pilot_regression": pilot_regression,
        "outputs": outputs,
        "limitations": [
            (
                "No independent truth labels: no validated "
                "sensitivity, specificity or FDR."
            ),
            (
                "Read-rank uniqueness is not "
                "sequence-mapping uniqueness."
            ),
            (
                "Evidence flags do not change "
                "the frozen P3 q<0.05 policy."
            ),
            (
                "Rescue refers to guide evidence, "
                "not DISS cut-site mismatches."
            ),
        ],
    }

    manifest_path = (
        outdir
        / "guide_rescue_evidence_audit_manifest.json"
    )

    manifest_path.write_text(
        json.dumps(
            manifest,
            indent=2,
        )
        + "\n"
    )

    return manifest_path


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    args = parse_args()

    source = args.assignments.resolve()
    matrix = args.matrix.resolve()
    outdir = args.outdir.resolve()

    require(
        source.is_file(),
        f"Missing input: {source}",
    )

    require(
        matrix.is_file(),
        f"Missing input: {matrix}",
    )

    require(
        source != matrix,
        "Assignments and matrix must be different files",
    )

    before = {
        str(path): sha256(path)
        for path in [source, matrix]
    }

    print("=" * 72)
    print("Stage16g: Guide-rescue evidence audit")
    print("=" * 72)

    print(f"Assignments: {source}")
    print(f"Matrix:      {matrix}")
    print(f"Output:      {outdir}")
    print(f"Pilot check: {args.pilot_regression}")

    a = pd.read_parquet(source)
    m = pd.read_parquet(matrix)

    # -------------------------------------------------------------
    # Validate inputs and frozen P3 calling policy
    # -------------------------------------------------------------

    counts = validate_inputs(
        a,
        m,
        args.pilot_regression,
    )

    print("\nAssignment status counts:")

    print(
        counts.to_string()
    )

    # -------------------------------------------------------------
    # Independent guide-read ranking
    # -------------------------------------------------------------

    a = calculate_guide_ranks(a, m)

    # -------------------------------------------------------------
    # Descriptive evidence flags
    # -------------------------------------------------------------

    a = calculate_evidence_flags(a)

    # -------------------------------------------------------------
    # Generate outputs
    # -------------------------------------------------------------

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    files = {}

    files[
        "guide_rescue_assignment_evidence.parquet"
    ] = save(
        a.sort_values(KEY),
        outdir,
        "guide_rescue_assignment_evidence.parquet",
    )

    files[
        "guide_rescue_evidence_summary.tsv"
    ] = save(
        make_evidence_summary(a),
        outdir,
        "guide_rescue_evidence_summary.tsv",
    )

    files[
        "guide_rescue_evidence_flags.tsv"
    ] = save(
        make_evidence_flag_summary(a),
        outdir,
        "guide_rescue_evidence_flags.tsv",
    )

    files[
        "rescued_guide_evidence_strata.tsv"
    ] = save(
        make_rescue_strata(
            a,
            args.pilot_regression,
        ),
        outdir,
        "rescued_guide_evidence_strata.tsv",
    )

    # -------------------------------------------------------------
    # Confirm input files have not changed
    # -------------------------------------------------------------

    for path, original_hash in before.items():

        require(
            sha256(Path(path)) == original_hash,
            f"Input changed: {path}",
        )

    # -------------------------------------------------------------
    # Write provenance manifest
    # -------------------------------------------------------------

    manifest_path = write_manifest(
        outdir=outdir,
        input_hashes=before,
        assignment_counts=counts,
        files=files,
        pilot_regression=args.pilot_regression,
    )

    print("\nManifest:", manifest_path)

    print(
        "\nPASS: assignment counts, read-origin conservation, "
        "evidence audit, unchanged input hashes"
    )

    print("Outputs:", outdir)


if __name__ == "__main__":
    main()
