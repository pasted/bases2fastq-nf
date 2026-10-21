#!/usr/bin/env python3

"""Stage16i: validated, provenance-preserving P3 guide assignment export.

Refactored from:
    20261015_07_stage05l_final_guide_assignment_export.py

Purpose
-------
Export final P3_SPECIFIC guide assignments while preserving the original
Stage05k P0_EXACT calls as provenance.

Scientific policy
-----------------
- P3_SPECIFIC is the primary guide-assignment policy.
- Primary guide calls require within-well BH-adjusted binomial_q < 0.05.
- Frozen Stage05k calls are retained without modifying their statistics.
- FROZEN_ONLY assignments are provenance records, not primary P3 calls.
- QC categories are descriptive, not additional statistical thresholds.
- The complete P3 cell-by-guide matrix is retained.

Validation
----------
General scientific invariants are always enforced.

Pilot-specific expected counts are checked only when --pilot-regression
is supplied.

All input files are treated as read-only and checked using SHA256 hashes.

The script is independent of the current working directory and is
compatible with standalone execution and Nextflow.
"""

import argparse
import hashlib
import json

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------

KEY = ["well_id", "cell_id", "guide"]

CELL = ["well_id", "cell_id"]

POLICIES = [
    "P0_EXACT",
    "P1_EXACT_D1",
    "P2_EXACT_D1_D2",
    "P3_SPECIFIC",
]

GUIDES = [
    "PTK2B_A",
    "PTK2B_B",
    "PTK2B_C",
    "PTK2B_D",
    "PLCG2_A",
    "PLCG2_B",
    "PLCG2_C",
]

EXPECTED = {
    "frozen": 7503,
    "p3": 11592,
    "consensus": 7327,
    "rescued": 4265,
    "frozen_only": 176,
    "union": 11768,
    "p3_cells": 10114,
    "cell_universe": 16133,
    "p3_matrix": 112910,
}

EXPECTED_FROZEN_MATRIX = 103250


# ---------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Stage16i: validated P3_SPECIFIC final guide assignment export"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--frozen-stage05k",
        type=Path,
        required=True,
        help="Frozen Stage05k cell_guide_evidence.tsv",
    )

    parser.add_argument(
        "--p0-evidence",
        type=Path,
        required=True,
        help="P0_EXACT_cell_construct_evidence.parquet",
    )

    parser.add_argument(
        "--p3-evidence",
        type=Path,
        required=True,
        help="P3_SPECIFIC_cell_construct_evidence.parquet",
    )

    parser.add_argument(
        "--qc-evidence",
        type=Path,
        required=True,
        help="Stage05l called_and_lost_call_evidence.parquet",
    )

    parser.add_argument(
        "--identity",
        type=Path,
        required=True,
        help="Stage05l call_identity_by_policy.tsv",
    )

    parser.add_argument(
        "--fixed-cell-universe",
        type=Path,
        required=True,
        help="Stage05l fixed_universe_cell_status.tsv",
    )

    parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Output directory for final guide assignment tables",
    )

    parser.add_argument(
        "--pilot-regression",
        action="store_true",
        help="Enforce exact expected counts from the Teton Atlas LO pilot",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------

def require(condition, message):

    if not bool(condition):
        raise RuntimeError(message)


def digest(path):

    h = hashlib.sha256()

    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)

    return h.hexdigest()


def write(df, name, outdir):

    path = outdir / name

    if path.suffix == ".parquet":
        df.to_parquet(path, index=False)
    else:
        df.to_csv(
            path,
            sep="\t",
            index=False,
            na_rep="",
        )

    print(
        f"{name:48s} {len(df):>9,d} rows",
        flush=True,
    )

    return path


def validate_keys(df, label):

    require(
        df[KEY].notna().all().all(),
        f"{label}: null assignment key",
    )

    require(
        not df.duplicated(KEY).any(),
        f"{label}: duplicate assignment key",
    )

    require(
        df.guide.isin(GUIDES).all(),
        f"{label}: unknown construct",
    )


def status(n):

    return np.select(
        [n == 0, n == 1],
        ["UNCALLED", "SINGLE"],
        default="MULTIPLE",
    )


def check_pilot_count(name, observed, expected, enabled):

    """Apply an exact pilot-specific count check when requested."""

    if not enabled:
        return

    require(
        int(observed) == int(expected),
        (
            f"Pilot regression failed: {name}: "
            f"observed={observed}, expected={expected}"
        ),
    )


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    args = parse_args()

    outdir = args.outdir.resolve()

    paths = {
        "frozen_stage05k": args.frozen_stage05k.resolve(),
        "p0_evidence": args.p0_evidence.resolve(),
        "p3_evidence": args.p3_evidence.resolve(),
        "qc_evidence": args.qc_evidence.resolve(),
        "identity": args.identity.resolve(),
        "fixed_cell_universe": args.fixed_cell_universe.resolve(),
    }

    print("=" * 80)
    print("Stage16i: Final P3_SPECIFIC guide assignment export")
    print("=" * 80)

    print(f"Output directory: {outdir}")
    print(f"Pilot regression: {args.pilot_regression}")

    print("\nInput files")

    for label, path in paths.items():

        require(
            path.is_file(),
            f"Missing {label}: {path}",
        )

        print(f"  {label:24s} {path}")

    # -------------------------------------------------------------
    # Input provenance
    # -------------------------------------------------------------

    input_hashes = {
        name: digest(path)
        for name, path in paths.items()
    }

    # -------------------------------------------------------------
    # Load input tables
    # -------------------------------------------------------------

    frozen = pd.read_csv(
        paths["frozen_stage05k"],
        sep="\t",
    )

    p0 = pd.read_parquet(paths["p0_evidence"])

    p3 = pd.read_parquet(paths["p3_evidence"])

    qc = pd.read_parquet(paths["qc_evidence"])

    identity = pd.read_csv(
        paths["identity"],
        sep="\t",
    )

    fixed = pd.read_csv(
        paths["fixed_cell_universe"],
        sep="\t",
    )

    print("\nInput rows")

    for name, df in [
        ("frozen", frozen),
        ("p0", p0),
        ("p3", p3),
        ("qc", qc),
        ("identity", identity),
        ("fixed", fixed),
    ]:

        print(f"  {name:16s} {len(df):>10,d}")

    # -------------------------------------------------------------
    # Validate assignment keys
    # -------------------------------------------------------------

    for name, df in [
        ("frozen", frozen),
        ("p0", p0),
        ("p3", p3),
        ("identity", identity),
    ]:

        validate_keys(df, name)

    require(
        not fixed.duplicated(CELL).any(),
        "Fixed universe has duplicate cells",
    )

    require(
        fixed[CELL].notna().all().all(),
        "Fixed universe contains null cell keys",
    )

    # -------------------------------------------------------------
    # Pilot-specific input checks
    # -------------------------------------------------------------

    check_pilot_count(
        "cell_universe",
        len(fixed),
        EXPECTED["cell_universe"],
        args.pilot_regression,
    )

    check_pilot_count(
        "p3_matrix",
        len(p3),
        EXPECTED["p3_matrix"],
        args.pilot_regression,
    )

    check_pilot_count(
        "frozen_matrix",
        len(frozen),
        EXPECTED_FROZEN_MATRIX,
        args.pilot_regression,
    )

    check_pilot_count(
        "p0_matrix",
        len(p0),
        EXPECTED_FROZEN_MATRIX,
        args.pilot_regression,
    )

    check_pilot_count(
        "frozen_calls",
        frozen.guide_called.sum(),
        EXPECTED["frozen"],
        args.pilot_regression,
    )

    check_pilot_count(
        "p0_calls",
        p0.guide_called.sum(),
        EXPECTED["frozen"],
        args.pilot_regression,
    )

    check_pilot_count(
        "p3_calls",
        p3.guide_called.sum(),
        EXPECTED["p3"],
        args.pilot_regression,
    )

    # -------------------------------------------------------------
    # General P3 scientific invariants
    # -------------------------------------------------------------

    require(
        (p3.groupby(CELL).size() == len(GUIDES)).all(),
        "Incomplete P3 seven-guide matrix",
    )

    require(
        p3.guide_called.notna().all(),
        "Null P3 calls",
    )

    require(
        (
            p3.guide_called
            == (p3.binomial_q < 0.05).fillna(False)
        ).all(),
        "P3 guide_called inconsistent with q<0.05",
    )

    # -------------------------------------------------------------
    # Frozen P0 reconciliation
    # -------------------------------------------------------------

    p0check = p0[
        KEY + [
            "guide_called",
            "guide_reads",
            "binomial_p",
            "binomial_q",
        ]
    ].merge(
        frozen[
            KEY + [
                "guide_called",
                "guide_reads",
                "binomial_p",
                "binomial_q",
            ]
        ],
        on=KEY,
        how="outer",
        validate="one_to_one",
        indicator=True,
        suffixes=("_p0", "_frozen"),
    )

    require(
        p0check._merge.eq("both").all(),
        "P0 and frozen keys differ",
    )

    for col in ["guide_called", "guide_reads"]:

        require(
            p0check[f"{col}_p0"].eq(
                p0check[f"{col}_frozen"]
            ).all(),
            f"P0 and frozen {col} differ",
        )

    for col in ["binomial_p", "binomial_q"]:

        require(
            np.allclose(
                p0check[f"{col}_p0"],
                p0check[f"{col}_frozen"],
                equal_nan=True,
                rtol=1e-8,
                atol=1e-12,
            ),
            f"P0 and frozen {col} differ",
        )

    print("\nPASS: Frozen Stage05k / P0 reconciliation")

    # -------------------------------------------------------------
    # Construct complete P3 guide matrix
    # -------------------------------------------------------------

    original = frozen[
        KEY + [
            "guide_called",
            "guide_reads",
            "binomial_p",
            "binomial_q",
            "signature_5_reads",
            "signature_3_reads",
            "total_cell_guide_reads",
            "background_fraction",
        ]
    ].rename(
        columns={
            c: f"stage05k_{c}"
            for c in frozen.columns
            if c in [
                "guide_called",
                "guide_reads",
                "binomial_p",
                "binomial_q",
                "signature_5_reads",
                "signature_3_reads",
                "total_cell_guide_reads",
                "background_fraction",
            ]
        }
    )

    matrix = p3.merge(
        original,
        on=KEY,
        how="left",
        validate="one_to_one",
    )

    matrix["stage05k_frozen_called"] = (
        matrix.stage05k_guide_called.eq(True)
    )

    matrix["stage05l_p3_called"] = (
        matrix.guide_called.astype(bool)
    )

    matrix["assignment_status"] = np.select(
        [
            (
                matrix.stage05k_frozen_called
                & matrix.stage05l_p3_called
            ),
            (
                ~matrix.stage05k_frozen_called
                & matrix.stage05l_p3_called
            ),
            (
                matrix.stage05k_frozen_called
                & ~matrix.stage05l_p3_called
            ),
        ],
        [
            "CONSENSUS_CALL",
            "RESCUED_CALL",
            "FROZEN_ONLY",
        ],
        default="NOT_CALLED",
    )

    matrix["primary_called"] = (
        matrix.stage05l_p3_called
    )

    matrix["primary_call_policy"] = "P3_SPECIFIC"

    # -------------------------------------------------------------
    # Merge read-level QC provenance
    # -------------------------------------------------------------

    require(
        not qc.duplicated(["policy"] + KEY).any(),
        "Duplicate Step05 QC keys",
    )

    qc = qc.loc[
        qc.policy.eq("P3_SPECIFIC")
    ].copy()

    qc_columns = [
        "exact_reads",
        "exact_5_reads",
        "exact_3_reads",
        "d1_reads",
        "d1_5_reads",
        "d1_3_reads",
        "d2_reads",
        "d2_5_reads",
        "d2_3_reads",
        "rescue_reads",
        "rescue_fraction",
        "d2_fraction",
        "exact_both_ends",
        "combined_both_ends",
        "newly_completed_both_ends",
        "is_unique_top",
        "guide_fraction_in_cell",
        "cell_n_called",
        "cell_multiple_called",
        "qc_category",
    ]

    require(
        set(qc_columns).issubset(qc.columns),
        (
            "Missing QC columns: "
            f"{sorted(set(qc_columns) - set(qc.columns))}"
        ),
    )

    matrix = matrix.merge(
        qc[KEY + qc_columns],
        on=KEY,
        how="left",
        validate="one_to_one",
        indicator="_qc_join",
    )

    # -------------------------------------------------------------
    # Validate read provenance
    # -------------------------------------------------------------

    eligible = matrix.assignment_status.ne("NOT_CALLED")

    require(
        matrix.loc[
            eligible, "_qc_join"
        ].eq("both").all(),
        "Some union assignments have no P3 quality evidence",
    )

    for col in [
        "exact_reads",
        "d1_reads",
        "d2_reads",
    ]:

        require(
            matrix.loc[
                eligible, col
            ].notna().all(),
            f"Missing {col} in union",
        )

    require(
        (
            matrix.loc[
                eligible,
                ["exact_reads", "d1_reads", "d2_reads"],
            ].sum(axis=1)
            == matrix.loc[eligible, "guide_reads"]
        ).all(),
        "Read provenance mismatch",
    )

    require(
        (
            matrix.loc[
                eligible,
                "combined_both_ends",
            ].astype(bool)
            == matrix.loc[
                eligible,
                "both_signatures",
            ].astype(bool)
        ).all(),
        "Guide-end evidence mismatch",
    )

    matrix = matrix.drop(columns=["_qc_join"])

    print("PASS: Read provenance reconciliation")

    # -------------------------------------------------------------
    # Assignment counts
    # -------------------------------------------------------------

    counts = matrix.assignment_status.value_counts()

    observed_counts = {
        "frozen": int(
            matrix.stage05k_frozen_called.sum()
        ),
        "p3": int(
            matrix.stage05l_p3_called.sum()
        ),
        "consensus": int(
            counts.get("CONSENSUS_CALL", 0)
        ),
        "rescued": int(
            counts.get("RESCUED_CALL", 0)
        ),
        "frozen_only": int(
            counts.get("FROZEN_ONLY", 0)
        ),
        "union": int(eligible.sum()),
    }

    # General accounting invariants.

    require(
        observed_counts["consensus"]
        + observed_counts["rescued"]
        == observed_counts["p3"],
        "P3 call accounting mismatch",
    )

    require(
        observed_counts["consensus"]
        + observed_counts["frozen_only"]
        == observed_counts["frozen"],
        "Frozen call accounting mismatch",
    )

    require(
        observed_counts["consensus"]
        + observed_counts["rescued"]
        + observed_counts["frozen_only"]
        == observed_counts["union"],
        "Union call accounting mismatch",
    )

    require(
        observed_counts["frozen"]
        == int(frozen.guide_called.sum()),
        "Frozen calls missing from P3 matrix",
    )

    require(
        observed_counts["p3"]
        == int(p3.guide_called.sum()),
        "P3 call mismatch",
    )

    # Pilot-specific category counts.

    for name in [
        "frozen",
        "p3",
        "consensus",
        "rescued",
        "frozen_only",
        "union",
    ]:

        check_pilot_count(
            name,
            observed_counts[name],
            EXPECTED[name],
            args.pilot_regression,
        )

    # -------------------------------------------------------------
    # P3-called cells
    # -------------------------------------------------------------

    p3_cells = matrix.loc[
        matrix.primary_called,
        CELL,
    ].drop_duplicates()

    observed_counts["p3_cells"] = len(p3_cells)

    check_pilot_count(
        "p3_cells",
        len(p3_cells),
        EXPECTED["p3_cells"],
        args.pilot_regression,
    )

    # -------------------------------------------------------------
    # Step06 identity validation
    # -------------------------------------------------------------

    check = matrix[
        KEY + [
            "stage05k_frozen_called",
            "stage05l_p3_called",
        ]
    ].merge(
        identity[
            KEY + [
                "P0_EXACT",
                "P3_SPECIFIC",
            ]
        ],
        on=KEY,
        how="right",
        validate="one_to_one",
        indicator=True,
    )

    require(
        check._merge.eq("both").all(),
        "Step06 identity key absent in P3 matrix",
    )

    require(
        (
            check.stage05k_frozen_called
            == check.P0_EXACT
        ).all(),
        "Frozen identity disagrees with Step06",
    )

    require(
        (
            check.stage05l_p3_called
            == check.P3_SPECIFIC
        ).all(),
        "P3 identity disagrees with Step06",
    )

    print("PASS: Step06 guide identity reconciliation")

    # -------------------------------------------------------------
    # Final assignment table
    # -------------------------------------------------------------

    matrix = matrix.sort_values(
        KEY
    ).reset_index(drop=True)

    assignments = matrix.loc[
        matrix.assignment_status.ne("NOT_CALLED")
    ].copy()

    validate_keys(
        assignments,
        "final assignments",
    )

    require(
        len(assignments) == observed_counts["union"],
        "Final assignment count mismatch",
    )

    # -------------------------------------------------------------
    # Fixed cell universe
    # -------------------------------------------------------------

    cells = fixed[CELL].copy()

    for policy, flag in [
        (
            "P0_EXACT",
            "stage05k_frozen_called",
        ),
        (
            "P3_SPECIFIC",
            "stage05l_p3_called",
        ),
    ]:

        sums = (
            matrix.loc[matrix[flag]]
            .groupby(CELL)
            .size()
            .rename(f"{policy}_n_called")
        )

        cells = cells.join(
            sums,
            on=CELL,
        )

        cells[f"{policy}_n_called"] = (
            cells[f"{policy}_n_called"]
            .fillna(0)
            .astype("int64")
        )

        cells[f"{policy}_status"] = status(
            cells[f"{policy}_n_called"]
        )

        require(
            (
                cells[f"{policy}_n_called"]
                == fixed.set_index(CELL).loc[
                    pd.MultiIndex.from_frame(cells[CELL]),
                    f"{policy}_n_called",
                ].to_numpy()
            ).all(),
            f"Fixed-universe {policy} counts disagree",
        )

    cells["primary_call_policy"] = "P3_SPECIFIC"

    cells["primary_n_called"] = (
        cells.P3_SPECIFIC_n_called
    )

    cells["primary_status"] = (
        cells.P3_SPECIFIC_status
    )

    cells["n_rescued_calls"] = cells.join(
        assignments.loc[
            assignments.assignment_status.eq(
                "RESCUED_CALL"
            )
        ]
        .groupby(CELL)
        .size()
        .rename("_n"),
        on=CELL,
    )["_n"].fillna(0).astype("int64")

    cells["n_frozen_only"] = cells.join(
        assignments.loc[
            assignments.assignment_status.eq(
                "FROZEN_ONLY"
            )
        ]
        .groupby(CELL)
        .size()
        .rename("_n"),
        on=CELL,
    )["_n"].fillna(0).astype("int64")

    require(
        len(cells) == len(fixed),
        "Final cell universe mismatch",
    )

    observed_counts["cell_universe"] = len(cells)

    observed_counts["p3_matrix"] = len(matrix)

    check_pilot_count(
        "cell_universe",
        len(cells),
        EXPECTED["cell_universe"],
        args.pilot_regression,
    )

    check_pilot_count(
        "p3_matrix",
        len(matrix),
        EXPECTED["p3_matrix"],
        args.pilot_regression,
    )

    print("PASS: Fixed cell universe reconciliation")

    # -------------------------------------------------------------
    # Summary tables
    # -------------------------------------------------------------

    summary = (
        assignments.groupby(
            [
                "assignment_status",
                "well_id",
                "guide",
            ]
        )
        .agg(
            assignments=("cell_id", "size"),
            both_ends=("combined_both_ends", "sum"),
            median_guide_reads=("guide_reads", "median"),
            median_p3_q=("binomial_q", "median"),
        )
        .reset_index()
    )

    summary["primary_call_policy"] = "P3_SPECIFIC"

    quality = (
        assignments.groupby(
            [
                "assignment_status",
                "qc_category",
            ]
        )
        .agg(
            assignments=("cell_id", "size"),
            with_exact_reads=(
                "exact_reads",
                lambda s: int((s > 0).sum()),
            ),
            with_d1=(
                "d1_reads",
                lambda s: int((s > 0).sum()),
            ),
            with_d2=(
                "d2_reads",
                lambda s: int((s > 0).sum()),
            ),
            both_ends=("combined_both_ends", "sum"),
            unique_top=("is_unique_top", "sum"),
            median_reads=("guide_reads", "median"),
            median_p3_q=("binomial_q", "median"),
        )
        .reset_index()
    )

    quality["fraction_within_status"] = (
        quality.assignments
        / quality.groupby(
            "assignment_status"
        ).assignments.transform("sum")
    )

    # -------------------------------------------------------------
    # Write outputs
    # -------------------------------------------------------------

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    outputs = {}

    for frame, filename in [
        (
            assignments,
            "final_cell_guide_assignments.parquet",
        ),
        (
            assignments,
            "final_cell_guide_assignments.tsv",
        ),
        (
            matrix,
            "final_cell_guide_matrix.parquet",
        ),
        (
            cells.sort_values(CELL),
            "final_cell_guide_status.tsv",
        ),
        (
            summary,
            "final_guide_assignment_summary.tsv",
        ),
        (
            quality,
            "final_guide_assignment_qc.tsv",
        ),
    ]:

        outputs[filename] = write(
            frame,
            filename,
            outdir,
        )

    # -------------------------------------------------------------
    # Verify input immutability
    # -------------------------------------------------------------

    for name, path in paths.items():

        require(
            digest(path) == input_hashes[name],
            f"Input changed during export: {name}",
        )

    print("PASS: Input SHA256 immutability checks")

    # -------------------------------------------------------------
    # Manifest
    # -------------------------------------------------------------

    manifest = {
        "stage": "stage05l_step07_final_guide_assignment_export",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "primary_policy": "P3_SPECIFIC",
        "primary_call_rule": (
            "within-well BH-adjusted binomial_q < 0.05"
        ),
        "rescue_rule": (
            "eligible spatially assigned distance-1/2 reads; "
            "unique guide within distance 2"
        ),
        "frozen_policy": (
            "P0_EXACT (Stage05k; unchanged)"
        ),
        "assignment_status_definitions": {
            "CONSENSUS_CALL": (
                "called in both frozen Stage05k and P3"
            ),
            "RESCUED_CALL": (
                "called in P3 but not frozen Stage05k"
            ),
            "FROZEN_ONLY": (
                "called in frozen Stage05k but not P3; "
                "not a primary call"
            ),
            "NOT_CALLED": (
                "neither frozen Stage05k nor P3 called"
            ),
        },
        "qc_category_rule": (
            "Descriptive Step05 categories; "
            "not an additional calling filter"
        ),
        "qc_strong_definition": (
            "P3 q<0.05, both guide ends, >=3 guide reads, "
            "unique top construct, >=0.5 fraction of cell "
            "guide reads, exactly one P3 call"
        ),
        "expected_counts": (
            EXPECTED if args.pilot_regression else None
        ),
        "observed_counts": observed_counts,
        "inputs": {
            name: {
                "path": str(path),
                "sha256": input_hashes[name],
            }
            for name, path in paths.items()
        },
        "outputs": {
            name: {
                "path": str(path),
                "rows": int(
                    len(
                        pd.read_parquet(path)
                        if path.suffix == ".parquet"
                        else pd.read_csv(path, sep="\t")
                    )
                ),
                "sha256": digest(path),
            }
            for name, path in outputs.items()
        },
        "notes": [
            (
                "P3 q-values and guide counts are never "
                "replaced with frozen values."
            ),
            (
                "Frozen statistics are stored in separate "
                "stage05k_* columns."
            ),
            (
                "FROZEN_ONLY rows are retained for provenance "
                "and excluded from primary P3 calls."
            ),
            (
                "The full P3 matrix includes NOT_CALLED rows; "
                "assignments export is the union."
            ),
            (
                "No existing Stage05k or Stage05l input "
                "is modified."
            ),
        ],
    }

    manifest_path = (
        outdir / "final_export_manifest.json"
    )

    manifest_path.write_text(
        json.dumps(
            manifest,
            indent=2,
            allow_nan=False,
        ) + "\n"
    )

    print(
        f"{manifest_path.name:48s} written",
        flush=True,
    )

    # -------------------------------------------------------------
    # Final report
    # -------------------------------------------------------------

    print("\n" + "=" * 80)

    print(
        "PASS: exact frozen reconciliation, P3 significance, "
        "identity overlap, read provenance, fixed cell universe, "
        "immutable input hashes"
    )

    print("\nFinal categories")

    for name, value in observed_counts.items():

        print(f"  {name:20s} {value:>10,d}")

    print(f"\nOutputs: {outdir}")

    print("=" * 80)


if __name__ == "__main__":
    main()
