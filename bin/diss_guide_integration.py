#!/usr/bin/env python3

"""Stage16h: lossless Stage05i DISS + Stage05l P3 guide integration.

Refactored from:
    20261015_08_stage05l_diss_guide_integration.py

Inputs are read-only.

Cell identifiers are joined on (well, cell_id).
Unassociated transcripts are never guide annotated.

P3 calls are primary.
FROZEN_ONLY is provenance, not a primary assignment.

The --pilot-regression flag enforces the original Teton Atlas pilot counts.
All structural, biological-policy and reconciliation checks remain active
without pilot regression.
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

KEY = ["well", "cell_id"]

EXPECTED = {
    "transcripts": 11522,
    "diss_cells": 3772,
    "guide_cells": 16133,
    "assignments": 11768,
    "p3_calls": 11592,
    "consensus": 7327,
    "rescued": 4265,
    "frozen_only": 176,
}

ASSIGNMENT_STATUSES = {
    "CONSENSUS_CALL",
    "RESCUED_CALL",
    "FROZEN_ONLY",
}


# ---------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Lossless DISS transcript/cell integration with "
            "Stage05l P3 guide assignments"
        )
    )

    parser.add_argument(
        "--diss-transcripts",
        type=Path,
        required=True,
        help="Stage05i DISS_transcripts_final.parquet",
    )

    parser.add_argument(
        "--diss-cell-metadata",
        type=Path,
        required=True,
        help="Stage05i DISS_cell_metadata.tsv.gz",
    )

    parser.add_argument(
        "--guide-assignments",
        type=Path,
        required=True,
        help="Stage05l final_cell_guide_assignments.parquet",
    )

    parser.add_argument(
        "--guide-cell-status",
        type=Path,
        required=True,
        help="Stage05l final_cell_guide_status.tsv",
    )

    parser.add_argument(
        "--guide-export-manifest",
        type=Path,
        required=True,
        help="Stage05l final_export_manifest.json",
    )

    parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Directory for integration outputs",
    )

    parser.add_argument(
        "--pilot-regression",
        action="store_true",
        help="Enforce original Teton Atlas pilot counts",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------
# Utilities
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


def unique_keys(frame, name):

    require(
        frame[KEY].notna().all().all(),
        f"{name}: null well/cell_id",
    )

    require(
        not frame.duplicated(KEY).any(),
        f"{name}: duplicate well/cell_id",
    )


def save(frame, filename, outdir):

    path = outdir / filename

    if filename.endswith(".parquet"):

        frame.to_parquet(
            path,
            index=False,
        )

    elif filename.endswith(".tsv.gz"):

        frame.to_csv(
            path,
            sep="\t",
            index=False,
            compression="gzip",
            na_rep="",
        )

    else:

        frame.to_csv(
            path,
            sep="\t",
            index=False,
            na_rep="",
        )

    print(
        f"{filename:48s} {len(frame):>9,d} rows",
        flush=True,
    )

    return path


def ordered_json(values):

    return json.dumps(
        list(values),
        ensure_ascii=False,
        separators=(",", ":"),
    )


# ---------------------------------------------------------------------
# Guide assignment summary
# ---------------------------------------------------------------------

def summarise_assignments(assignments):

    """One row per cell; JSON arrays retain all guides without multiplication."""

    rows = []

    for (well, cell_id), block in assignments.groupby(
        KEY,
        sort=True,
    ):

        primary = (
            block.loc[
                block.stage05l_p3_called.astype(bool)
            ]
            .sort_values("guide")
        )

        frozen_only = (
            block.loc[
                block.assignment_status.eq("FROZEN_ONLY")
            ]
            .sort_values("guide")
        )

        consensus = (
            block.loc[
                block.assignment_status.eq("CONSENSUS_CALL")
            ]
            .sort_values("guide")
        )

        rescued = (
            block.loc[
                block.assignment_status.eq("RESCUED_CALL")
            ]
            .sort_values("guide")
        )

        rows.append({
            "well": well,
            "cell_id": cell_id,

            "guide_primary_n_called": len(primary),

            "guide_primary_guides_json": ordered_json(
                primary.guide
            ),

            "guide_primary_qc_categories_json": ordered_json(
                primary.qc_category.fillna("")
            ),

            "guide_primary_assignment_statuses_json": ordered_json(
                primary.assignment_status
            ),

            "guide_primary_q_values_json": ordered_json(
                primary.binomial_q.astype(float)
            ),

            "guide_primary_read_counts_json": ordered_json(
                primary.guide_reads.astype(int)
            ),

            "guide_consensus_n": len(consensus),

            "guide_rescued_n": len(rescued),

            "guide_frozen_only_n": len(frozen_only),

            "guide_frozen_only_guides_json": ordered_json(
                frozen_only.guide
            ),

            "guide_primary_any_strong": bool(
                primary.qc_category.eq(
                    "STRONG_DUAL_END_SINGLE_DOMINANT"
                ).any()
            ),
        })

    result = pd.DataFrame(rows)

    unique_keys(
        result,
        "assignment summary",
    )

    return result


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    args = parse_args()

    outdir = args.outdir.resolve()

    paths = {
        "diss_transcripts": args.diss_transcripts.resolve(),
        "diss_cell_metadata": args.diss_cell_metadata.resolve(),
        "guide_assignments": args.guide_assignments.resolve(),
        "guide_cell_status": args.guide_cell_status.resolve(),
        "guide_export_manifest": args.guide_export_manifest.resolve(),
    }

    print("=" * 80)
    print("Stage16h: DISS/guide integration")
    print("=" * 80)

    for label, path in paths.items():

        require(
            path.is_file(),
            f"Missing {label}: {path}",
        )

        print(f"{label:24s} {path}")

    print(f"\nOutput directory: {outdir}")
    print(f"Pilot regression: {args.pilot_regression}")

    # -----------------------------------------------------------------
    # Record input hashes
    # -----------------------------------------------------------------

    hashes = {
        name: sha256(path)
        for name, path in paths.items()
    }

    # -----------------------------------------------------------------
    # Load inputs
    # -----------------------------------------------------------------

    transcripts = pd.read_parquet(
        paths["diss_transcripts"]
    )

    diss_cells = pd.read_csv(
        paths["diss_cell_metadata"],
        sep="\t",
    )

    assignments = pd.read_parquet(
        paths["guide_assignments"]
    )

    guide_cells = (
        pd.read_csv(
            paths["guide_cell_status"],
            sep="\t",
        )
        .rename(columns={"well_id": "well"})
    )

    # -----------------------------------------------------------------
    # Validate required columns
    # -----------------------------------------------------------------

    for name, frame, needed in [

        (
            "transcripts",
            transcripts,
            KEY + [
                "transcript_id",
                "cell_associated",
                "well_cell_id",
            ],
        ),

        (
            "DISS cells",
            diss_cells,
            KEY + [
                "well_cell_id",
                "total_diss_reads",
            ],
        ),

        (
            "guide cells",
            guide_cells,
            KEY + [
                "primary_n_called",
                "primary_status",
                "P3_SPECIFIC_n_called",
                "n_frozen_only",
                "n_rescued_calls",
            ],
        ),

        (
            "assignments",
            assignments,
            [
                "well_id",
                "cell_id",
                "guide",
                "assignment_status",
                "stage05l_p3_called",
                "qc_category",
                "binomial_q",
                "guide_reads",
            ],
        ),
    ]:

        missing = sorted(
            set(needed) - set(frame.columns)
        )

        require(
            not missing,
            f"{name}: missing {missing}",
        )

    # -----------------------------------------------------------------
    # Observed input counts
    # -----------------------------------------------------------------

    observed_input_counts = {
        "transcripts": len(transcripts),
        "diss_cells": len(diss_cells),
        "guide_cells": len(guide_cells),
        "assignments": len(assignments),

        "p3_calls": int(
            assignments.stage05l_p3_called.sum()
        ),

        "consensus": int(
            assignments.assignment_status.eq(
                "CONSENSUS_CALL"
            ).sum()
        ),

        "rescued": int(
            assignments.assignment_status.eq(
                "RESCUED_CALL"
            ).sum()
        ),

        "frozen_only": int(
            assignments.assignment_status.eq(
                "FROZEN_ONLY"
            ).sum()
        ),
    }

    print("\nObserved input counts:")

    for name, count in observed_input_counts.items():
        print(f"{name:24s} {count:>9,d}")

    # -----------------------------------------------------------------
    # Optional pilot regression
    # -----------------------------------------------------------------

    if args.pilot_regression:

        for name, expected in EXPECTED.items():

            observed = observed_input_counts[name]

            require(
                observed == expected,
                (
                    f"Pilot regression mismatch: {name}; "
                    f"expected {expected}, observed {observed}"
                ),
            )

        print("\nPASS: Original pilot input counts")

    # -----------------------------------------------------------------
    # Structural validation
    # -----------------------------------------------------------------

    unique_keys(
        diss_cells,
        "DISS cells",
    )

    unique_keys(
        guide_cells,
        "guide cells",
    )

    require(
        not transcripts.transcript_id.duplicated().any(),
        "Duplicate DISS transcript_id",
    )

    require(
        transcripts.cell_associated.notna().all(),
        "Null cell_associated flag",
    )

    require(
        (
            transcripts.loc[
                transcripts.cell_associated,
                "cell_id",
            ] > 0
        ).all(),
        "Cell-associated transcript with invalid cell_id",
    )

    require(
        transcripts.loc[
            ~transcripts.cell_associated,
            "well_cell_id",
        ].isna().all(),
        "Unassociated transcript has non-null well_cell_id",
    )

    require(
        (
            transcripts.loc[
                transcripts.cell_associated,
                "well_cell_id",
            ]
            ==
            transcripts.loc[
                transcripts.cell_associated,
                "well",
            ].astype(str)
            + ":"
            + transcripts.loc[
                transcripts.cell_associated,
                "cell_id",
            ].astype(str)
        ).all(),
        "DISS transcript well_cell_id mismatch",
    )

    require(
        (
            diss_cells.well_cell_id
            ==
            diss_cells.well.astype(str)
            + ":"
            + diss_cells.cell_id.astype(str)
        ).all(),
        "DISS cell well_cell_id mismatch",
    )

    require(
        not assignments.duplicated(
            ["well_id", "cell_id", "guide"]
        ).any(),
        "Duplicate guide assignment key",
    )

    # -----------------------------------------------------------------
    # Validate assignment statuses and P3 policy
    # -----------------------------------------------------------------

    require(
        assignments.assignment_status.notna().all(),
        "Null guide assignment status",
    )

    require(
        assignments.assignment_status.isin(
            ASSIGNMENT_STATUSES
        ).all(),
        "Unexpected guide assignment status",
    )

    require(
        (
            observed_input_counts["consensus"]
            + observed_input_counts["rescued"]
            + observed_input_counts["frozen_only"]
        )
        == observed_input_counts["assignments"],
        "Guide assignment category counts do not reconcile",
    )

    require(
        (
            assignments.stage05l_p3_called.astype(bool)
            ==
            assignments.assignment_status.isin([
                "CONSENSUS_CALL",
                "RESCUED_CALL",
            ])
        ).all(),
        "Frozen-only incorrectly treated as primary",
    )

    require(
        observed_input_counts["p3_calls"]
        ==
        (
            observed_input_counts["consensus"]
            + observed_input_counts["rescued"]
        ),
        "P3 primary assignment count mismatch",
    )

    require(
        (
            guide_cells.primary_n_called
            ==
            guide_cells.P3_SPECIFIC_n_called
        ).all(),
        "Guide-cell primary call count differs from P3",
    )

    # -----------------------------------------------------------------
    # Validate transcript/cell membership
    # -----------------------------------------------------------------

    transcript_cells = (
        transcripts.loc[
            transcripts.cell_associated,
            KEY,
        ]
        .drop_duplicates()
    )

    transcript_check = transcript_cells.merge(
        diss_cells[KEY],
        on=KEY,
        how="left",
        validate="one_to_one",
        indicator=True,
    )

    require(
        transcript_check._merge.eq("both").all(),
        "Associated transcript cell missing from DISS metadata",
    )

    require(
        len(transcript_cells) == len(diss_cells),
        "DISS metadata and associated transcript cell sets differ",
    )

    # -----------------------------------------------------------------
    # Summarise guide assignments by cell
    # -----------------------------------------------------------------

    a = assignments.rename(
        columns={"well_id": "well"}
    ).copy()

    summary = summarise_assignments(a)

    require(
        summary.guide_primary_n_called.sum()
        == observed_input_counts["p3_calls"],
        "Aggregated primary calls do not reconcile",
    )

    require(
        summary.guide_frozen_only_n.sum()
        == observed_input_counts["frozen_only"],
        "Aggregated frozen-only calls do not reconcile",
    )

    # -----------------------------------------------------------------
    # Merge assignment summaries into guide-cell universe
    # -----------------------------------------------------------------

    guide_cells = guide_cells.merge(
        summary,
        on=KEY,
        how="left",
        validate="one_to_one",
    )

    require(
        (
            guide_cells.primary_n_called
            ==
            guide_cells.guide_primary_n_called.fillna(0).astype("int64")
        ).all(),
        "Cell-level primary call counts disagree with Step07",
    )

    require(
        (
            guide_cells.n_frozen_only
            ==
            guide_cells.guide_frozen_only_n.fillna(0).astype("int64")
        ).all(),
        "Cell-level frozen-only counts disagree with Step07",
    )

    require(
        (
            guide_cells.n_rescued_calls
            ==
            guide_cells.guide_rescued_n.fillna(0).astype("int64")
        ).all(),
        "Cell-level rescued counts disagree with Step07",
    )

    # Populate empty lists for guide-positive but uncalled cells.

    json_cols = [
        c
        for c in guide_cells.columns
        if c.endswith("_json")
    ]

    for col in json_cols:
        guide_cells[col] = guide_cells[col].fillna("[]")

    for col in [
        "guide_primary_n_called",
        "guide_consensus_n",
        "guide_rescued_n",
        "guide_frozen_only_n",
    ]:

        guide_cells[col] = (
            guide_cells[col]
            .fillna(0)
            .astype("int64")
        )

    guide_cells["guide_primary_any_strong"] = (
        guide_cells.guide_primary_any_strong
        .fillna(False)
        .astype(bool)
    )

    guide_cells["guide_evidence_present"] = True

    # -----------------------------------------------------------------
    # Annotate all DISS cells
    # -----------------------------------------------------------------

    diss_annotated = diss_cells.merge(
        guide_cells.drop(
            columns=["well_cell_id"],
            errors="ignore",
        ),
        on=KEY,
        how="left",
        validate="one_to_one",
        indicator="_guide_join",
    )

    diss_annotated["guide_evidence_present"] = (
        diss_annotated._guide_join.eq("both")
    )

    diss_annotated["guide_annotation_state"] = np.select(
        [
            ~diss_annotated.guide_evidence_present,

            diss_annotated.primary_n_called.fillna(0).gt(0),
        ],
        [
            "NO_GUIDE_EVIDENCE",
            "P3_CALLED",
        ],
        default="GUIDE_EVIDENCE_UNCALLED",
    )

    diss_annotated = diss_annotated.drop(
        columns=["_guide_join"]
    )

    require(
        len(diss_annotated) == len(diss_cells),
        "DISS cell rows multiplied",
    )

    # -----------------------------------------------------------------
    # Annotate full guide-cell universe with DISS membership
    # -----------------------------------------------------------------

    diss_extra = diss_cells.rename(
        columns={
            c: f"diss_{c}"
            for c in diss_cells.columns
            if c not in KEY
        }
    )

    guide_annotated = guide_cells.merge(
        diss_extra,
        on=KEY,
        how="left",
        validate="one_to_one",
        indicator="_diss_join",
    )

    guide_annotated["diss_cell_present"] = (
        guide_annotated._diss_join.eq("both")
    )

    guide_annotated = guide_annotated.drop(
        columns=["_diss_join"]
    )

    require(
        len(guide_annotated) == len(guide_cells),
        "Guide cell rows multiplied",
    )

    # -----------------------------------------------------------------
    # Transcript-level annotation
    # -----------------------------------------------------------------

    annotation_columns = [
        c
        for c in diss_annotated.columns
        if c.startswith("guide_")
        or c in [
            "primary_call_policy",
            "primary_n_called",
            "primary_status",
            "P0_EXACT_n_called",
            "P0_EXACT_status",
            "P3_SPECIFIC_n_called",
            "P3_SPECIFIC_status",
            "n_rescued_calls",
            "n_frozen_only",
        ]
    ]

    annotated = transcripts.merge(
        diss_annotated[KEY + annotation_columns],
        on=KEY,
        how="left",
        validate="many_to_one",
        indicator="_cell_join",
    )

    require(
        len(annotated) == len(transcripts),
        "Transcript rows multiplied",
    )

    require(
        annotated.transcript_id.equals(
            transcripts.transcript_id
        ),
        "Transcript order or identity changed",
    )

    require(
        annotated.loc[
            annotated.cell_associated,
            "_cell_join",
        ].eq("both").all(),
        "Associated transcript failed to join DISS cell metadata",
    )

    require(
        annotated.loc[
            ~annotated.cell_associated,
            "_cell_join",
        ].eq("left_only").all(),
        "Unassociated transcript unexpectedly joined cell metadata",
    )

    annotated = annotated.drop(
        columns=["_cell_join"]
    )

    require(
        annotated.loc[
            ~annotated.cell_associated,
            "guide_annotation_state",
        ].isna().all(),
        "Unassociated transcripts acquired guide annotation",
    )

    # -----------------------------------------------------------------
    # Full guide assignment table with DISS membership
    # -----------------------------------------------------------------

    assignment_diss = a.merge(
        diss_extra[
            KEY + [
                "diss_well_cell_id",
                "diss_total_diss_reads",
                "diss_unique_targets",
            ]
        ],
        on=KEY,
        how="left",
        validate="many_to_one",
    )

    assignment_diss["diss_cell_present"] = (
        assignment_diss.diss_well_cell_id.notna()
    )

    require(
        len(assignment_diss)
        == observed_input_counts["assignments"],
        "Guide assignments multiplied or lost",
    )

    # -----------------------------------------------------------------
    # By-well reconciliation
    # -----------------------------------------------------------------

    by_well = []

    wells = sorted(
        set(transcripts.well)
        | set(guide_cells.well)
    )

    for well in wells:

        d = diss_annotated.loc[
            diss_annotated.well.eq(well)
        ]

        g = guide_annotated.loc[
            guide_annotated.well.eq(well)
        ]

        t = transcripts.loc[
            transcripts.well.eq(well)
        ]

        by_well.append({
            "well": well,

            "diss_transcripts": len(t),

            "cell_associated_transcripts": int(
                t.cell_associated.sum()
            ),

            "diss_cells": len(d),

            "guide_evidence_cells": len(g),

            "diss_cells_with_guide_evidence": int(
                d.guide_evidence_present.sum()
            ),

            "diss_cells_with_p3_calls": int(
                d.guide_annotation_state.eq("P3_CALLED").sum()
            ),

            "diss_cells_with_guide_evidence_uncalled": int(
                d.guide_annotation_state.eq(
                    "GUIDE_EVIDENCE_UNCALLED"
                ).sum()
            ),

            "diss_cells_without_guide_evidence": int(
                d.guide_annotation_state.eq(
                    "NO_GUIDE_EVIDENCE"
                ).sum()
            ),

            "guide_cells_with_diss": int(
                g.diss_cell_present.sum()
            ),

            "p3_assignments_all_cells": int(
                g.primary_n_called.sum()
            ),

            "p3_assignments_diss_cells": int(
                g.loc[
                    g.diss_cell_present,
                    "primary_n_called",
                ].sum()
            ),
        })

    by_well = pd.DataFrame(by_well)

    require(
        by_well.diss_transcripts.sum()
        == observed_input_counts["transcripts"],
        "By-well transcript mismatch",
    )

    require(
        by_well.diss_cells.sum()
        == observed_input_counts["diss_cells"],
        "By-well DISS cell mismatch",
    )

    require(
        by_well.guide_evidence_cells.sum()
        == observed_input_counts["guide_cells"],
        "By-well guide cell mismatch",
    )

    require(
        by_well.p3_assignments_all_cells.sum()
        == observed_input_counts["p3_calls"],
        "By-well P3 mismatch",
    )

    # -----------------------------------------------------------------
    # Save outputs
    # -----------------------------------------------------------------

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    outputs = {}

    for frame, filename in [

        (
            annotated,
            "DISS_transcripts_with_guides.parquet",
        ),

        (
            annotated,
            "DISS_transcripts_with_guides.tsv.gz",
        ),

        (
            diss_annotated,
            "DISS_cells_with_guides.parquet",
        ),

        (
            diss_annotated,
            "DISS_cells_with_guides.tsv.gz",
        ),

        (
            guide_annotated,
            "guide_cells_with_DISS.parquet",
        ),

        (
            assignment_diss,
            "DISS_cell_guide_assignments.parquet",
        ),

        (
            by_well,
            "DISS_guide_integration_by_well.tsv",
        ),
    ]:

        outputs[filename] = save(
            frame,
            filename,
            outdir,
        )

    # -----------------------------------------------------------------
    # Verify inputs have not changed
    # -----------------------------------------------------------------

    for label, path in paths.items():

        require(
            sha256(path) == hashes[label],
            f"Input changed: {label}",
        )

    # -----------------------------------------------------------------
    # Manifest
    # -----------------------------------------------------------------

    observed = {
        "transcripts": len(annotated),

        "diss_cells": len(diss_annotated),

        "guide_cells": len(guide_annotated),

        "assignments": len(assignment_diss),

        "p3_calls": int(
            assignment_diss.stage05l_p3_called.sum()
        ),

        "diss_cells_with_guide_evidence": int(
            diss_annotated.guide_evidence_present.sum()
        ),

        "diss_cells_with_p3_calls": int(
            diss_annotated.guide_annotation_state.eq(
                "P3_CALLED"
            ).sum()
        ),

        "associated_transcripts": int(
            transcripts.cell_associated.sum()
        ),
    }

    manifest = {
        "stage": "stage05l_step08_diss_guide_integration",

        "created_utc": datetime.now(
            timezone.utc
        ).isoformat(),

        "join_key": KEY,

        "primary_policy": "P3_SPECIFIC",

        "primary_call_definition": (
            "Stage05l P3 q<0.05; frozen-only is not primary"
        ),

        "cell_annotation_states": [
            "NO_GUIDE_EVIDENCE",
            "GUIDE_EVIDENCE_UNCALLED",
            "P3_CALLED",
        ],

        "transcript_policy": (
            "Preserve every transcript; annotate associated cells only; "
            "no row multiplication"
        ),

        "multi_guide_policy": (
            "JSON arrays ordered by guide; "
            "full assignments provided separately"
        ),

        "pilot_regression": args.pilot_regression,

        "expected_counts": (
            EXPECTED if args.pilot_regression else None
        ),

        "input_counts": observed_input_counts,

        "observed_counts": observed,

        "inputs": {
            name: {
                "path": str(path),
                "sha256": hashes[name],
            }
            for name, path in paths.items()
        },

        "outputs": {
            name: {
                "path": str(path),
                "sha256": sha256(path),
            }
            for name, path in outputs.items()
        },

        "notes": [
            "No source files are modified.",
            (
                "Guide evidence absence differs from guide evidence "
                "without significant P3 calls."
            ),
            (
                "DISS-negative guide cells and guide-negative "
                "DISS cells are both retained."
            ),
            (
                "No spatial coordinate inference or "
                "nearest-neighbour matching is performed."
            ),
        ],
    }

    manifest_path = (
        outdir / "DISS_guide_integration_manifest.json"
    )

    manifest_path.write_text(
        json.dumps(
            manifest,
            indent=2,
            allow_nan=False,
        ) + "\n"
    )

    print(
        f"{manifest_path.name:48s} written"
    )

    print("\nDISS cell annotation states:")

    print(
        diss_annotated.guide_annotation_state
        .value_counts()
        .to_string()
    )

    print(
        "\nPASS: lossless DISS transcript and cell joins; "
        "full guide universe; P3 call counts; "
        "frozen-only provenance; unchanged input hashes"
    )

    print(f"Outputs: {outdir}", flush=True)


if __name__ == "__main__":
    main()
