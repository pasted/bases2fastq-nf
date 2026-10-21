#!/usr/bin/env python3
"""Stage 15: lossless, validated DISS/guide export bundle.

Run from scripts/, or pass --analysis-root. Requires pandas and pyarrow.
No guide calls or editing statuses are inferred or changed.
"""
import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def check_unique(df, keys, name):
    require(set(keys).issubset(df.columns), f'{name}: missing keys {set(keys) - set(df.columns)}')
    require(not df[keys].isna().any().any(), f'{name}: null keys')
    require(not df.duplicated(keys).any(), f'{name}: duplicate keys {keys}')


def same_values(left, right):
    """Compare aligned series, allowing nulls and minor float serialization differences."""
    if pd.api.types.is_numeric_dtype(left) and pd.api.types.is_numeric_dtype(right):
        return bool(np.isclose(left.to_numpy(dtype=float, na_value=np.nan),
                               right.to_numpy(dtype=float, na_value=np.nan),
                               rtol=1e-10, atol=1e-12, equal_nan=True).all())
    return bool(left.astype('string').fillna('<NULL>').eq(
        right.astype('string').fillna('<NULL>')).all())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis-root', type=Path,
                        default=Path('../results/Teton_Atlas_LO_pilot_B02/analysis'))
    parser.add_argument('--outdir', type=Path, default=None)
    parser.add_argument('--pilot-regression', action='store_true',
                        help='Enforce frozen pilot counts (do not use for future runs).')
    args = parser.parse_args()
    root = args.analysis_root.resolve()
    out = (args.outdir or root / 'stage05l_final_interoperability_exports').resolve()
    require(out != root, 'Output directory must differ from analysis root')
    inputs = {
        'transcripts': root / 'stage05l_diss_guide_integration/DISS_transcripts_with_guides.parquet',
        'cells': root / 'stage05l_diss_guide_integration/DISS_cells_with_guides.parquet',
        'assignments': root / 'stage05l_diss_guide_integration/DISS_cell_guide_assignments.parquet',
        'guide_cells': root / 'stage05l_diss_guide_integration/guide_cells_with_DISS.parquet',
        'evidence': root / 'stage05l_guide_rescue_evidence_audit/guide_rescue_assignment_evidence.parquet',
        'family': root / 'stage05l_guide_specific_cutsite_detection/diss_cells_guide_family_classes.parquet',
        'features': root / 'stage05i_interoperability_exports/DISS_feature_metadata.tsv',
        'cell_metadata': root / 'stage05i_interoperability_exports/DISS_cell_metadata.tsv.gz',
        'target_counts': root / 'stage05i_interoperability_exports/DISS_cell_target_counts.tsv.gz',
        'allele_counts': root / 'stage05i_interoperability_exports/DISS_cell_allele_counts.tsv.gz',
    }
    require(all(p.is_file() for p in inputs.values()),
            'Missing inputs: ' + ', '.join(k for k, p in inputs.items() if not p.is_file()))
    hashes_before = {key: sha256(path) for key, path in inputs.items()}
    t = pd.read_parquet(inputs['transcripts'])
    c = pd.read_parquet(inputs['cells'])
    a = pd.read_parquet(inputs['assignments'])
    g = pd.read_parquet(inputs['guide_cells'])
    e = pd.read_parquet(inputs['evidence'])
    f = pd.read_parquet(inputs['family'])
    feat = pd.read_csv(inputs['features'], sep='\t')
    cm = pd.read_csv(inputs['cell_metadata'], sep='\t')
    checks = []

    def validate(name, condition, detail=''):
        require(condition, f'FAIL {name}: {detail}')
        checks.append({'check': name, 'status': 'PASS', 'detail': str(detail)})

    check_unique(t, ['transcript_id'], 'transcripts')
    check_unique(c, ['well', 'cell_id'], 'cells')
    check_unique(c, ['well_cell_id'], 'cells')
    check_unique(f, ['well', 'cell_id'], 'family')
    check_unique(cm, ['well', 'cell_id'], 'cell_metadata')
    check_unique(a, ['well', 'cell_id', 'guide'], 'assignments')
    check_unique(e, ['well_id', 'cell_id', 'guide'], 'evidence')
    validate('assignment_evidence_row_count', len(a) == len(e), len(a))
    validate('family_cell_row_count', len(f) == len(c), len(f))
    validate('cell_metadata_row_count', len(cm) == len(c), len(cm))
    validate('transcript_cell_association',
             t.loc[t.cell_associated, 'well_cell_id'].notna().all() and
             t.loc[~t.cell_associated, 'well_cell_id'].isna().all(),
             f'associated={int(t.cell_associated.sum())}')
    guide_cols = [col for col in t.columns if col.startswith('guide_') or col in
                  ['P0_EXACT_n_called', 'P0_EXACT_status', 'P3_SPECIFIC_n_called',
                   'P3_SPECIFIC_status', 'primary_call_policy', 'primary_n_called',
                   'primary_status', 'n_rescued_calls', 'n_frozen_only']]
    validate('unassociated_transcripts_not_guide_annotated',
             t.loc[~t.cell_associated, guide_cols].isna().all().all())
    cell_keys = pd.MultiIndex.from_frame(c[['well', 'cell_id']])
    associated_keys = pd.MultiIndex.from_frame(t.loc[t.cell_associated, ['well', 'cell_id']])
    validate('associated_transcripts_have_cells', associated_keys.isin(cell_keys).all())
    validate('cell_family_keys_match',
             set(map(tuple, c[['well', 'cell_id']].itertuples(index=False, name=None))) ==
             set(map(tuple, f[['well', 'cell_id']].itertuples(index=False, name=None))))
    validate('cell_metadata_keys_match',
             set(map(tuple, c[['well', 'cell_id']].itertuples(index=False, name=None))) ==
             set(map(tuple, cm[['well', 'cell_id']].itertuples(index=False, name=None))))

    # Join evidence on three-part key; only append fields absent from assignment table.
    # Explicitly harmonize well_id -> well; never merge on cell_id alone.
    e = e.rename(columns={'well_id': 'well'})
    key = ['well', 'cell_id', 'guide']
    common = [col for col in a.columns if col in e.columns and col not in key]
    paired = a[key + common].merge(e[key + common], on=key, how='outer',
                                   validate='one_to_one', indicator=True,
                                   suffixes=('_assignment', '_evidence'))
    validate('assignment_evidence_keys_match', paired['_merge'].eq('both').all())
    for col in common:
        validate(f'evidence_shared_{col}',
                 same_values(paired[col + '_assignment'], paired[col + '_evidence']))
    extra = [col for col in e.columns if col not in a.columns]
    final_a = a.merge(e[key + extra], on=key, how='left', validate='one_to_one', indicator=True)
    validate('assignment_evidence_join_complete', final_a['_merge'].eq('both').all())
    final_a = final_a.drop(columns='_merge')
    final_a['well_cell_id'] = final_a.well.astype(str) + ':' + final_a.cell_id.astype(str)
    # Independent flags; preserve the exclusive evidence_flag unchanged.
    flags = {
        'flag_no_both_end_support': ~final_a.both_end_support.astype(bool),
        'flag_not_unique_read_top': ~final_a.unique_read_top.astype(bool),
        'flag_mismatch_only_read_support': final_a.mismatch_only_support.astype(bool),
        'flag_q_near_threshold': final_a.near_q_threshold.astype(bool),
    }
    for name, values in flags.items():
        final_a[name] = values
    final_a['n_independent_flags'] = final_a[list(flags)].sum(axis=1)
    final_a['any_independent_flag'] = final_a.n_independent_flags.gt(0)
    priority = np.select([flags[name] for name in flags],
                         ['NO_BOTH_END_SUPPORT', 'NOT_UNIQUE_READ_TOP',
                          'MISMATCH_ONLY_READ_SUPPORT', 'Q_NEAR_THRESHOLD'],
                         default='NO_LISTED_FLAG')
    validate('exclusive_priority_flag_reproduced',
             np.array_equal(priority, final_a.evidence_flag.to_numpy()))
    validate('read_origin_conservation',
             np.isclose(final_a.guide_reads, final_a.exact_reads + final_a.rescue_reads).all())
    validate('p3_and_frozen_status',
             final_a.loc[final_a.assignment_status.eq('FROZEN_ONLY'), 'stage05l_p3_called'].eq(False).all() and
             final_a.loc[~final_a.assignment_status.eq('FROZEN_ONLY'), 'stage05l_p3_called'].eq(True).all())

    # Add only the derived guide-family fields to cell export; preserve original cell fields.
    derived = ['n_p3', 'n_plcg2', 'n_ptk2b', 'guide_list', 'guide_family_class']
    validate('family_derived_columns_available', set(derived).issubset(f.columns))
    final_c = c.merge(f[['well', 'cell_id'] + derived], on=['well', 'cell_id'],
                      how='left', validate='one_to_one', indicator=True)
    validate('cell_family_join_complete', final_c['_merge'].eq('both').all())
    final_c = final_c.drop(columns='_merge')
    validate('feature_ids_unique', not feat.feature_id.duplicated().any())
    validate('transcript_features_defined',
             t.feature_id.dropna().isin(feat.feature_id).all())
    feat['editing_status'] = 'NOT_ASSESSED'
    validate('reference_provenance_present',
             feat[['production_reference_source', 'empirical_reference_override']].notna().all().all())
    if args.pilot_regression:
        expected = {'transcripts': 11522, 'associated': 5640, 'cells': 3772,
                    'guide_cells': 16133, 'assignments': 11768,
                    'CONSENSUS_CALL': 7327, 'RESCUED_CALL': 4265, 'FROZEN_ONLY': 176}
        actual = {'transcripts': len(t), 'associated': int(t.cell_associated.sum()),
                  'cells': len(c), 'guide_cells': len(g), 'assignments': len(a),
                  **a.assignment_status.value_counts().to_dict()}
        validate('pilot_regression_counts', all(actual.get(k) == v for k, v in expected.items()),
                 json.dumps(actual, sort_keys=True))

    out.mkdir(parents=True, exist_ok=True)
    outputs = {}

    def write_parquet(name, frame):
        path = out / name
        frame.to_parquet(path, index=False)
        outputs[name] = path
        print(f'{name:48s} {len(frame):>9,d} rows')

    def write_tsv(name, frame):
        path = out / name
        frame.to_csv(path, sep='\t', index=False)
        outputs[name] = path
        print(f'{name:48s} {len(frame):>9,d} rows')

    write_parquet('DISS_transcripts_final.parquet', t)
    write_parquet('DISS_cells_final.parquet', final_c)
    write_parquet('DISS_guide_assignments_final.parquet', final_a)
    write_parquet('DISS_guide_cells_final.parquet', g)
    write_tsv('DISS_features_final.tsv', feat)
    for src, dest in [('target_counts', 'DISS_cell_target_counts_final.tsv.gz'),
                      ('allele_counts', 'DISS_cell_allele_counts_final.tsv.gz')]:
        path = out / dest
        shutil.copyfile(inputs[src], path)
        validate(f'{src}_byte_identical', sha256(path) == hashes_before[src])
        outputs[dest] = path
        print(f'{dest:48s} copied losslessly')

    # Include a schema of all delivered columns and their observed storage dtypes.
    frames = {'DISS_transcripts_final.parquet': t, 'DISS_cells_final.parquet': final_c,
              'DISS_guide_assignments_final.parquet': final_a,
              'DISS_guide_cells_final.parquet': g, 'DISS_features_final.tsv': feat}
    schema = {
        'schema_version': '1.0.0',
        'join_keys': {'cell': ['well', 'cell_id'],
                      'guide_assignment': ['well', 'cell_id', 'guide'],
                      'transcript': ['transcript_id']},
        'semantics': {'editing_status': 'NOT_ASSESSED; mismatch count is not editing evidence',
                      'evidence_flag': 'mutually exclusive priority category',
                      'flag_*': 'independent overlapping QC indicators',
                      'guide_call_policy': 'P3_SPECIFIC; FROZEN_ONLY retained separately',
                      'unassociated_transcripts': 'guide annotation fields null'},
        'files': {name: {'rows': len(frame),
                         'columns': {col: str(dtype) for col, dtype in frame.dtypes.items()}}
                  for name, frame in frames.items()},
    }
    schema_path = out / 'DISS_export_schema.json'
    schema_path.write_text(json.dumps(schema, indent=2, sort_keys=True) + '\n')
    outputs[schema_path.name] = schema_path
    write_tsv('DISS_export_validation.tsv', pd.DataFrame(checks))
    hashes_after = {key: sha256(path) for key, path in inputs.items()}
    require(hashes_before == hashes_after, 'Upstream input hashes changed')
    manifest = {
        'stage': 'stage05l_final_interoperability_exports',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'schema_version': '1.0.0',
        'pilot_regression': args.pilot_regression,
        'inputs': {key: {'path': str(path), 'sha256': hashes_before[key]}
                   for key, path in inputs.items()},
        'outputs': {name: {'sha256': sha256(path), 'bytes': path.stat().st_size}
                    for name, path in outputs.items()},
        'validation_count': len(checks),
        'counts': {'transcripts': len(t), 'associated_transcripts': int(t.cell_associated.sum()),
                   'unassociated_transcripts': int((~t.cell_associated).sum()),
                   'cells': len(c), 'guide_cells': len(g), 'assignments': len(final_a),
                   'assignment_status': a.assignment_status.value_counts().to_dict()},
        'notes': ['No guide calls changed', 'No editing inference',
                  'Guide-evidence cell universe is distinct from DISS cell universe'],
    }
    manifest_path = out / 'DISS_export_manifest.json'
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    print(f'{manifest_path.name:48s} written')
    print(f'PASS: {len(checks)} checks; upstream hashes unchanged')
    print(f'Outputs: {out}')


if __name__ == '__main__':
    main()
