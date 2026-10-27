#!/usr/bin/env python3
"""Recalculate Stage05k guide enrichment after eligible Stage05l mismatch rescue.

Run from bases2fastq-nf/scripts. Never modifies Stage05k inputs.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import binom

EXACT = None
FROZEN = None
RESCUE = None
SPEC = None
OUT = None
GUIDES = ['PTK2B_A', 'PTK2B_B', 'PTK2B_C', 'PTK2B_D', 'PLCG2_A', 'PLCG2_B', 'PLCG2_C']
WELLS = ['A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'D1', 'D2', 'E1', 'E2', 'F1', 'F2']
SIGNATURES = {f'{g}_{end}': g for g in GUIDES for end in ('5', '3')}
POLICIES = [('P0_EXACT', 0, False), ('P1_EXACT_D1', 1, False),
            ('P2_EXACT_D1_D2', 2, False), ('P3_SPECIFIC', 2, True)]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def bh_adjust(values):
    values = np.asarray(values, dtype=float)
    if not len(values):
        return values
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.minimum.accumulate((ranked * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    result = np.empty(len(values), dtype=float)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def load_exact():
    cols = ['read_id', 'target_name', 'label_Cell_global', 'well_id', 'is_nuclear']
    raw = pq.read_table(EXACT, columns=cols).to_pandas(ignore_metadata=True)
    raw = raw.loc[raw['target_name'].isin(SIGNATURES)].copy()
    raw['guide'] = raw['target_name'].map(SIGNATURES)
    raw['guide_end'] = raw['target_name'].str[-1]
    raw['cell_id'] = raw['label_Cell_global'].fillna(0).astype('int64')
    require(len(raw) == 364323, f'P0 exact read count {len(raw)} != 364323')
    require(int((raw.cell_id > 0).sum()) == 194717, 'P0 cellular count mismatch')
    require(int((raw.cell_id <= 0).sum()) == 169606, 'P0 extracellular count mismatch')
    require(raw.loc[raw.cell_id > 0, 'cell_id'].nunique() == 14750, 'P0 cell count mismatch')
    require(not raw[['well_id', 'read_id']].duplicated().any(),
            'Exact reads contain repeated well/read IDs; inspect before proceeding')
    return raw[['well_id', 'read_id', 'guide', 'guide_end', 'cell_id', 'is_nuclear']].assign(mismatch_count=0)


def load_rescue(exact):
    cols = ['well_id', 'read_id', 'guide', 'construct', 'guide_end', 'mismatch_count',
            'cell_id', 'is_nuclear', 'primary_rescue', 'rescue_status', 'spatial_match_count']
    rescue = pq.read_table(RESCUE, columns=cols).to_pandas(ignore_metadata=True)
    rescue = rescue.loc[rescue.mismatch_count.isin([1, 2]) & rescue.primary_rescue.fillna(False)].copy()
    require(len(rescue) == 232386, f'Eligible rescue count {len(rescue)} != 232386')
    require((rescue.rescue_status == 'ELIGIBLE_UNALIGNED').all(), 'Unexpected rescue status')
    require((rescue.spatial_match_count == 1).all(), 'Unexpected spatial multiplicity')
    require(not rescue[['well_id', 'read_id']].duplicated().any(), 'Duplicate rescued read IDs')
    require(not rescue[['well_id', 'read_id']].merge(exact[['well_id', 'read_id']],
                 on=['well_id', 'read_id'], how='inner').shape[0], 'Rescue overlaps exact input')
    require((rescue.guide == rescue.construct + '_' + rescue.guide_end.astype(str)).all(),
            'Guide identity/end inconsistent with construct')
    spec = pq.read_table(SPEC, columns=['well_id', 'read_id', 'guide', 'mismatch_count',
                                         'unique_within_2', 'distance_margin']).to_pandas(ignore_metadata=True)
    require(not spec[['well_id', 'read_id']].duplicated().any(), 'Duplicate specificity IDs')
    rescue = rescue.merge(spec, on=['well_id', 'read_id'], how='left',
                          validate='one_to_one', suffixes=('', '_spec'), indicator=True)
    require((rescue._merge == 'both').all(), 'Missing cross-guide specificity results')
    require((rescue.guide == rescue.guide_spec).all(), 'Cross-guide identity mismatch')
    require((rescue.mismatch_count == rescue.mismatch_count_spec).all(), 'Mismatch distance disagreement')
    rescue['cell_id'] = rescue.cell_id.fillna(0).astype('int64')
    rescue['guide_end'] = rescue.guide_end.astype(str)
    rescue['guide'] = rescue['construct']
    return rescue


def make_matrix(reads):
    require(reads.guide.isin(GUIDES).all(), 'Unexpected construct')
    cellular = reads.loc[reads.cell_id > 0]
    extracellular = reads.loc[reads.cell_id <= 0]
    bg_count = extracellular.groupby(['well_id', 'guide']).size().rename('background_reads')
    bg_total = extracellular.groupby('well_id').size().rename('background_total_reads')
    bg = pd.MultiIndex.from_product([WELLS, GUIDES], names=['well_id', 'guide']).to_frame(index=False)
    bg = bg.join(bg_count, on=['well_id', 'guide']).join(bg_total, on='well_id')
    bg[['background_reads', 'background_total_reads']] = bg[['background_reads', 'background_total_reads']].fillna(0).astype('int64')
    bg['background_fraction'] = np.where(bg.background_total_reads > 0,
                                          bg.background_reads / bg.background_total_reads, np.nan)
    cells = cellular[['well_id', 'cell_id']].drop_duplicates()
    require(not cells.empty, 'No guide-positive cells')
    matrix = cells.merge(pd.DataFrame({'guide': GUIDES}), how='cross')
    observed = (cellular.groupby(['well_id', 'cell_id', 'guide'])
                .agg(guide_reads=('read_id', 'size'), nuclear_reads=('is_nuclear', 'sum'),
                     signature_5_reads=('guide_end', lambda x: int((x == '5').sum())),
                     signature_3_reads=('guide_end', lambda x: int((x == '3').sum())))
                .reset_index())
    matrix = matrix.merge(observed, on=['well_id', 'cell_id', 'guide'], how='left', validate='one_to_one')
    count_cols = ['guide_reads', 'nuclear_reads', 'signature_5_reads', 'signature_3_reads']
    matrix[count_cols] = matrix[count_cols].fillna(0).astype('int64')
    totals = cellular.groupby(['well_id', 'cell_id']).size().rename('total_cell_guide_reads')
    matrix = matrix.join(totals, on=['well_id', 'cell_id'])
    matrix = matrix.merge(bg, on=['well_id', 'guide'], how='left', validate='many_to_one')
    matrix['observed'] = matrix.guide_reads > 0
    matrix['observed_fraction'] = matrix.guide_reads / matrix.total_cell_guide_reads
    matrix['expected_background_reads'] = matrix.total_cell_guide_reads * matrix.background_fraction
    matrix['fold_enrichment'] = np.where(matrix.background_fraction > 0,
                                          matrix.observed_fraction / matrix.background_fraction, np.nan)
    matrix['both_signatures'] = (matrix.signature_5_reads > 0) & (matrix.signature_3_reads > 0)
    # scipy.stats.binom.sf(k-1,n,p) is the exact one-sided binomial tail used by binomtest(..., alternative='greater').
    matrix['binomial_p'] = np.ones(len(matrix), dtype=float)
    missing = matrix.background_fraction.isna()
    positive_zero_bg = (matrix.guide_reads > 0) & (matrix.background_fraction == 0)
    valid_positive = (matrix.guide_reads > 0) & (matrix.background_fraction > 0)
    matrix.loc[missing | positive_zero_bg, 'binomial_p'] = np.nan
    matrix.loc[valid_positive, 'binomial_p'] = binom.sf(
        matrix.loc[valid_positive, 'guide_reads'].to_numpy() - 1,
        matrix.loc[valid_positive, 'total_cell_guide_reads'].to_numpy(),
        matrix.loc[valid_positive, 'background_fraction'].to_numpy())
    matrix['binomial_q'] = np.nan
    for _, indices in matrix.groupby('well_id', sort=False).groups.items():
        indices = matrix.loc[indices].index[matrix.loc[indices, 'binomial_p'].notna()]
        matrix.loc[indices, 'binomial_q'] = bh_adjust(matrix.loc[indices, 'binomial_p'].to_numpy())
    matrix['guide_called'] = matrix.binomial_q.notna() & (matrix.binomial_q < 0.05)
    require((matrix.groupby(['well_id', 'cell_id']).size() == 7).all(), 'Incomplete seven-construct matrix')
    require(int(matrix.guide_reads.sum()) == len(cellular), 'Cellular read accounting mismatch')
    require(int(bg.drop_duplicates('well_id').background_total_reads.sum()) == len(extracellular),
            'Extracellular read accounting mismatch')
    return matrix, bg


def main():
    for path in [EXACT, FROZEN, RESCUE, SPEC]:
        require(path.exists(), f'Missing input: {path}')
    OUT.mkdir(parents=True, exist_ok=True)
    exact = load_exact()
    rescue = load_rescue(exact)
    frozen = pd.read_csv(FROZEN, sep='\t')
    require(len(frozen) == 103250, 'Unexpected frozen evidence size')
    require(int(frozen.guide_called.sum()) == 7503, 'Unexpected frozen call count')
    frozen_keys = frozen[['well_id', 'cell_id', 'guide', 'guide_reads', 'binomial_p',
                          'binomial_q', 'guide_called']]
    summary = []
    comparison = []
    backgrounds = []
    for policy, max_distance, specific in POLICIES:
        print(f'\n{policy}', flush=True)
        if max_distance == 0:
            reads = exact
        else:
            add = rescue.loc[rescue.mismatch_count <= max_distance]
            if specific:
                add = add.loc[add.unique_within_2.fillna(False)]
            add = add[['well_id', 'read_id', 'guide', 'guide_end', 'cell_id', 'is_nuclear', 'mismatch_count']]
            reads = pd.concat([exact, add], ignore_index=True)
        matrix, bg = make_matrix(reads)
        matrix.insert(0, 'policy', policy)
        bg.insert(0, 'policy', policy)
        backgrounds.append(bg)
        if max_distance == 0:
            check = matrix.merge(frozen_keys, on=['well_id', 'cell_id', 'guide'],
                                 how='outer', validate='one_to_one', indicator=True,
                                 suffixes=('', '_frozen'))
            require((check._merge == 'both').all(), 'P0 cell/construct matrix differs from frozen')
            require((check.guide_reads == check.guide_reads_frozen).all(), 'P0 read counts differ')
            require((check.guide_called == check.guide_called_frozen).all(), 'P0 individual calls differ')
            for col in ['binomial_p', 'binomial_q']:
                a = check[col].to_numpy(dtype=float)
                b = check[f'{col}_frozen'].to_numpy(dtype=float)
                require(np.allclose(a, b, rtol=1e-8, atol=1e-12, equal_nan=True),
                        f'P0 {col} values differ from frozen')
            print('PASS: P0 matches frozen counts, p-values, q-values and every individual call')
        merged = matrix[['well_id', 'cell_id', 'guide', 'guide_reads', 'binomial_q',
                         'guide_called']].merge(
            frozen_keys[['well_id', 'cell_id', 'guide', 'guide_called']],
            on=['well_id', 'cell_id', 'guide'], how='left', validate='one_to_one',
            suffixes=('', '_frozen'))
        merged['guide_called_frozen'] = merged.guide_called_frozen.fillna(False).astype(bool)
        merged['call_change'] = np.select(
            [merged.guide_called & ~merged.guide_called_frozen,
             ~merged.guide_called & merged.guide_called_frozen,
             merged.guide_called & merged.guide_called_frozen],
            ['NEW', 'LOST', 'RETAINED'], default='NOT_CALLED')
        merged.insert(0, 'policy', policy)
        comparison.append(merged)
        n_cells = matrix[['well_id', 'cell_id']].drop_duplicates().shape[0]
        called_cells = matrix.loc[matrix.guide_called, ['well_id', 'cell_id']].drop_duplicates().shape[0]
        stats = dict(policy=policy, total_reads=len(reads), cellular_reads=int((reads.cell_id > 0).sum()),
                     extracellular_reads=int((reads.cell_id <= 0).sum()), guide_positive_cells=n_cells,
                     tested_pairs=len(matrix), called_pairs=int(matrix.guide_called.sum()),
                     called_cells=called_cells, new_calls=int((merged.call_change == 'NEW').sum()),
                     lost_calls=int((merged.call_change == 'LOST').sum()),
                     retained_calls=int((merged.call_change == 'RETAINED').sum()))
        summary.append(stats)
        print(stats, flush=True)
        matrix.to_parquet(OUT / f'{policy}_cell_construct_evidence.parquet', index=False)
    combined = pd.concat(comparison, ignore_index=True)
    combined.to_csv(OUT / 'policy_call_comparison.tsv', sep='\t', index=False)
    combined.loc[combined.call_change == 'NEW'].to_csv(
        OUT / 'newly_called_cell_constructs.tsv', sep='\t', index=False)
    combined.loc[combined.call_change == 'LOST'].to_csv(
        OUT / 'lost_cell_construct_calls.tsv', sep='\t', index=False)
    pd.concat(backgrounds, ignore_index=True).to_csv(
        OUT / 'policy_background_fractions.tsv', sep='\t', index=False)
    pd.DataFrame(summary).to_csv(OUT / 'stage05l_calling_summary.tsv', sep='\t', index=False)
    print(f'\nComplete. Outputs: {OUT}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Background-adjusted DISS guide calling'
    )
    parser.add_argument('--exact', type=Path, required=True)
    parser.add_argument('--frozen', type=Path, required=True)
    parser.add_argument('--rescue', type=Path, required=True)
    parser.add_argument('--specificity', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)

    args = parser.parse_args()

    EXACT = args.exact.resolve()
    FROZEN = args.frozen.resolve()
    RESCUE = args.rescue.resolve()
    SPEC = args.specificity.resolve()
    OUT = args.outdir.resolve()

    main()
