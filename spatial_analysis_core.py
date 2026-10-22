#!/usr/bin/env python3
"""Reusable, study-agnostic spatial transcriptomics analysis engine.

Input expression table: rows = features, columns = segment IDs (CSV or TSV).
Input metadata table: one row per segment with segment_id, sample_id, group;
optional columns: compartment, region, site, batch, qc_pass.

Example:
  python spatial_analysis_core.py --expression normalized.csv --metadata metadata.csv \
      --out results --contrasts GroupA:GroupB --contrasts GroupB:GroupC

Inference defaults to subject-level Welch tests, so repeated segments from a
single biological sample are never counted as independent replicates.
Expression must already be normalized and log-transformed for meaningful
model effect estimates. For raw GeoMx counts, use the companion adapter.
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

LOG = logging.getLogger('spatial_analysis')


@dataclass(frozen=True)
class Settings:
    alpha: float = 0.05
    min_samples_per_group: int = 3
    max_heatmap_features: int = 50
    max_embedding_features: int = 1500
    seed: int = 42
    image_dpi: int = 600


def sanitize(value: object) -> str:
    return re.sub(r'[^A-Za-z0-9._-]+', '_', str(value)).strip('_')[:110]


def read_table(path: Path, **kwargs) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {'.xlsx', '.xls'}:
        return pd.read_excel(path, **kwargs)
    if suffix in {'.csv', '.tsv', '.txt'}:
        return pd.read_csv(path, sep='\t' if suffix in {'.tsv', '.txt'} else ',', **kwargs)
    raise ValueError(f'Unsupported table extension: {path}')


def validate_inputs(expr: pd.DataFrame, metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {'segment_id', 'sample_id', 'group'}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f'Missing metadata columns: {sorted(missing)}')
    metadata = metadata.copy()
    for column in required:
        metadata[column] = metadata[column].astype('string').str.strip()
    if metadata['segment_id'].isna().any() or metadata['segment_id'].duplicated().any():
        raise ValueError('Metadata segment IDs must be unique and nonmissing')
    if metadata[['sample_id', 'group']].isna().any().any():
        raise ValueError('Sample identity and experimental group must be available for every segment')
    if expr.index.has_duplicates or expr.columns.has_duplicates:
        raise ValueError('Expression feature and segment identifiers must be unique')
    expr.index = expr.index.astype(str).str.strip()
    expr.columns = expr.columns.astype(str).str.strip()
    metadata = metadata.set_index('segment_id', drop=False)
    common = metadata.index.intersection(expr.columns, sort=False)
    if len(common) < 4:
        raise ValueError('Fewer than four aligned expression segments')
    missing_expr = metadata.index.difference(expr.columns)
    if len(missing_expr):
        LOG.warning('%d metadata segments have no expression measurements', len(missing_expr))
    extra_expr = expr.columns.difference(metadata.index)
    if len(extra_expr):
        LOG.warning('%d expression segments have no metadata', len(extra_expr))
    metadata = metadata.loc[common].copy()
    expr = expr.loc[:, common].apply(pd.to_numeric, errors='coerce')
    if not np.isfinite(expr.to_numpy()).all():
        raise ValueError('Expression matrix includes missing or infinite values')
    if 'qc_pass' in metadata:
        keep = metadata.qc_pass.astype(str).str.lower().isin(['true', '1', 'pass', 'yes'])
        LOG.info('Excluding %d QC-failed segments', int((~keep).sum()))
        metadata = metadata.loc[keep].copy()
        expr = expr.loc[:, metadata.index]
    if expr.shape[0] < 2 or expr.shape[1] < 4:
        raise ValueError('Insufficient features or QC-passing segments')
    collisions = metadata.groupby('sample_id').group.nunique()
    if (collisions > 1).any():
        raise ValueError('A biological sample has more than one disease/group label')
    return expr, metadata


def write_table(df: pd.DataFrame, out: Path, name: str, index: bool = False) -> None:
    df.to_csv(out / f'{name}.csv', index=index)
    try:
        df.to_excel(out / f'{name}.xlsx', index=index)
    except ImportError:
        LOG.warning('openpyxl unavailable; XLSX export skipped')


def save_plot(fig: plt.Figure, out: Path, name: str, dpi: int) -> None:
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(out / f'{name}.{ext}', dpi=dpi, bbox_inches='tight')
    plt.close(fig)


def subject_averages(expr: pd.DataFrame, metadata: pd.DataFrame, stratify: tuple[str, ...] = ()) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Collapse repeated AOIs within each independent sample and stratum."""
    fields = ['sample_id', 'group'] + [x for x in stratify if x in metadata.columns]
    design = metadata[fields].copy()
    merged = expr.T.join(design)
    grouped = merged.groupby(fields, dropna=False, observed=True, sort=False)
    averaged = grouped[expr.index.tolist()].mean().reset_index()
    return averaged[expr.index.tolist()], averaged[fields]


def differential_expression(expr: pd.DataFrame, metadata: pd.DataFrame,
                            group_a: str, group_b: str, settings: Settings,
                            stratify: dict[str, str] | None = None) -> pd.DataFrame:
    """Biological-replicate-level differential expression with BH adjustment.

    A Welch comparison uses log2-normalized values averaged per subject.
    The estimate is group A minus group B on the provided expression scale.
    This is exploratory and does not replicate a vendor mixed-effect model.
    """
    stratify = stratify or {}
    selected = metadata.group.isin([group_a, group_b])
    for col, value in stratify.items():
        if col not in metadata:
            raise ValueError(f'Unknown stratum field: {col}')
        selected &= metadata[col].astype(str).eq(str(value))
    if not selected.any():
        return pd.DataFrame()
    aggregated_expr, design = subject_averages(expr.loc[:, selected], metadata.loc[selected])
    xa = aggregated_expr.loc[design.group.eq(group_a)].to_numpy(dtype=float)
    xb = aggregated_expr.loc[design.group.eq(group_b)].to_numpy(dtype=float)
    if min(len(xa), len(xb)) < settings.min_samples_per_group:
        LOG.warning('Skipping contrast %s vs %s, stratum=%s: independent n=%d/%d',
                    group_a, group_b, stratify, len(xa), len(xb))
        return pd.DataFrame()
    with np.errstate(invalid='ignore', divide='ignore'):
        test = stats.ttest_ind(xa, xb, equal_var=False, axis=0, nan_policy='omit')
    effect = xa.mean(axis=0) - xb.mean(axis=0)
    p = np.asarray(test.pvalue, dtype=float)
    good = np.isfinite(p)
    q = np.full(len(p), np.nan)
    if good.any():
        q[good] = multipletests(p[good], method='fdr_bh')[1]
    out = pd.DataFrame({'feature': aggregated_expr.columns, 'estimate_log2': effect,
                        'p_value': p, 'FDR': q,
                        'n_samples_a': len(xa), 'n_samples_b': len(xb),
                        'group_a': group_a, 'group_b': group_b,
                        'method': 'subject_mean_Welch'})
    for col, value in stratify.items():
        out[col] = value
    return out.sort_values('p_value', na_position='last').reset_index(drop=True)


def volcano(df: pd.DataFrame, destination: Path, name: str, settings: Settings) -> None:
    valid = df.dropna(subset=['estimate_log2', 'p_value']).copy()
    if valid.empty:
        return
    fig, ax = plt.subplots(figsize=(7, 5.4))
    hits = valid.FDR.lt(settings.alpha)
    ax.scatter(valid.loc[~hits, 'estimate_log2'],
               -np.log10(valid.loc[~hits, 'p_value'].clip(lower=1e-300)),
               s=13, alpha=0.55, c='#8a94a6', rasterized=True)
    ax.scatter(valid.loc[hits, 'estimate_log2'],
               -np.log10(valid.loc[hits, 'p_value'].clip(lower=1e-300)),
               s=17, alpha=0.85, c='#b83b5e', rasterized=True)
    for _, row in valid.loc[hits].nsmallest(10, 'FDR').iterrows():
        ax.annotate(row['feature'], (row['estimate_log2'], -np.log10(max(row['p_value'], 1e-300))),
                    fontsize=7)
    ax.axhline(-np.log10(settings.alpha), color='#666666', ls='--', lw=.8)
    ax.set(xlabel='Mean difference in log2 expression (A - B)',
           ylabel='-log10(nominal p-value)', title=name.replace('_', ' '))
    save_plot(fig, destination, name, settings.image_dpi)


def exploratory_visualizations(expr: pd.DataFrame, meta: pd.DataFrame,
                               dest: Path, settings: Settings) -> None:
    sns.set_theme(style='ticks', context='paper')
    variance = expr.var(axis=1).sort_values(ascending=False)
    selected = variance.head(settings.max_embedding_features).index
    x = expr.loc[selected].T
    scale = x.std(axis=0).replace(0, 1)
    x = (x - x.mean(axis=0)) / scale
    pca = PCA(n_components=2, random_state=settings.seed)
    coords = pca.fit_transform(x)
    plotdata = meta[['segment_id', 'sample_id', 'group']].copy()
    for field in ('compartment', 'region', 'site', 'batch'):
        if field in meta:
            plotdata[field] = meta[field]
    plotdata['PC1'], plotdata['PC2'] = coords[:, 0], coords[:, 1]
    write_table(plotdata.reset_index(drop=True), dest, 'pca_coordinates')
    fig, ax = plt.subplots(figsize=(7, 5.5))
    sns.scatterplot(data=plotdata, x='PC1', y='PC2', hue='group',
                    style='compartment' if 'compartment' in plotdata else None,
                    s=45, ax=ax)
    ax.set(xlabel=f'PC1 ({pca.explained_variance_ratio_[0]:.1%})',
           ylabel=f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
    save_plot(fig, dest, 'pca_by_group', settings.image_dpi)
    if len(meta) >= 8:
        from sklearn.manifold import TSNE
        perplexity = min(30., max(2., (len(meta)-1)/3.))
        tsne = TSNE(n_components=2, perplexity=perplexity, init='pca',
                    learning_rate='auto', random_state=settings.seed).fit_transform(x)
        plotdata['TSNE1'], plotdata['TSNE2'] = tsne[:, 0], tsne[:, 1]
        fig, ax = plt.subplots(figsize=(7, 5.5))
        sns.scatterplot(data=plotdata, x='TSNE1', y='TSNE2', hue='group',
                        style='compartment' if 'compartment' in plotdata else None, ax=ax)
        save_plot(fig, dest, 'tsne_by_group', settings.image_dpi)
        write_table(plotdata.reset_index(drop=True), dest, 'embeddings')
    means, design = subject_averages(expr, meta)
    n_features = min(settings.max_heatmap_features, expr.shape[0])
    ranked = expr.var(axis=1).nlargest(n_features).index
    heat = means.loc[:, ranked].T
    heat.index.name = 'feature'
    heat.columns = (design.sample_id.astype(str) + ' / ' + design.group.astype(str)).tolist()
    normalized = heat.sub(heat.mean(axis=1), axis=0).div(heat.std(axis=1).replace(0, 1), axis=0).fillna(0)
    grid = sns.clustermap(normalized, cmap='vlag', center=0, figsize=(max(8, len(design)*.35), 11),
                          xticklabels=True, yticklabels=True)
    save_plot(grid.fig, dest, 'top_variable_feature_heatmap', settings.image_dpi)
    if len(means) >= 3:
        correlations = means.T.corr().clip(-1, 1)
        correlations.index = heat.columns
        correlations.columns = heat.columns
        write_table(correlations.reset_index(names='sample'), dest, 'sample_correlations')
        off = correlations.to_numpy()[~np.eye(len(correlations), dtype=bool)]
        low, high = np.nanmin(off), np.nanmax(off)
        if abs(high - low) < 1e-10:
            low, high = low - .01, high + .01
        grid = sns.clustermap(correlations, cmap='RdBu_r', vmin=low, vmax=high,
                              center=(low+high)/2, figsize=(10, 9))
        save_plot(grid.fig, dest, 'sample_correlation_heatmap', settings.image_dpi)


def analyze(expr: pd.DataFrame, metadata: pd.DataFrame, out: Path,
            contrasts: list[tuple[str, str]], settings: Settings,
            strata: list[str] | None = None) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    expr, metadata = validate_inputs(expr, metadata)
    if not np.isfinite(expr.values).all():
        raise ValueError('Nonfinite expression values')
    write_table(metadata.reset_index(drop=True), out, 'validated_metadata')
    subject_summary = metadata.groupby('group').agg(
        n_segments=('segment_id', 'nunique'), n_independent_samples=('sample_id', 'nunique')).reset_index()
    write_table(subject_summary, out, 'study_summary')
    exploratory_visualizations(expr, metadata, out, settings)
    results = []
    for a, b in contrasts:
        if a == b:
            continue
        combinations = [({}, 'all')]
        for stratum in strata or []:
            if stratum not in metadata:
                LOG.warning('Unknown stratum ignored: %s', stratum)
                continue
            combinations.extend([({stratum: str(x)}, f'{stratum}_{sanitize(x)}')
                                 for x in metadata[stratum].dropna().unique()])
        for selector, label in combinations:
            result = differential_expression(expr, metadata, a, b, settings, selector)
            if result.empty:
                continue
            name = f'{sanitize(a)}_vs_{sanitize(b)}_{sanitize(label)}'
            write_table(result, out, f'deg_{name}')
            volcano(result, out, f'volcano_{name}', settings)
            results.append({'comparison': name, 'n_tested': int(result.p_value.notna().sum()),
                            'n_fdr_significant': int(result.FDR.lt(settings.alpha).sum())})
    report = {'expression_features': int(expr.shape[0]), 'segments': int(expr.shape[1]),
              'independent_samples': int(metadata.sample_id.nunique()),
              'comparisons': results, 'settings': settings.__dict__,
              'python_version': sys.version.split()[0],
              'limitations': ['Exploratory group comparisons',
                              'Subject-level replication used to avoid AOI pseudoreplication',
                              'Welch tests are not equivalent to mixed models or vendor DE results',
                              'Input expression must already be log2 normalized']}
    (out / 'analysis_report.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    LOG.info('Completed: %s', out.resolve())
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expression', type=Path, required=True)
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--out', type=Path, default=Path('spatial_analysis_results'))
    parser.add_argument('--contrasts', action='append', default=[], metavar='GROUP_A:GROUP_B')
    parser.add_argument('--stratify', action='append', default=[], metavar='COLUMN')
    parser.add_argument('--min-samples', type=int, default=3)
    parser.add_argument('--alpha', type=float, default=.05)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    expr = read_table(args.expression, index_col=0)
    meta = read_table(args.metadata)
    pairs = [tuple(pair.split(':', 1)) for pair in args.contrasts]
    if not pairs:
        groups = sorted(meta.group.dropna().astype(str).unique())
        from itertools import combinations
        pairs = list(combinations(groups, 2))
    analyze(expr, meta, args.out, pairs,
            Settings(alpha=args.alpha, min_samples_per_group=args.min_samples), args.stratify)


if __name__ == '__main__':
    main()
