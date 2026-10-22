#!/usr/bin/env python3
"""GeoMx DSP spatial transcriptomics: reusable ingestion, QC and analysis.

Designed for GeoMx-style workbooks with SegmentProperties, TargetProperties,
and TargetCountMatrix sheets. Discovers inputs without relying on particular
laboratory names, study groups, organism, or fixed sample counts.

Quick start:
  python geomx_spatial_pipeline.py --root INPUT_DIRECTORY --out OUTPUT_DIRECTORY

Review the generated segment_metadata_template.csv, add sample_id and group,
then run validated analyses with:
  python geomx_spatial_pipeline.py --root /path/to/study \
      --metadata /path/to/verified_metadata.csv --out /path/to/results

Only verified biological sample mapping should be supplied for inference.
If metadata is not supplied, QC, counts, and normalization are still produced.
This program does not claim to reproduce all proprietary vendor algorithms.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from spatial_analysis_core import Settings, analyze, read_table, write_table, save_plot

LOG = logging.getLogger('geomx_pipeline')

REQUIRED_SHEETS = {'SegmentProperties', 'TargetProperties', 'TargetCountMatrix'}
METADATA_KEYS = {'segment_id', 'sample_id', 'group'}


def discover(root: Path) -> tuple[list[Path], list[Path]]:
    workbooks, vendor_tables = [], []
    for p in root.rglob('*'):
        if not p.is_file() or p.name.startswith('~$'):
            continue
        if p.suffix.lower() == '.xlsx':
            workbooks.append(p)
        if p.suffix.lower() in {'.csv', '.tsv'}:
            vendor_tables.append(p)
    eligible = []
    for book in workbooks:
        try:
            sheets = set(pd.ExcelFile(book).sheet_names)
            if REQUIRED_SHEETS.issubset(sheets):
                eligible.append(book)
        except Exception as exc:
            LOG.warning('Skipping unreadable workbook %s: %s', book, exc)
    return eligible, vendor_tables


def read_workbook(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    with pd.ExcelFile(path) as book:
        meta = pd.read_excel(book, 'SegmentProperties')
        properties = pd.read_excel(book, 'TargetProperties')
        counts = pd.read_excel(book, 'TargetCountMatrix')
    for frame in (meta, properties, counts):
        frame.columns = frame.columns.astype(str).str.strip()
    required_meta = {'SegmentDisplayName'}
    if not required_meta.issubset(meta.columns):
        raise ValueError(f'Workbook missing segment metadata: {required_meta - set(meta.columns)}')
    if 'TargetName' not in counts or 'TargetName' not in properties:
        raise ValueError('TargetName is required in both the count and target-property sheets')
    return meta, properties, counts


def build_segment_template(source: pd.DataFrame) -> pd.DataFrame:
    """Never infer independent sample identity from a slide or ROI label."""
    template = pd.DataFrame({'segment_id': source.SegmentDisplayName.astype(str).str.strip(),
                             'sample_id': pd.Series(pd.NA, index=source.index, dtype='string'),
                             'group': pd.Series(pd.NA, index=source.index, dtype='string')})
    available = {'SlideName': 'slide', 'ROILabel': 'roi', 'SegmentLabel': 'compartment',
                 'ROIComments': 'annotation_note', 'QCFlags': 'vendor_qc_flags'}
    for source_name, output_name in available.items():
        if source_name in source:
            template[output_name] = source[source_name].values
    for optional in ('region', 'site', 'batch'):
        if optional not in template:
            template[optional] = pd.NA
    return template


def verify_metadata(metadata: pd.DataFrame, available: pd.Index) -> pd.DataFrame:
    absent = METADATA_KEYS - set(metadata.columns)
    if absent:
        raise ValueError(f'Metadata missing required fields: {sorted(absent)}')
    if metadata.segment_id.duplicated().any():
        raise ValueError('Verified metadata contains duplicate segment IDs')
    missing = available.difference(metadata.segment_id.astype(str))
    if len(missing):
        raise ValueError(f'Missing metadata for {len(missing)} segments; examples: {missing[:4].tolist()}')
    for name in ('sample_id', 'group'):
        if metadata[name].isna().any() or metadata[name].astype(str).str.strip().eq('').any():
            raise ValueError(f'All entries in {name} must be independently verified and nonempty')
    return metadata


def normalize_and_filter(raw: pd.DataFrame, properties: pd.DataFrame,
                         source_meta: pd.DataFrame, prevalence: float,
                         q3_percentile: float) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Endogenous-target Q3 normalization with per-segment LOQ filtering audit.

    Input LOQ is taken from the instrument workbook; it is never silently
    reconstructed. Q3 factors are calculated from all endogenous targets,
    before gene prevalence filtering. No imputation is performed.
    """
    if not (0 < prevalence <= 1 and 0 < q3_percentile < 100):
        raise ValueError('Prevalence and quantile must be within valid ranges')
    raw = raw.copy().set_index('TargetName')
    if raw.index.has_duplicates:
        raise ValueError('Duplicate feature IDs are not supported; inspect the workbook')
    raw.columns = raw.columns.astype(str).str.strip()
    segments = pd.Index(source_meta.SegmentDisplayName.astype(str).str.strip())
    if segments.has_duplicates:
        raise ValueError('Duplicate segment display names detected')
    absent = segments.difference(raw.columns)
    if len(absent):
        raise ValueError(f'Count matrix lacks {len(absent)} segments')
    raw = raw.loc[:, segments].apply(pd.to_numeric, errors='coerce')
    if raw.isna().any().any() or (raw < 0).any().any():
        raise ValueError('Count matrix contains missing or negative values')
    classes = properties.set_index('TargetName')['CodeClass'].astype(str).str.casefold()
    targets = classes.index[classes.eq('endogenous')].intersection(raw.index)
    if len(targets) < 10:
        raise ValueError('Fewer than ten endogenous targets identified')
    raw = raw.loc[targets]
    loq_columns = [name for name in source_meta if str(name).strip().startswith('LOQ (')]
    if len(loq_columns) != 1:
        raise ValueError(f'Exactly one LOQ column expected, found {len(loq_columns)}')
    source_meta = source_meta.copy()
    source_meta['segment_id'] = source_meta.SegmentDisplayName.astype(str).str.strip()
    loq = pd.to_numeric(source_meta.set_index('segment_id').loc[segments, loq_columns[0]],
                        errors='coerce')
    if loq.isna().any():
        raise ValueError('LOQ measurements missing; do not replace with arbitrary thresholds')
    detected = raw.ge(loq, axis='columns')
    target_prevalence = detected.mean(axis=1)
    selected = target_prevalence.ge(prevalence)
    q3 = raw.quantile(q3_percentile / 100, axis='index')
    if q3.isna().any() or (q3 <= 0).any():
        raise ValueError('Nonpositive Q3 normalization factors')
    # Scaling to median Q3 retains a comparable count-like expression scale.
    normalized = raw.loc[selected].div(q3, axis='columns') * q3.median()
    log_expression = np.log2(normalized + 1)
    qc = pd.DataFrame({'segment_id': segments, 'LOQ': loq.to_numpy(),
                       'endogenous_targets_detected': detected.sum(axis=0).to_numpy(),
                       'Q3_factor': q3.to_numpy(),
                       'endogenous_targets_total': len(raw),
                       'retained_targets': int(selected.sum())})
    audit = pd.DataFrame({'feature': target_prevalence.index,
                          'fraction_segments_above_LOQ': target_prevalence.values,
                          'retained': selected.values})
    return log_expression, qc, audit


def qc_figures(source_meta: pd.DataFrame, qc: pd.DataFrame,
               audit: pd.DataFrame, out: Path, settings: Settings) -> None:
    sns.set_theme(style='ticks', context='paper')
    combined = source_meta.copy()
    combined['segment_id'] = combined.SegmentDisplayName.astype(str).str.strip()
    combined = combined.merge(qc, on='segment_id', validate='one_to_one')
    numeric = ['RawReads', 'AlignedReads', 'SequencingSaturation',
               'AOINucleiCount', 'AOISurfaceArea', 'Q3_factor',
               'endogenous_targets_detected']
    numeric = [key for key in numeric if key in combined]
    if 'SegmentLabel' in combined:
        for col in numeric:
            combined[col] = pd.to_numeric(combined[col], errors='coerce')
            fig, ax = plt.subplots(figsize=(7, 4.8))
            sns.boxplot(data=combined, x='SegmentLabel', y=col, ax=ax)
            sns.stripplot(data=combined, x='SegmentLabel', y=col, ax=ax,
                          color='#333333', size=2, alpha=.35)
            ax.set(xlabel='Segment type', title=f'Distribution of {col}')
            save_plot(fig, out, f'qc_{col}', settings.image_dpi)
    fig, ax = plt.subplots(figsize=(7, 4.8))
    sns.histplot(audit.fraction_segments_above_LOQ, bins=25, ax=ax, color='#537895')
    ax.set(xlabel='Fraction of segments above LOQ', ylabel='Number of targets',
           title='Feature detection prevalence')
    save_plot(fig, out, 'feature_detection_prevalence', settings.image_dpi)
    write_table(combined, out, 'qc_metrics_per_segment')
    if 'QCFlags' in combined:
        flags = combined.QCFlags.fillna('No flags').astype(str).value_counts().rename_axis('QCFlags').reset_index(name='n_segments')
        write_table(flags, out, 'vendor_qc_flag_summary')


def export_vendor_inventory(vendor_files: list[Path], root: Path, out: Path) -> None:
    records = []
    for file in vendor_files:
        relative = str(file.relative_to(root))
        row = {'relative_path': relative, 'size_bytes': file.stat().st_size,
               'format': file.suffix.lower(), 'columns': '', 'n_rows': None,
               'has_estimate': False, 'has_p_value': False, 'has_fdr': False}
        try:
            data = read_table(file, nrows=8)
            columns = [str(x) for x in data.columns]
            row['columns'] = ' | '.join(columns)
            tokens = [x.lower() for x in columns]
            row['has_estimate'] = any('estimate' in x or 'logfc' in x or 'log2fc' in x for x in tokens)
            row['has_p_value'] = any('p.value' in x or 'pvalue' in x or 'pr(>' in x or x == 'p' for x in tokens)
            row['has_fdr'] = any('fdr' in x or 'adj.p' in x or 'padj' in x or 'qvalue' in x for x in tokens)
        except Exception as exc:
            row['read_error'] = str(exc)[:150]
        records.append(row)
    write_table(pd.DataFrame(records, columns=['relative_path', 'size_bytes', 'format',
                'columns', 'n_rows', 'has_estimate', 'has_p_value', 'has_fdr', 'read_error']),
                out, 'vendor_table_inventory')


def run(args: argparse.Namespace) -> None:
    root = args.root.expanduser().resolve()
    output = args.out.expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(root)
    if root == output or root in output.parents and len(output.parts) <= len(root.parts):
        raise ValueError('Output directory must not replace source directory')
    output.mkdir(parents=True, exist_ok=True)
    workbooks, vendor_tables = discover(root)
    export_vendor_inventory([p for p in vendor_tables if output not in p.parents], root, output)
    if args.workbook:
        workbook = args.workbook.expanduser().resolve()
    elif len(workbooks) == 1:
        workbook = workbooks[0]
    elif not workbooks:
        raise FileNotFoundError('No workbook with all required GeoMx sheets was discovered')
    else:
        raise ValueError('Multiple eligible GeoMx workbooks; supply --workbook explicitly:\n' +
                         '\n'.join(map(str, workbooks)))
    LOG.info('Reading %s', workbook)
    source_meta, target_properties, raw = read_workbook(workbook)
    template = build_segment_template(source_meta)
    write_table(template, output, 'segment_metadata_template')
    settings = Settings(alpha=args.alpha, min_samples_per_group=args.min_samples,
                        image_dpi=args.dpi)
    expression, qc, audit = normalize_and_filter(raw, target_properties, source_meta,
                                                 args.min_prevalence, args.q3_percentile)
    write_table(audit, output, 'target_detection_audit')
    write_table(qc, output, 'normalization_factor_audit')
    expression.to_csv(output / 'log2_q3_normalized_expression.csv', index=True)
    qc_figures(source_meta, qc, audit, output, settings)
    results = {'source_format': 'GeoMx-style Excel workbook',
               'workbook_sha256': hashlib.sha256(workbook.read_bytes()).hexdigest(),
               'n_segments': int(expression.shape[1]),
               'n_endogenous_features_retained': int(expression.shape[0]),
               'q3_percentile': args.q3_percentile,
               'min_prevalence': args.min_prevalence,
               'metadata_verified_externally': bool(args.metadata),
               'statistical_inference': 'not_run'}
    if args.metadata:
        verified = verify_metadata(read_table(args.metadata), expression.columns)
        joined = verified.merge(template.drop(columns=['sample_id', 'group']),
                                on='segment_id', how='left', suffixes=('', '_source'),
                                validate='one_to_one')
        for key in ('compartment', 'region', 'site', 'batch'):
            backup = f'{key}_source'
            if backup in joined:
                joined[key] = joined[key].fillna(joined[backup]) if key in joined else joined[backup]
        if 'vendor_qc_flags' in joined and args.exclude_vendor_flags:
            joined['qc_pass'] = joined.vendor_qc_flags.isna() | joined.vendor_qc_flags.astype(str).str.strip().eq('')
        from itertools import combinations
        groups = sorted(joined.group.astype(str).unique())
        contrasts = ([tuple(text.split(':', 1)) for text in args.contrasts]
                     if args.contrasts else list(combinations(groups, 2)))
        results['analysis'] = analyze(expression, joined, output / 'statistical_analysis',
                                      contrasts, settings,
                                      [col for col in args.stratify if col in joined])
        results['statistical_inference'] = 'exploratory_subject_level'
    else:
        LOG.warning('Metadata not supplied. QC and normalization completed; inference skipped.')
        (output / 'METADATA_VERIFICATION_REQUIRED.txt').write_text(
            'Assign sample_id and group independently for every segment in the '
            'segment_metadata_template.csv file. Review the biological codebook '
            'and all ROI identifiers before rerunning with --metadata. '
            'Do not use the slide name as the biological sample identifier.\n',
            encoding='utf-8')
    (output / 'run_manifest.json').write_text(json.dumps(results, indent=2), encoding='utf8')
    LOG.info('Saved all outputs to %s', output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--out', type=Path, default=Path('geomx_results'))
    parser.add_argument('--workbook', type=Path)
    parser.add_argument('--metadata', type=Path, help='Independent sample mapping validated against a codebook')
    parser.add_argument('--contrasts', action='append', default=[], metavar='GROUP_A:GROUP_B')
    parser.add_argument('--stratify', action='append', default=['compartment', 'region'])
    parser.add_argument('--min-prevalence', type=float, default=.10)
    parser.add_argument('--q3-percentile', type=float, default=75.)
    parser.add_argument('--min-samples', type=int, default=3)
    parser.add_argument('--alpha', type=float, default=.05)
    parser.add_argument('--dpi', type=int, default=600)
    parser.add_argument('--exclude-vendor-flags', action='store_true')
    arguments = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    try:
        run(arguments)
    except Exception:
        LOG.exception('Analysis failed')
        sys.exit(1)


if __name__ == '__main__':
    main()
