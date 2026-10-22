"""Optional, study-independent quality control and statistical extensions.

Expression matrices have features in rows and segment identifiers in columns.
Statistical inputs must be normalized, log2-scale continuous expression.
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.formula.api import mixedlm
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


def validate_annotations(expression, metadata, required=('segment_id', 'sample_id', 'group'), nested=None):
    """Check unique IDs, alignment, required fields, and nesting consistency."""
    errors = []
    meta = metadata.copy()
    if not set(required).issubset(meta.columns):
        errors.append(f'Missing metadata fields: {sorted(set(required)-set(meta.columns))}')
    if expression.index.astype(str).duplicated().any() or expression.columns.astype(str).duplicated().any():
        errors.append('Duplicate feature or segment identifiers')
    if 'segment_id' in meta:
        meta['segment_id'] = meta.segment_id.astype('string').str.strip()
        if meta.segment_id.isna().any() or meta.segment_id.eq('').any() or meta.segment_id.duplicated().any():
            errors.append('Missing or duplicate segment IDs')
        absent = pd.Index(expression.columns.astype(str)).difference(meta.segment_id.dropna().astype(str))
        if len(absent):
            errors.append(f'{len(absent)} expression segments lack annotations')
    for key in required:
        if key in meta and (meta[key].isna().any() or meta[key].astype(str).str.strip().eq('').any()):
            errors.append(f'Missing/blank values in {key}')
    for child, parent in nested or [('segment_id','sample_id'), ('sample_id','group')]:
        if child in meta and parent in meta and meta.groupby(child, dropna=False)[parent].nunique(dropna=False).gt(1).any():
            errors.append(f'Inconsistent nesting: {child} maps to multiple {parent} values')
    numeric = expression.apply(pd.to_numeric, errors='coerce')
    if not np.isfinite(numeric.to_numpy()).all():
        errors.append('Nonfinite or nonnumeric expression values')
    return {'is_valid': not errors, 'issues': errors,
            'matching_segments': int(len(pd.Index(expression.columns.astype(str)).intersection(meta.segment_id.astype(str)))) if 'segment_id' in meta else 0}


def segment_qc(matrix, metadata, min_total=0, min_detected=0, max_mito_percent=None,
               mito_prefix='MT-', min_aligned_percent=None, min_nuclei=None,
               metrics=None, input_scale='counts'):
    """Audit configurable filters; total/detection metrics require nonnegative raw counts."""
    if input_scale != 'counts':
        raise ValueError('Segment count-based QC requires raw, nonnegative counts')
    if matrix.empty or (matrix.to_numpy() < 0).any() or not np.isfinite(matrix.to_numpy()).all():
        raise ValueError('Raw expression must be finite, nonempty, and nonnegative')
    if not matrix.columns.isin(metadata.index).all():
        raise ValueError('Metadata is missing segment identifiers')
    if min_total < 0 or min_detected < 0:
        raise ValueError('QC thresholds must be nonnegative')
    meta = metadata.loc[matrix.columns]
    audit = pd.DataFrame(index=matrix.columns)
    audit['total_counts'] = matrix.sum(axis=0)
    audit['detected_features'] = matrix.gt(0).sum(axis=0)
    passed = audit.total_counts.ge(min_total) & audit.detected_features.ge(min_detected)
    if max_mito_percent is not None:
        mito = matrix.index.astype(str).str.upper().str.startswith(mito_prefix.upper())
        if not mito.any():
            raise ValueError('Mitochondrial QC requested but no matching features exist')
        audit['mito_percent'] = matrix.loc[mito].sum(axis=0).div(audit.total_counts.replace(0,np.nan)).mul(100)
        passed &= audit.mito_percent.le(max_mito_percent)
    for threshold, key in ((min_aligned_percent, 'aligned_percent'), (min_nuclei,'nuclei')):
        if threshold is not None:
            if metrics is None or key not in metrics:
                raise ValueError(f'{key} threshold requires an external metrics table with this column')
            audit[key] = pd.to_numeric(metrics.reindex(matrix.columns)[key], errors='coerce')
            passed &= audit[key].ge(threshold)
    audit['qc_pass'] = passed
    audit['failed_criteria'] = ''
    criteria = {'min_total':audit.total_counts.ge(min_total), 'min_detected':audit.detected_features.ge(min_detected)}
    if max_mito_percent is not None: criteria['max_mito_percent'] = audit.mito_percent.le(max_mito_percent)
    if min_aligned_percent is not None: criteria['min_aligned_percent'] = audit.aligned_percent.ge(min_aligned_percent)
    if min_nuclei is not None: criteria['min_nuclei'] = audit.nuclei.ge(min_nuclei)
    audit['failed_criteria'] = [';'.join(k for k,v in criteria.items() if not bool(v.loc[s])) for s in audit.index]
    return matrix.loc[:, passed], meta.loc[passed], audit


def aggregate_probes(probes, annotations, probe_id='probe_id', gene_id='gene_id',
                     min_total=0, min_detected_fraction=0.0, method='sum',
                     exclusion_column=None):
    """Filter raw probes, map identifiers uniquely, and aggregate to gene counts.

    Sum is appropriate for counts when probes target nonoverlapping molecules;
    selecting the correct aggregation depends on the assay and annotation.
    """
    if method not in ('sum','mean','median') or not 0 <= min_detected_fraction <= 1:
        raise ValueError('Invalid aggregation method or detection fraction')
    if probes.index.has_duplicates or annotations[probe_id].duplicated().any():
        raise ValueError('Probe IDs must be unique')
    if not {probe_id,gene_id}.issubset(annotations):
        raise ValueError('Probe mapping lacks required columns')
    if (probes.to_numpy() < 0).any() or not np.isfinite(probes.to_numpy()).all():
        raise ValueError('Probe matrix must contain nonnegative finite counts')
    ann = annotations.set_index(probe_id)
    ann.index = ann.index.astype(str)
    keep = probes.sum(axis=1).ge(min_total) & probes.gt(0).mean(axis=1).ge(min_detected_fraction)
    if exclusion_column:
        if exclusion_column not in ann: raise ValueError('Missing exclusion annotation')
        excluded = ann[exclusion_column].reindex(probes.index).fillna(False).astype(str).str.lower().isin(('true','1','yes','fail'))
        keep &= ~excluded
    mapping = ann[gene_id].reindex(probes.index)
    keep &= mapping.notna() & mapping.astype(str).str.strip().ne('')
    qc = pd.DataFrame({'probe_id':probes.index, 'total_counts':probes.sum(axis=1).values,
                       'detected_fraction':probes.gt(0).mean(axis=1).values,
                       'gene_id':mapping.values,'retained':keep.values})
    filtered = probes.loc[keep].copy()
    if filtered.empty: raise ValueError('No probes passed QC')
    filtered['__target__'] = mapping.loc[filtered.index].astype(str).to_numpy()
    grouped = filtered.groupby('__target__', sort=True)
    result = getattr(grouped,method)(numeric_only=True)
    result.index.name = 'gene_id'
    return result, qc


def spatial_contrasts(expr, meta, contrast_column, pairs, sample_column='sample_id',
                      strata=None, min_samples=3, paired=False):
    """Welch or paired subject-level tests on log2 expression; no segment pseudoreplication."""
    if contrast_column not in meta or sample_column not in meta:
        raise ValueError('Contrast and biological sample columns are required')
    if not expr.columns.isin(meta.index).all(): raise ValueError('Unmatched segments')
    output = []
    for stratum in (strata or [{}]):
        mask = pd.Series(True,index=meta.index)
        for col,value in stratum.items():
            if col not in meta: raise ValueError(f'Unknown stratum: {col}')
            mask &= meta[col].eq(value)
        m = meta.loc[mask & meta[contrast_column].isin([v for pair in pairs for v in pair])]
        for a,b in pairs:
            subset = m.loc[m[contrast_column].isin((a,b))]
            if subset.empty: continue
            frame = expr[subset.index].T.join(subset[[sample_column,contrast_column]])
            collapsed = frame.groupby([sample_column,contrast_column],observed=True)[expr.index.tolist()].mean()
            groups = collapsed.index.get_level_values(contrast_column)
            aa = collapsed.loc[groups == a]; bb = collapsed.loc[groups == b]
            if min(len(aa),len(bb)) < min_samples: continue
            if paired:
                shared = aa.index.get_level_values(sample_column).intersection(bb.index.get_level_values(sample_column))
                if len(shared) < min_samples: continue
                xa = aa.droplevel(contrast_column).loc[shared].to_numpy()
                xb = bb.droplevel(contrast_column).loc[shared].to_numpy()
                p = stats.ttest_rel(xa,xb,axis=0).pvalue
            else:
                xa,xb = aa.to_numpy(),bb.to_numpy()
                p = stats.ttest_ind(xa,xb,axis=0,equal_var=False).pvalue
            p = np.asarray(p,dtype=float)
            q = np.full(len(p),np.nan)
            valid = np.isfinite(p)
            if valid.any(): q[valid] = multipletests(p[valid],method='fdr_bh')[1]
            out = pd.DataFrame({'feature':expr.index,'contrast':f'{a}_vs_{b}',
                     'estimate_log2':xa.mean(axis=0)-xb.mean(axis=0),
                     'p_value':p,'FDR':q,'n_a':len(xa),'n_b':len(xb),
                     'method':'paired_subject_t' if paired else 'subject_Welch'})
            for k,v in stratum.items():out[k]=v
            output.append(out)
    return pd.concat(output,ignore_index=True) if output else pd.DataFrame(columns=['feature','contrast','estimate_log2','p_value','FDR'])


def fit_mixed_effects(expr, meta, formula='expression ~ C(group)', subject_column='sample_id',
                      nested_columns=(), features=None, min_subjects=4):
    """Per-feature Gaussian random-intercept LMM with optional variance components.

    The subject random intercept handles repeated measures. Optional nested IDs
    become variance components and must identify units globally or be composite IDs.
    All fit errors and convergence flags are explicitly reported.
    """
    if not formula.lstrip().startswith('expression ~'):
        raise ValueError('The formula response must be expression')
    needed = [subject_column,*nested_columns]
    if any(k not in meta for k in needed): raise ValueError('Missing random-effect metadata')
    if meta[subject_column].nunique() < min_subjects:
        raise ValueError('Insufficient independent samples for mixed-effects analysis')
    vc = {f'vc_{i}':f'0 + C(Q("{col}"))' for i,col in enumerate(nested_columns)}
    rows=[]
    for feature in (features if features is not None else list(expr.index)):
        data=meta.loc[expr.columns].copy()
        data['expression'] = expr.loc[feature].to_numpy(dtype=float)
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                fitted = mixedlm(formula,data=data,groups=data[subject_column],
                                 vc_formula=vc or None,re_formula='1').fit(reml=False,method='lbfgs',disp=False)
            converged=bool(fitted.converged)
            for term in fitted.fe_params.index:
                if term=='Intercept':continue
                rows.append({'feature':feature,'term':term,'estimate':float(fitted.fe_params[term]),
                             'std_error':float(fitted.bse[term]),'p_value':float(fitted.pvalues[term]),
                             'converged':converged,'warning':' | '.join(str(w.message) for w in caught)[:500],
                             'error':''})
        except Exception as exc:
            rows.append({'feature':feature,'term':'','estimate':np.nan,'std_error':np.nan,
                         'p_value':np.nan,'converged':False,'warning':'','error':str(exc)[:500]})
    out=pd.DataFrame(rows,columns=['feature','term','estimate','std_error','p_value','converged','warning','error'])
    out['FDR']=np.nan
    for term,idx in out.loc[out.converged & out.p_value.notna()].groupby('term').groups.items():
        out.loc[idx,'FDR']=multipletests(out.loc[idx,'p_value'],method='fdr_bh')[1]
    return out


def compute_umap(expr, n_neighbors=15, min_dist=.1, seed=42, max_features=2000):
    """Compute genuine UMAP; do not silently substitute PCA."""
    try:
        import umap
    except ImportError as exc:
        raise ImportError('Install umap-learn to enable UMAP') from exc
    if expr.shape[1]<4: raise ValueError('UMAP requires at least four segments')
    if n_neighbors < 2: raise ValueError('n_neighbors must be at least two')
    selected=expr.var(axis=1).nlargest(max_features).index
    x=StandardScaler().fit_transform(expr.loc[selected].T)
    coords=umap.UMAP(n_neighbors=min(n_neighbors,expr.shape[1]-1),min_dist=min_dist,random_state=seed).fit_transform(x)
    return pd.DataFrame(coords,index=expr.columns,columns=['UMAP1','UMAP2'])


def deg_overlap(contrasts, mixed, alpha=.05, minimum_effect=0.0, term=None):
    """Overlap only matched feature sets; an LMM term filter is recommended."""
    if contrasts.empty or mixed.empty: return {'n_contrast':0,'n_mixed':0,'n_overlap':0,'features':[]}
    if term is None and mixed.loc[mixed.converged,'term'].nunique()>1:
        raise ValueError('Specify one LMM term to avoid combining unrelated coefficients')
    selected=mixed if term is None else mixed.loc[mixed.term.eq(term)]
    a=set(contrasts.loc[contrasts.FDR.le(alpha) & contrasts.estimate_log2.abs().ge(minimum_effect),'feature'])
    b=set(selected.loc[selected.converged & selected.FDR.le(alpha),'feature'])
    overlap=sorted(a & b)
    return {'n_contrast':len(a),'n_mixed':len(b),'n_overlap':len(overlap),'features':overlap,'lmm_term':term}
