#!/usr/bin/env python3
"""General-purpose extension runner for normalized spatial transcriptomics data.

Input expression: feature rows by segment columns, CSV. Input metadata: CSV
with segment_id, sample_id, group and optional covariates and spatial labels.
Do not apply count-based QC to normalized log expression.
"""
import argparse
import json
from pathlib import Path
import pandas as pd
from spatial_extensions import (validate_annotations, segment_qc, aggregate_probes,
    spatial_contrasts, fit_mixed_effects, compute_umap, deg_overlap)


def run(args):
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    expr=pd.read_csv(args.expression,index_col=0)
    expr.index=expr.index.astype(str)
    expr.columns=expr.columns.astype(str)
    meta=pd.read_csv(args.metadata)
    result=validate_annotations(expr,meta)
    (out/'annotation_validation.json').write_text(json.dumps(result,indent=2))
    if not result['is_valid']:raise ValueError('; '.join(result['issues']))
    meta=meta.set_index('segment_id',drop=False).loc[expr.columns]
    if args.input_scale=='counts':
        expr,meta,audit=segment_qc(expr,meta,min_total=args.min_total,
                                   min_detected=args.min_detected,
                                   max_mito_percent=args.max_mito_percent)
        audit.to_csv(out/'segment_qc_audit.csv',index_label='segment_id')
    if args.probe_annotations:
        if args.input_scale!='counts':raise ValueError('Probe aggregation must precede log transformation')
        annotation=pd.read_csv(args.probe_annotations)
        expr,audit=aggregate_probes(expr,annotation,probe_id=args.probe_id,
             gene_id=args.gene_id,min_total=args.min_probe_total,
             min_detected_fraction=args.min_probe_detection,method=args.aggregation)
        audit.to_csv(out/'probe_qc_audit.csv',index=False)
    expr.to_csv(out/'processed_expression.csv')
    meta.to_csv(out/'processed_metadata.csv',index=False)
    report={'segments':int(expr.shape[1]),'features':int(expr.shape[0]),'input_scale':args.input_scale,
            'statistics_run':False}
    if args.run_statistics:
        if args.input_scale !='log2':
            raise ValueError('Statistical modules require normalized log2 expression. Process raw counts separately.')
        pairs=[tuple(value.split(':',1)) for value in args.contrast]
        strata=[dict([value.split('=',1)]) for value in args.stratum]
        contrast=spatial_contrasts(expr,meta,args.condition,pairs,strata=strata or None,
                                  min_samples=args.min_samples,paired=args.paired) if pairs else pd.DataFrame()
        contrast.to_csv(out/'spatial_contrasts.csv',index=False)
        mixed=pd.DataFrame()
        if args.mixed_formula:
            mixed=fit_mixed_effects(expr,meta,formula=args.mixed_formula,
                       subject_column=args.subject_column,nested_columns=args.nested_column)
            mixed.to_csv(out/'mixed_effects.csv',index=False)
        if args.overlap and not contrast.empty and not mixed.empty:
            overlap=deg_overlap(contrast,mixed,alpha=args.alpha,minimum_effect=args.min_effect,
                                term=args.overlap_term)
            (out/'deg_overlap.json').write_text(json.dumps(overlap,indent=2))
        if args.umap:
            coordinates=compute_umap(expr,n_neighbors=args.umap_neighbors,seed=args.seed)
            coordinates.to_csv(out/'umap_coordinates.csv',index_label='segment_id')
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            fig,ax=plt.subplots(figsize=(7,5))
            for group, ids in meta.groupby(args.condition).groups.items():
                part=coordinates.loc[list(ids)]
                ax.scatter(part.UMAP1,part.UMAP2,label=str(group),s=25,alpha=.8)
            ax.set(xlabel='UMAP1',ylabel='UMAP2',title='UMAP of expression profiles')
            ax.legend(title=args.condition,bbox_to_anchor=(1.02,1),loc='upper left')
            fig.savefig(out/'umap.png',dpi=300,bbox_inches='tight')
            plt.close(fig)
        report.update(statistics_run=True,n_contrast_rows=int(len(contrast)),n_lmm_rows=int(len(mixed)))
    (out/'extended_run_manifest.json').write_text(json.dumps(report,indent=2))
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--expression',required=True)
    p.add_argument('--metadata',required=True)
    p.add_argument('--out',default='extended_results')
    p.add_argument('--input-scale',choices=['counts','log2'],required=True)
    p.add_argument('--min-total',type=float,default=0)
    p.add_argument('--min-detected',type=int,default=0)
    p.add_argument('--max-mito-percent',type=float)
    p.add_argument('--probe-annotations')
    p.add_argument('--probe-id',default='probe_id')
    p.add_argument('--gene-id',default='gene_id')
    p.add_argument('--min-probe-total',type=float,default=0)
    p.add_argument('--min-probe-detection',type=float,default=0)
    p.add_argument('--aggregation',choices=['sum','mean','median'],default='sum')
    p.add_argument('--run-statistics',action='store_true')
    p.add_argument('--condition',default='group')
    p.add_argument('--contrast',action='append',default=[],metavar='A:B')
    p.add_argument('--stratum',action='append',default=[],metavar='COLUMN=VALUE')
    p.add_argument('--min-samples',type=int,default=3)
    p.add_argument('--paired',action='store_true')
    p.add_argument('--mixed-formula',help='e.g., expression ~ C(group) + C(region)')
    p.add_argument('--subject-column',default='sample_id')
    p.add_argument('--nested-column',action='append',default=[])
    p.add_argument('--umap',action='store_true')
    p.add_argument('--umap-neighbors',type=int,default=15)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--overlap',action='store_true')
    p.add_argument('--overlap-term')
    p.add_argument('--alpha',type=float,default=.05)
    p.add_argument('--min-effect',type=float,default=0)
    run(p.parse_args())

if __name__=='__main__':main()
