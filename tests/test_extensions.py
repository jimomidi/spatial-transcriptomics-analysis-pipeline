"""Synthetic regression tests for generic extension modules."""
import unittest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from spatial_extensions import (validate_annotations,segment_qc,aggregate_probes,
                                spatial_contrasts,fit_mixed_effects,deg_overlap)

class TestExtensions(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(7)
        ids=[f's_{i}_{r}' for i in range(12) for r in range(2)]
        self.meta=pd.DataFrame({'segment_id':ids,'sample_id':[f's_{i}' for i in range(12) for _ in range(2)],
                                'group':['A' if i<6 else 'B' for i in range(12) for _ in range(2)],
                                'region':['X','Y']*12}).set_index('segment_id',drop=False)
        self.expr=pd.DataFrame(rng.normal(5,1,(6,24)),index=[f'feature{i}' for i in range(6)],columns=ids)
        self.expr.loc['feature0',self.meta.group.eq('B')]+=2
    def test_validation_rejects_mismatched(self):
        self.assertTrue(validate_annotations(self.expr,self.meta.reset_index(drop=True))['is_valid'])
        m=self.meta.reset_index(drop=True).iloc[:-1]
        self.assertFalse(validate_annotations(self.expr,m)['is_valid'])
    def test_segment_qc_audit(self):
        counts=self.expr.clip(lower=0).round()
        counts.iloc[:,0]=0
        mat,meta,audit=segment_qc(counts,self.meta,min_total=1,min_detected=1)
        self.assertEqual(mat.shape[1],23)
        self.assertIn('min_total',audit.iloc[0].failed_criteria)
        self.assertEqual(len(meta),23)
    def test_probe_aggregation(self):
        counts=self.expr.clip(lower=0).round()
        ann=pd.DataFrame({'probe_id':counts.index,'gene_id':['G1']*3+['G2']*3})
        mat,qc=aggregate_probes(counts,ann,method='sum')
        self.assertEqual(mat.shape,(2,24))
        self.assertEqual(len(qc),6)
    def test_contrasts_use_subjects(self):
        res=spatial_contrasts(self.expr,self.meta,'group',[('A','B')],min_samples=3)
        self.assertEqual(int(res.n_a.iloc[0]),6)
        self.assertEqual(int(res.n_b.iloc[0]),6)
        self.assertLess(res.loc[res.feature.eq('feature0'),'p_value'].iloc[0],.05)
    def test_paired_region(self):
        res=spatial_contrasts(self.expr,self.meta,'region',[('X','Y')],paired=True,min_samples=3)
        self.assertEqual(int(res.n_a.iloc[0]),12)
    def test_mixed_output(self):
        result=fit_mixed_effects(self.expr,self.meta,'expression ~ C(group)',features=['feature0'])
        self.assertIn('converged',result)
        self.assertIn('FDR',result)
        self.assertEqual(result.feature.iloc[0],'feature0')
    def test_overlap(self):
        a=pd.DataFrame({'feature':['f1'],'FDR':[.01],'estimate_log2':[1.]})
        b=pd.DataFrame({'feature':['f1'],'term':['C(group)[T.B]'],'FDR':[.01],'converged':[True]})
        self.assertEqual(deg_overlap(a,b,term='C(group)[T.B]')['n_overlap'],1)

if __name__=='__main__':unittest.main()
