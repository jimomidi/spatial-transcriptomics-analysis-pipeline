"""Synthetic-data checks for the reusable analysis pipeline."""
import tempfile
import unittest
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from spatial_analysis_core import Settings, validate_inputs, differential_expression
from geomx_spatial_pipeline import normalize_and_filter, build_segment_template, discover


class TestAnalysis(unittest.TestCase):
    def setUp(self):
        self.metadata = pd.DataFrame({
            'segment_id': [f'unit_{i}' for i in range(12)],
            'sample_id': [f'subject_{i // 2}' for i in range(12)],
            'group': ['A'] * 6 + ['B'] * 6,
            'compartment': ['X', 'Y'] * 6,
        })
        samples = self.metadata.segment_id.tolist()
        self.expr = pd.DataFrame(np.arange(48, dtype=float).reshape(4, 12) / 10 + 1,
                                 index=['G1', 'G2', 'G3', 'G4'], columns=samples)

    def test_sample_level_replicates(self):
        expr, meta = validate_inputs(self.expr, self.metadata)
        result = differential_expression(expr, meta, 'A', 'B', Settings())
        self.assertEqual(int(result.n_samples_a.iloc[0]), 3)
        self.assertEqual(int(result.n_samples_b.iloc[0]), 3)
        self.assertIn('FDR', result.columns)

    def test_reject_group_collision(self):
        invalid = self.metadata.copy()
        invalid.loc[0, 'group'] = 'B'
        with self.assertRaises(ValueError):
            validate_inputs(self.expr, invalid)

    def test_geomx_normalization(self):
        ids = ['segment_1', 'segment_2', 'segment_3', 'segment_4']
        source = pd.DataFrame({'SegmentDisplayName': ids,
                               'LOQ (Assay)': [2.0] * 4,
                               'SegmentLabel': ['C1'] * 4})
        props = pd.DataFrame({'TargetName': [f'G{i}' for i in range(12)],
                              'CodeClass': ['Endogenous'] * 12})
        counts = pd.DataFrame({'TargetName': props.TargetName,
                               **{sid: np.arange(1, 13) + j for j, sid in enumerate(ids)}})
        expression, qc, audit = normalize_and_filter(counts, props, source, .1, 75)
        self.assertEqual(expression.shape, (12, 4))
        self.assertEqual(qc.shape[0], 4)
        self.assertTrue(audit.retained.all())
        template = build_segment_template(source)
        self.assertTrue(template.sample_id.isna().all())

    def test_discover_compatible_workbook(self):
        with tempfile.TemporaryDirectory() as tmp:
            book = Path(tmp) / 'input.xlsx'
            with pd.ExcelWriter(book) as writer:
                for name in ('SegmentProperties', 'TargetProperties', 'TargetCountMatrix'):
                    pd.DataFrame({'example': [1]}).to_excel(writer, sheet_name=name, index=False)
            books, tables = discover(Path(tmp))
            self.assertEqual(books, [book])
            self.assertEqual(tables, [])


if __name__ == '__main__':
    unittest.main()
