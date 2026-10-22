# Spatial Transcriptomics Analysis Pipeline

A reusable, study-independent Python workflow for **spatial transcriptomics**, with a generic expression-analysis engine and a companion adapter for **GeoMx DSP-style** Excel workbooks. The repository includes data checks, descriptive quality control, LOQ-based feature detection, Q3 normalization, sample-level exploratory statistics, and publication-oriented plots.

The scripts do **not** contain laboratory-specific file paths, sample names, disease labels, organism assumptions, or hard-coded sample counts. Supply all data locations and analysis design using command-line options and verified metadata.

## Workflow

| Component | Script | Role |
|---|---|---|
| General analysis | `spatial_analysis_core.py` | Validate expression and sample metadata, summarize replication, PCA/t-SNE, correlation and heatmaps, subject-level differential expression, BH FDR, volcano plots |
| GeoMx ingestion | `geomx_spatial_pipeline.py` | Discover compatible workbooks, audit measured LOQ and endogenous features, Q3-normalize, export QC visualizations and a metadata template; call general analysis after metadata verification |

The GeoMx adapter **imports** the general engine. Keep both scripts in the same directory.

## Setup

- Python 3.10 or newer
- Install dependencies in a clean Python environment:

```bash
python -m pip install -r requirements.txt
```

## Repository layout

```text
.
├── spatial_analysis_core.py
├── geomx_spatial_pipeline.py
├── tests/
│   └── test_pipeline.py
├── data/
│   └── .gitkeep
├── results/
│   └── .gitkeep
├── requirements.txt
├── README.md
├── LICENSE
└── .gitignore
```

Input data and analysis results are excluded from version control. This repository ships with **no original patient, animal, or vendor data**.

## Option A: Generic expression analysis

Input expression matrix: CSV, TSV, or Excel; **features in rows**, **segment IDs in columns**, **feature IDs in the first column**. Provide already log2-normalized expression values, with no missing or infinite cells.

Input metadata: CSV, TSV, or Excel; **one row per segment**:

| Column | Requirement | Meaning |
|---|---|---|
| `segment_id` | Required | Exactly matches expression-matrix column ID |
| `sample_id` | Required | Independent biological sample/subject ID, not an AOI or slide ID |
| `group` | Required | Experimental group label |
| `compartment` | Optional | Cell type or compartment for faceting/stratification |
| `region` | Optional | Histological region |
| `site` | Optional | Sampling site |
| `batch` | Optional | Instrument batch or slide |
| `qc_pass` | Optional | `true`/`false`, `1`/`0`, or `pass`/`fail` |

Run from the repository root using **your own** input arguments:

```bash
python spatial_analysis_core.py \
  --expression INPUT_EXPRESSION_TABLE \
  --metadata VERIFIED_METADATA_TABLE \
  --out OUTPUT_DIRECTORY \
  --contrasts GroupA:GroupB \
  --stratify compartment
```

Leave `--contrasts` unset for all pairwise comparisons. `--stratify` can be supplied more than once. Options include `--min-samples` (default 3) and `--alpha` (default 0.05).

### Generic outputs

- Validated segment metadata and group/sample counts (CSV and Excel)
- PCA and, with enough segments, t-SNE visualizations
- Highly variable feature heatmap and sample-level correlation matrix
- Pairwise gene/feature comparison tables with `estimate_log2`, `p_value`, `FDR`, and independent replicate counts
- Volcano plots and a machine-readable JSON analysis report
- Plot exports in PNG, PDF, and SVG (600 DPI raster setting)

## Option B: GeoMx DSP workbook workflow

The input directory must contain a compatible Excel workbook with these worksheets:

- `SegmentProperties`: includes `SegmentDisplayName`, QC fields, and one `LOQ (` threshold column
- `TargetProperties`: includes `TargetName` and `CodeClass`
- `TargetCountMatrix`: first column `TargetName`; remaining columns are segment IDs

The program searches recursively for compatible workbooks. When multiple compatible workbooks exist, provide `--workbook` to make the selection explicit. It also lists CSV/TSV analysis exports as a **vendor-table inventory**; it does not blindly reinterpret them as differential-expression results.

**Step 1.** Run QC and normalization, without inferential group analysis:

```bash
python geomx_spatial_pipeline.py \
  --root INPUT_DIRECTORY \
  --out OUTPUT_DIRECTORY
```

**Step 2.** Review the generated `segment_metadata_template.csv`; independently fill in `sample_id` and `group`, and correct optional `compartment`, `region` and `site` labels from your verified study codebook. Save a copy outside the results directory or with another appropriate name.

**Step 3.** Supply verified metadata for exploratory analysis:

```bash
python geomx_spatial_pipeline.py \
  --root INPUT_DIRECTORY \
  --metadata VERIFIED_METADATA_TABLE \
  --out OUTPUT_DIRECTORY \
  --contrasts GroupA:GroupB \
  --stratify compartment
```

Additional options: `--min-prevalence` (default 0.10), `--q3-percentile` (default 75), `--min-samples` (default 3), `--alpha` (default 0.05), `--dpi` (default 600), `--exclude-vendor-flags` (optional exclusion of flagged segments), and `--workbook` (manual workbook selection).

### GeoMx outputs

- Segment metadata template for independent verification
- Detection-prevalence and LOQ-retention audit
- Per-segment Q3 factors and log2(Q3-normalized expression + 1)
- Read-depth, sequencing saturation, nuclei and other available QC plots
- Vendor QC flag and supplied table inventories
- If verified metadata is provided: results generated by the general analysis engine in a statistical-analysis subdirectory
- JSON run manifest with parameters and source-file checksum

## Statistical interpretation and limitations

1. **Independent biological units:** multiple AOIs from the same subject are averaged before a two-group Welch test. This avoids simple AOI pseudoreplication, but does not replace a study-specific hierarchical model.
2. **Effect estimate:** `estimate_log2` is the difference in mean log2-scale expression (first group minus second group). It is not a fitted vendor coefficient and should not be described as one.
3. **Multiplicity:** Benjamini–Hochberg correction is applied within each completed comparison. Interpret small-group results cautiously; the workflow is exploratory.
4. **LOQ:** supplied instrument LOQ values are used; missing LOQs result in an error. A target must be detected in at least the specified fraction of segments to pass filtering.
5. **Normalization:** Q3 factors are computed from all endogenous targets before the prevalence filter; this is a documented generic implementation and may differ from vendor software details.
6. **Provenance:** sample identities and biological meaning of ROI labels **cannot** be reconstructed safely from position, slide names or naming conventions alone. No group inference is run before complete verified metadata is supplied.
7. **Generalizability:** the scripts do not replace vendor probe QC, gene-specific design adjustment, negative-probe normalization, spatial autocorrelation models, or external biological validation. Compatibility is focused on the documented workbook structure, not every GeoMx software export variant.
8. **Privacy:** source data are not bundled. Before publicly sharing outputs, review filenames, segment identifiers, report fields and metadata for restricted or identifiable information.

## Testing

```bash
python -m unittest discover -s tests -v
```

Tests use synthetic data only and check metadata validity, sample-level replication, normalization, and input discovery. Run a study-specific benchmark against reference outputs before relying on findings for publication.

## License

Distributed under the [MIT License](LICENSE). See `LICENSE` for terms.

## Optional six-module extension runner

`extended_pipeline.py` adds validated annotations, configurable segment QC,
probe-to-gene aggregation, subject-aware contrasts, Gaussian linear mixed-effects
models, real UMAP embeddings and optional DEG overlap analysis. It uses the
independent module `spatial_extensions.py`. The original core and GeoMx workbook
workflow remain available unchanged.

**Input:** Expression CSV has features in rows and segment IDs in columns;
metadata CSV includes `segment_id`, `sample_id`, `group`. An optional probe
annotation CSV maps `probe_id` to `gene_id`. All identifiers, groups, paths and
spatial strata are user-configurable, with no study-specific values.

For raw counts (QC and optional probe aggregation only):

```bash
python extended_pipeline.py --expression INPUT_COUNTS.csv --metadata INPUT_METADATA.csv \
  --input-scale counts --min-total 100 --min-detected 10 \
  --probe-annotations INPUT_PROBE_MAPPING.csv --out OUTPUT_DIRECTORY
```

For **previously normalized, log2-transformed** gene expression:

```bash
python extended_pipeline.py --expression INPUT_LOG2_EXPRESSION.csv \
  --metadata INPUT_METADATA.csv --input-scale log2 --run-statistics \
  --condition group --contrast CategoryA:CategoryB \
  --mixed-formula 'expression ~ C(group)' --subject-column sample_id \
  --umap --overlap --overlap-term 'C(group)[T.CategoryB]' --out OUTPUT_DIRECTORY
```

To run within-region comparisons, set `--condition region`, specify desired
`--contrast RegionA:RegionB`, and optionally `--paired` to match repeated
observations by biological sample. `--stratum column=value` restricts comparisons
to one metadata level. `--nested-column roi_id` adds a variance component to
mixed models. Composite nested identifiers should be globally unique.

**Statistical boundaries:** The contrast engine computes independent-sample
Welch or paired t-tests on subject-level means. Mixed models use Gaussian
likelihood on normalized log expression with a subject random intercept;
rank-deficient or nonconvergent fits are explicitly reported, not interpreted
as significant. FDR is adjusted separately for each mixed-model coefficient,
across features; contrasts use within-contrast FDR. Overlap of different
contrasts/terms should never be interpreted as a single matched hypothesis.
Probe aggregation on raw counts is assay-dependent; the supplied probe map
and filters must be validated independently. UMAP requires `umap-learn` and
never silently substitutes PCA. These analyses do not replace instrument-specific
probe QC, distributional count models, pathology review, or independent biological
validation. Sample size, biological design, normalization and covariates must be
checked before interpreting results.

Run extension tests with `python -m unittest discover -s tests -v`.
