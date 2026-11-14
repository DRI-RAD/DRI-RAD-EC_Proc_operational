# Level 5: curated AmeriFlux submission export

The exporter uses an explicit allowlist of AmeriFlux physical quantities plus the PI-requested final tower-team products. It does not export every Level 4 column. The PI reports Data Team agreement to include filled, EBC-corrected and LI-710 products with comments. At the PI's request, this exporter explicitly labels tower estimates with `_PI`; `_PI` and custom `_CORR` naming should be reviewed with the Data Team (the published qualifier guide reserves `_PI` for network use).

`run_pipeline()` exports the complete saved L4 after the L4 stage, even when no new rows are pending. `export_ameriflux = FALSE` skips this. To regenerate only L5:

```r
source("run_level5_ameriflux.R")
files <- run_level5_ameriflux(base_dir = data_root)
```

| Export | Source / processing |
|---|---|
| LE_1_1_1, H_1_1_1, FC | Raw EC, L3 despike/physical-range and fetch masks only; no SSITC mask, storage or filling |
| LE_1_1_2, H_1_1_2 | LE_710_filtered, H_710_filtered; L3 diagnostic, despike, long-run and flow QC applied in L4 |
| LE_PI_F_1_1_1, H_PI_F_1_1_1, FC_PI_F | Final L4 filled EC products; includes available storage, so FC_PI_F includes SC where available |
| LE_PI_F_1_1_2, H_PI_F_1_1_2 | Final L4 filled LI-710 products; no EC storage addition |
| LE_PI_CORR_1_1_1, H_PI_CORR_1_1_1 | Final EC filled and EBR-corrected products |
| LE_PI_CORR_1_1_2, H_PI_CORR_1_1_2 | Final LI-710 filled and EBR-corrected products |
| SLE_1_1_1, SH_1_1_1, SC | Separate storage estimates; never add again to filled fluxes |
| LE_SSITC_TEST_1_1_1, H_SSITC_TEST_1_1_1, FC_SSITC_TEST, TAU_SSITC_TEST | Standard SSITC tests, separate from flux masks |
| G_1_1_1, G_1_1_2, ... | Individual L1 QAQC heat-flux plates; source suffix retained |
| G_PI | L2 G_PI before gap filling; selected surface estimate including above-plate storage |
| SG_PI | L2 weighted SG; do not add again to G_PI |
| TS_1_1_4 through TS_1_1_6 | TS_CS65X_1_1_1 through TS_CS65X_1_1_3; CS65x source recorded in metadata |
| SECA_1_1_1 through SECA_1_1_3 | cs65x_ec_1_1_1 through cs65x_ec_1_1_3; apparent soil conductivity, dS m-1 |
| WS_1_1_1, WD_1_1_1, WS_MAX_1_1_1 | EC wind quantities |
| WS_1_1_2, WD_1_1_2, WD_SIGMA_1_1_2, WS_MAX_1_1_2 | RM Young 05103 observations; maximum is the maximum five-second wind speed |
| TA_1_1_4, RH_1_1_4, PA_1_1_2 | Retained LI-710 meteorology; no filling |
| PA_1_1_1 | Retained primary pressure |
| P_1_1_1, P_1_1_2 | Bulk and tipping-bucket precipitation before filling |
| P_PI_F_1_1_1 | Final selected/filled precipitation, including zero fallback |
| TA_PI_F, RH_PI_F, PA_PI_F, VPD_PI_F, SW_IN_PI_F, LW_IN_PI_F | L2 filled/derived meteorology |
| NETRAD_PI_F | L4 MDS-filled net radiation |

Other explicitly mapped radiation, indexed TA/RH/TS/SWC, turbulence, gas concentrations and footprint variables are retained. Units and source scripts for **every exported column** are in `metadata/variable_mapping.csv`. Fixed sensor slots do not shift when another sensor is absent; conflicts fail explicitly. Sensor positional indices must be confirmed against site sensor metadata; the PI has confirmed EC and LI-710 flux sensors are collocated at the same height.

Excluded: sun/shade and calorimetric/Leuning method components, weighted plate mean, G_PI_F, raw LI-710 duplicate fluxes, filtered/storage-added intermediate series, internal despike/flow/long-run QC and filling-method codes, MDS prediction columns, VBR/RE/ET products, and instrument diagnostics. No verified standard mapping was found for the retained minimum gas signal-strength channels, so these are omitted. None of these exclusions deletes L4 data or changes science calculations.

Each run creates a new export folder containing the data CSV, `submission_comments.txt`, `README.txt`, metadata (source/target mapping and omissions, coverage, source-file checksums), and optional diagnostics. Missing data are -9999. Timestamps are local standard time, half-hour interval starts/ends. No data are uploaded. `include_pi = FALSE` suppresses optional PI estimates but retains final EC/LI-710 filled products for backward compatibility. Explicit `extra_mapping` remains available for reviewed exceptions; no automatic numeric-column passthrough occurs.

References: https://ameriflux.lbl.gov/data/aboutdata/data-variables/ and https://ameriflux.lbl.gov/data/data-variable-qualifier-examples/

## Automatic submission diagnostics

By default `run_level5_ameriflux(..., make_plots = TRUE)` generates diagnostics
after the submission CSV is safely saved. The automatic L5 step in
`run_pipeline()` also enables them; use `ameriflux_plots = FALSE` to skip just
these diagnostics while retaining the CSV export and existing L4 figures.
`dry_run = TRUE` writes neither data nor plots.

Each export folder contains:

```text
diagnostics/
  index.html
  variable_summary.csv
  figures/
    001_LE_1_1_1.png
    ...one PNG for every submitted data column...
```

Open `diagnostics/index.html` in a browser. The report works offline with no
CDN or external libraries. Search by variable, original source, unit or processing
description; filter all-missing/partially-missing/complete variables; expand
individual plots or all visible results. Click a plot to open the full-resolution
PNG. Keep the HTML and figures folder together when sharing the report.

Each PNG contains the actual exported half-hour values (all points, full range),
a histogram or discrete QC/method-code frequency chart, and monthly valid-data
percentages. No new filling, filtering, clipping or thinning is applied.
Monthly completeness uses intervals represented in the export, including partial
months; midnight interval ends belong to the preceding month. It is not a claim
of complete calendar-month coverage. Timestamp columns provide the time axis,
not separate numeric plots. Entirely missing columns still receive a placeholder
figure. Units, source, processing notes, missing percentage, min/max/mean/median,
and first/last valid timestamps are saved in the report/summary.

```r
# Default: export CSV, comments, metadata, PNGs and HTML.
files <- run_level5_ameriflux(base_dir = data_root)
attr(files, "diagnostics")

# Export quickly without plots.
files <- run_level5_ameriflux(base_dir = data_root, make_plots = FALSE)

# Regenerate diagnostics from an existing submission CSV, without L1-L4 or L5 export.
plot_level5_ameriflux(files[["ECDP"]])

# Pipeline: retain the automatic export but skip its plots.
run_pipeline(base_dir = data_root, ameriflux_plots = FALSE)
```

Regeneration overwrites the report and its named PNGs, never the input CSV.
Figure creation uses base R graphics, so there are no additional R packages to
install beyond the pipeline's existing readr dependency. Generating a PNG for
every variable adds runtime and disk space. Individual plot errors are recorded
in the summary/HTML; a report-level error produces a warning and
`diagnostic_error.txt` while preserving the committed CSV.
