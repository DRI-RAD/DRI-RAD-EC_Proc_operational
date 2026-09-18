# Level 5: AmeriFlux preparation

`run_level5_ameriflux.R` reads existing cumulative Level 4 files and prepares a
local review CSV. It does not rerun Levels 1-4, change observations, or upload
anything. All six operational sites are selected by default. It is a separate
export step, not an extra calculation in the normal `run_pipeline()` sequence.

## Run

Use the same R environment and installed packages as the main pipeline. The
exporter itself uses base R, readr and dplyr through the shared I/O helpers.

```r
source("run_level5_ameriflux.R")
data_root <- "Z:/NWI/Task_3_Water_Use_ET_and_Meteoroligical_Monitoring/networks/eddy_stations"

# Inspect the mapping and counts for one site without writing files.
preview <- run_level5_ameriflux(
  sites = "ECSM", base_dir = data_root, dry_run = TRUE
)
preview$ECSM$audit
preview$ECSM$summary

# Export every site's saved Level 4 period.
files <- run_level5_ameriflux(base_dir = data_root)

# Optional period: bounds are inclusive interval-END clock readings.
# Local standard time was confirmed by the PI; there is no UTC/DST conversion.
files <- run_level5_ameriflux(
  sites = "ECSM", base_dir = data_root,
  start = "2025-01-01 00:30:00", end = "2026-01-01 00:00:00"
)
```

With no bounds, this exports the complete saved Level 4 time span. It does not
re-read EddyPro to set an end. Missing internal half-hours are inserted with
missing measurements. Timestamps are sorted and must be unique half-hour ends.
Duplicate or off-grid timestamps stop export. Calendar leap days are preserved.

## Files

Each run writes a new directory; older exports and all Level 4 files remain:

```text
<site output>/AmeriFlux/<execution timestamp>_<internal site>/
  ECSM_Level_5_ameriflux_ready.csv
  README.txt
  metadata/
    variable_mapping.csv
    coverage.csv
    source_files.csv
```

An optional `output_dir` redirects these exports to `<output_dir>/<site>/...`.
The data file has one header, no units row or quotes, and two leading clock
columns: `TIMESTAMP_START` and `TIMESTAMP_END` in `YYYYMMDDHHMM` format.
The start is exactly 30 minutes before the existing end. Missing/nonfinite
measurements are `-9999`; valid timestamps are retained for missing intervals.

`variable_mapping.csv` records source/target names, units and conversion notes,
including omitted variables. `coverage.csv` counts observations and gaps in each
exported column. `source_files.csv` records source paths and MD5 checksums.
Metadata files are for review and should not be submitted as flux data.

## PI-selected flux policy

| Output | Level 4 source and handling |
| --- | --- |
| LE, H, FC | Retained turbulent LE/H/FC masked only where the respective `_QC_despike >= 1` or `_QC_fetch >= 1` |
| LE/H/FC_SSITC_TEST | Supplied separately; NOT applied to the turbulent flux |
| SLE, SH, SC | Separate Level 4 storage estimates; not added to base LE/H/FC |
| G | G_plate, the averaged plate-depth soil heat flux |
| SG | Soil storage above plates, provided separately |
| G_PI | Calorimetric surface heat flux, already including soil storage |
| LE_PI_F, H_PI_F, FC_PI_F | Existing processed and filled Level 4 series, unchanged |
| LE_PI_CORR, H_PI_CORR | Existing EBR-corrected Level 4 series, unchanged |
| NETRAD_PI_F, G_PI_F | Existing filled radiation/surface soil heat-flux estimates |

Physical-range rejection is already part of the despike flag. SSITC, USTAR,
signal-strength and long-run flags do not filter the base Level 5 LE/H/FC.
No new despiking is calculated. Missing QC flags do not reject a value, matching
the pipeline's unassessed-despike behavior; their counts appear in mapping notes.
This is not equivalent to selecting Level 4 `*_filtered1` or `*_filtered2`.

The PI flux series retain Level 4's more extensive filtering, storage handling,
and MDS estimates. In particular, `FC_PI_F` is an NEE-like series: its storage
term is included when available, with the existing turbulent-flux fallback when
storage is absent. It must not be mistaken for base turbulent FC. Likewise,
the PI LE/H series may include storage. Do not add SG again to G_PI or G_PI_F.
SLE/SH/SC retain any Level 4 storage despiking already performed.

## Other measurements and conversions

The automatic mapping includes available radiation, precipitation P, pressure,
TA/RH/TS/SWC sensors already named with standard numeric position suffixes,
wind statistics, gas mole fractions and standard deviations, sonic temperature,
footprint distances, and SSITC tests. Existing filled TA/RH/PA/VPD/SW_IN/LW_IN/P
series are also retained under their `_PI_F` names. Sensor position suffixes
are preserved; this code cannot verify their heights/depths or BADM definitions.

Unit handling is explicit: pressure is converted to kPa, VPD to hPa, temperature
to deg C, fractions to percent when the source specifies a fraction, and fluxes
to W m-2 or umolCO2 m-2 s-1. Temperature standard deviations have no offset.
Unknown conversions omit the optional variable and record the reason; an unknown
LE/H/FC unit stops export. No unit is inferred from the magnitude of observations.

Two source-code-specific cases require more than trusting a header:

- L3 mixes EddyPro sonic temperature in K with EasyFlux values while retaining
  an EasyFlux header. `processing` determines which rows need the K-to-deg-C
  conversion; EasyFlux rows use the source unit. Unknown provenance with finite
  observations causes that optional variable to be omitted.
- L1 explicitly recomputes ALB as percent,
  so that calculated unit takes precedence over a legacy logger header.

For a reviewed metadata error, override the SOURCE unit (not the desired unit):

```r
preview <- run_level5_ameriflux(
  sites = "ECSM", base_dir = data_root, dry_run = TRUE,
  unit_overrides = c(PA = "kPa")
)
```

Use overrides only after confirming instrument/calculation units. Optional
`extra_mapping` is a data frame with `source`, `target`, `unit` (target unit), and
`note` columns. It can map additional sensors after their identities and positions
are confirmed. It is not a substitute for checking official variable names.
Automatic export deliberately does not invent positions for CS65X/Teros/LI-710
channels or standardized names for VBR/RE estimates. All unmapped columns appear
in the audit, so they can be reviewed and added explicitly.

## Before actual submission

The PI requested internal filenames and `_PI` products for local review. The
`ameriflux_ready` filename therefore does NOT certify submission compliance.
Registered site IDs are not needed to make this local file. Before uploading:

- Use the registered site ID and official filename convention.
- Resolve the local `_PI` labels: the upload instructions reserve `_PI` for the
  network; submitted gap-filled variables use `_F`. Do not blindly relabel a
  storage-inclusive PI flux as a turbulent flux. FC versus NEE and storage
  conventions need particular care.
- Review omitted variables, unit overrides, sign conventions, the sensor-position
  mapping, and the matching site/timezone and variable-processing metadata.

The exporter has synthetic tests for the requested masks, storage separation,
PI preservation, unit conversions, missing intervals and leap day, errors and
CSV round-trip validation. It has not run against operational Z: observations
in this environment because that path was inaccessible. It has not been validated
by the AmeriFlux submission service.

## References

- [AmeriFlux variable definitions](https://ameriflux.lbl.gov/data/aboutdata/data-variables/)
- [Half-hourly/hourly upload instructions](https://ameriflux.lbl.gov/data/uploading-half-hourly-hourly-data/)
- [Submission tips, including _PI and USTAR guidance](https://ameriflux.lbl.gov/data/tips-ameriflux-flux-met-data-processing-pipeline/)
- [EddyPro full-output units](https://bio.licor.com/env/support/EddyPro/topics/output-files-full-output.html)

Reviewed September 18, 2026. Website submission requirements take precedence
over the local review layout. Run `Rscript tests/test_ameriflux.R` from this
project's root to repeat the synthetic checks.
