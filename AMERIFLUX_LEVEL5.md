# Level 5: AmeriFlux observation export

`run_pipeline()` automatically exports each site's complete saved Level 4 period
following its L4 step, including when there are no pending L4 observations.
`dry_run = TRUE` writes nothing; runs without L4 do not export.
`export_ameriflux = FALSE` disables the automatic step. Export paths are returned
in `attr(result, "ameriflux")`. Export failures retain committed L4 files.

## Default: base observations and explicitly labeled filled products

The default exporter uses `include_pi = FALSE`. Base SW_IN, TA/RH/TS/SWC sensor
columns, PA and radiation come from L1 QC, BEFORE L2 gap filling. Base precipitation
comes from L1 P_Bulk and P, while P_F_1_1_1 is the separate L2 filled series.
Filled LI-710 LE/H are also exported under their explicit `_F` names; these
never replace the unfilled LE/H columns. Existing EC LE_PI_F, H_PI_F and FC_PI_F
are also included by default as local review columns. `include_pi` controls
additional PI products; it does not suppress these three EC MDS products.
L5 never fills measurement gaps; missing internal clock intervals get missing
measurements. Missing/nonfinite values are written as -9999.

Derived soil storage is an exception to a strictly observation-only description:
SG uses L2 calculations with interpolated temperature derivatives and median
SWC substitution in the heat-capacity calculation. G and SG are spatially
weighted derived products, not individual plate observations.

## EC and LI-710 naming and QC

| Export variable | Level 4 source | Processing |
| --- | --- | --- |
| LE_1_1_1, H_1_1_1 | LE, H | L3 EddyPro/EasyFlux selection; L5 applies despike/physical-range and fetch flags only. No SSITC filter, storage addition, MDS or EBR correction |
| FC | FC | Same EC mask policy; no storage addition or filling |
| LE_PI_F, H_PI_F | Same-name columns | Existing L4 QC, storage addition where available (turbulent fallback otherwise), then MDS filling; no EBR correction. Local review names, not approved upload labels |
| FC_PI_F | FC_PI_F | Existing L4 storage-inclusive CO2 series (turbulent fallback where storage is missing), filled as NEE by MDS. Local review name; not a filled turbulent FC series |
| LE_1_1_2, H_1_1_2 | LE_710_filtered, H_710_filtered | Unfilled LI-710 series: L3 diagnostic, despike/physical-range, long-run and flow flags applied in L4. No MDS or EBR correction |
| LE_F_1_1_2 | LE_710_PI_F | Existing L4 MDS-filled LI-710 LE after QC; no storage or EBR correction |
| H_F_1_1_2 | H_710_PI_F | Existing L4 MDS-filled LI-710 H after QC; no storage or EBR correction |
| LE_SSITC_TEST_1_1_1, H_SSITC_TEST_1_1_1 | LE_SSITC_TEST, H_SSITC_TEST | L3 selected EC SSITC flags provided separately |
| FC_SSITC_TEST, TAU_SSITC_TEST | Same-name columns | L3 selected EC flags |
| SLE_1_1_1, SH_1_1_1, SC | SLE, SH, SC | Separate EddyPro storage retained after L4 storage despiking; not added to base EC flux |

The requested `_1_1_1` / `_1_1_2` mapping distinguishes EC and LI-710 locally.
The three indices describe horizontal position, vertical position and replicate,
not instrument models. Confirm their actual locations/heights and corresponding
BADM records before submission. If the systems have different heights or
positions, adjust the mapping accordingly. No LI-710 SSITC flag is invented.
If the LI-710 source column is absent, its export column is absent; if present
but entirely missing, it remains an all--9999 column. Filled LI-710 values are
never substituted into unfilled columns. Base EC can use the existing L3 EasyFlux fallback where no
EddyPro fields are available; this is not MDS gap filling.

## Other default variables

All immediate inputs are the cumulative L4 CSV. Optional columns are included
only when present and their units can be converted. `variable_mapping.csv`
records the exact mapping, omissions and unit-conversion reasons for each run.

| Export variable(s) | Unit | Upstream origin and processing |
| --- | --- | --- |
| TIMESTAMP_START, TIMESTAMP_END | YYYYMMDDHHMM | L5: local-standard interval end, start = end minus 30 minutes |
| NETRAD, SW_IN, SW_OUT, LW_IN, LW_OUT | W m-2 | L1 QC; NETRAD recomputed from the four screened radiation components; no L2 filling |
| PPFD_IN, ALB | umolPhoton m-2 s-1; % | L1 PPFD screening; ALB = SW_OUT/SW_IN*100 (zero when SW_IN = 0) |
| PA | kPa | L1 pressure QC; no L2 filling |
| P_1_1_1 | mm per interval | L1 P_Bulk, from Pluvio or corrected SG400; no L2 gap filling |
| P_1_1_2 | mm per interval | L1 P, tipping-bucket precipitation; no L2 gap filling |
| P_F_1_1_1 | mm per interval | L2 P_PI_F: selected precipitation source, regression filling and final zero fallback. May use tipping-bucket source where bulk is unavailable |
| TA, RH, TS, SWC and numeric _i_j_k variants | deg C; RH/SWC: % | L1 variable-specific QC; no L2 filling |
| G | W m-2 | L2 G_plate: sun/shade-weighted plate-depth flux from L1-screened plates |
| SG | W m-2 | L2 weighted storage above plates: w*SG_sun + (1-w)*SG_shade, w from meta_G.csv G_sun_weight; includes interpolated temperature derivatives and missing-SWC substitution |
| WS, WS_MAX, USTAR, U_SIGMA, V_SIGMA, W_SIGMA | m s-1 | L3 selected EC turbulence statistics; no MDS or additional L5 flux mask |
| WD, ZL, MO_LENGTH, TAU | degrees; 1; m; kg m-1 s-2 | L3 selected EC output; no MDS or additional L5 flux mask |
| CO2, CO2_SIGMA, H2O, H2O_SIGMA | umolCO2 mol-1; mmolH2O mol-1 | L3 concentrations and sigmas; EddyPro sigmas calculated from variances and concentration/density conversion |
| T_SONIC, T_SONIC_SIGMA | deg C | L3 sonic statistics; L5 converts EddyPro temperature from K by row provenance; sigma has no temperature offset |
| FETCH_90, FETCH_70 | m | L3 footprint estimates; EasyFlux has no FETCH_70 mapping |

Sensor suffixes are preserved. Nonstandard CS65X/Teros names, VBR, RE, ET-depth,
MDS prediction and internal QC/method columns are not mapped automatically.
Unknown optional units cause an audited omission; unknown base EC flux units
stop export. Missing EC despike/fetch flags do not reject a value and their
counts are recorded. Other variables do not inherit the EC flux masks.

## Optional PI products for local review

The PI chose to preserve existing calculations and include EC filled exports
under the local names LE_PI_F, H_PI_F and FC_PI_F by default. These are the
existing MDS filled series, not quality-only predictions. No LE_F_1_1_1,
H_F_1_1_1 or FC_F_1_1_1 is generated: these results can include storage and are
not simply filled versions of the exported turbulent flux. No new MDS run,
storage subtraction or additional filtering is applied. Do not add SLE/SH/SC
to these products again. `include_storage` controls separate storage columns only.
LI-710 does not measure CO2 flux, so no LI-710 FC_F column is invented.

EBR-corrected CORR products remain excluded by default. Sensor indices must not
be repurposed to conceal a correction method. The official definitions and
qualifier examples were checked in a browser on October 22, 2026. They do not
define a tower-upload CORR qualifier. Upload instructions section 4.2 allow
nonstandard variables alongside standard variables, with the Data Team contacted
to agree on names and units; this does not guarantee publication under an
arbitrary PI_CORR label. For PI energy-balance closure products, provide the
correction method and storage/gap-filling provenance to the Data Team and agree
on the submission mapping. The existing LE_PI_CORR and H_PI_CORR are available
via `include_pi = TRUE` for local review; they remain distinct from EC PI_F.

`run_level5_ameriflux(..., include_pi = TRUE)` additionally retains G_PI and
existing LE_PI_CORR, H_PI_CORR, NETRAD_PI_F, G_PI_F,
TA_PI_F, RH_PI_F, PA_PI_F, VPD_PI_F, SW_IN_PI_F, LW_IN_PI_F and P_PI_F.
This is a local review option, not a claim that these names are upload-ready.

L2 fills the meteorological PI series by variable-specific interpolation,
regression and mean-diurnal methods; VPD derives from filled TA/RH. P_PI_F sets
remaining missing precipitation to zero. L4 MDS fills PI fluxes, NETRAD and G_PI;
PI fluxes can include storage and additional QC, and PI_CORR also includes EBR
correction. FC_PI_F is storage-inclusive where storage is available. The current
L3 USTAR flags are all zero, so the L4 FC USTAR condition does not actively reject
observations. G_PI already includes SG: never add SG to it again.

The documented variable vocabulary reserves `_PI` for the network team; `_F`
is available for tower-team gap-filled submissions. Combining local qualifiers
such as `_PI_F` or `_PI_CORR` does not make them accepted submission names.
Do not simply rename a storage-inclusive PI flux as turbulent LE/H/FC.
The default local export now contains the three requested EC PI_F products;
their inclusion does not certify these column names for upload. Upload instructions
section 5 accept only `_F` as a general qualifier. The network-only `_PI` label
does not by itself authorize energy-balance closure or storage corrections.

## Independent export and files

```r
source("run_level5_ameriflux.R")
# Complete L4 data: base products, filled LI-710/P, and local EC PI_F products.
files <- run_level5_ameriflux(base_dir = data_root)
# Optional local review of PI products:
review <- run_level5_ameriflux(base_dir = data_root, include_pi = TRUE)
# Preview values and mapping without writing:
preview <- run_level5_ameriflux(sites = "ECSM", base_dir = data_root, dry_run = TRUE)
preview$ECSM$audit
```

Each execution creates a new folder; prior exports and L4 data are retained:

```text
<site output>/AmeriFlux/<execution timestamp>_<site>/
  <site>_Level_5_ameriflux_ready.csv
  README.txt
  metadata/variable_mapping.csv
  metadata/coverage.csv
  metadata/source_files.csv
```

The data CSV has one header and no unit row. Metadata contain mapping, finite/
missing counts, input paths and MD5 hashes. Optional start/end bounds select
interval-end timestamps for standalone exports. `unit_overrides` changes the
assumed SOURCE unit only after review; `extra_mapping` explicitly adds reviewed
source/target/unit/note mappings. Defaults export seven operational sites.

Internal site filenames still require registered AmeriFlux IDs and the applicable
submission filename convention. The local `ameriflux_ready` name is not a
certification by the submission service. Official pages were verified through
the browser on October 22, 2026; the upload service itself was not tested.

## EddyPro random uncertainty in L3/L4 (separate from this export)

L3 retains `LE_rand_err`, `H_rand_err`, `FC_rand_err`, `TAU_rand_err` from EddyPro
`rand_err_LE`, `rand_err_H`, `rand_err_co2_flux`, `rand_err_Tau`. L4 carries these
columns and units unchanged. They describe the selected unfilled EddyPro
flux, not an uncertainty propagated through storage, MDS or EBR correction.
EasyFlux fallback rows and missing optional uncertainty inputs remain NA.
These columns are not automatically mapped into the AmeriFlux export.

The reader supports only this four-column additive schema change for old L3/L4
files, padding unknown historical uncertainty with NA. Other schema differences
still stop processing. Previously saved files are not rewritten on read. To
populate uncertainty over an existing period, reprocess L3_EC and L4 (or all
stages) with `reprocess = TRUE`; an ordinary incremental run does not backfill it.

## Script references

- L1: R/Level_1_met_QAQC_new.R
- L2: R/Level_2_met_PI_vars_new.R
- L3 EC: R/Level_3_EC_data_QAQC_new.R
- L3 LI-710: R/Level_3_LI710_QAQC_new.R
- L4: R/Level_4_gapfilling_ET_calculation_new.R
- L5: R/Level_5_AmeriFlux_export.R

[AmeriFlux variable definitions](https://ameriflux.lbl.gov/data/aboutdata/data-variables/)
[Half-hourly upload instructions](https://ameriflux.lbl.gov/data/uploading-half-hourly-hourly-data/)

Tests cover unfilled EC/LI-710 mapping, PI opt-in, unit conversions, automatic
exports, uncertainty propagation and compatible history merging. Production
Z: data and the AmeriFlux upload service were not exercised in this environment.
