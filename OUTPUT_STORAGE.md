# Output storage and period diagnostics

One current cumulative CSV per stage lives directly in the configured site
output directory. Each update merges pending rows with prior observations.
The filename date is the local execution date, not the observation period.

```text
Columbus_Salt_Marsh/
  ECSM_Level_4_post_processed_data_2026-09-16.csv
  ECSM_Level_4_post_processed_data_2026-09-16.csv.period.csv
  archive/
    ... every previous Level 4 CSV and its audit ...
    ... only the previous version of each Level 1-3 product ...
  figures/
    ECSM_Level_4_Figure_1a_energy_balance_full_period.jpeg
    ... the other four cumulative Level 4 figures ...
    period_20260716T1700_to_20260801T0000/
      ECSM_Level_1_Figure_met_QAQC.pdf
      ECSM_Level_3_LE_flag.pdf
      ECSM_Level_3_LE_710_flag.pdf
      ECSM_Level_4_Figure_1a_energy_balance.jpeg
      ... remaining named figures for this period ...
    ... older files/directories left untouched ...
```

## Rotation and figure updates

- Rotation is per successful stage. All original L4 CSVs and accompanying audits
  move to archive; no L4 archive version expires. Same-day filename collisions
  receive a time/unique suffix. Multiple legacy L4 root versions are preserved
  individually. Levels 1-3 retain their previous-version-only policy.
- With `end = NULL`, period directories use one common envelope per site/run:
  the earliest pending timestamp across selected stages through that site's last
  EddyPro observation. Stages share this folder even if their checkpoints differ.
  With an explicit `end`, directories retain each stage's first/last pending
  timestamps. Interior holes are possible; folder bounds do not imply that every
  interval was updated. Each figure still uses only its stage's pending rows.
  Reprocessing the same period overwrites only matching figure names there.
- All five cumulative L4 figures use the entire merged output, including prior
  rows outside the update. Their `_full_period` suffix is reserved for automatic
  replacement; it avoids overwriting existing period figures in the root.
- Cumulative figures are rendered from already processed values, without rerunning
  historical QAQC/MDS. Missing LI710 values remain missing and its panels can be
  empty. Period figures exclude training context; `figure_period = "new"` remains
  the supported argument.
- No automatic plot PDFs or permanent staging input copies are generated.
  Existing figures, figure archives, and legacy run folders are not moved,
  deleted or migrated. They can be organized manually.
- Skipped stages do not rotate CSVs or regenerate figures. Corrected/late source
  observations at processed timestamps require `reprocess = TRUE`. Bounds are
  optional: omitted start uses the earliest available interval and omitted end
  uses the site's latest EddyPro observation. Rows outside that range are retained.

The new CSV is prepared before publication. Figure generation/copy failures
roll back current/archived CSV changes and overwritten figures. Forced process
termination can leave a hidden `.publish-*` recovery directory, which blocks
another publication. Preserve its `recovery.rds` and `installed.rds` files until
an operator has restored the backup or confirmed publication completed. The site
lock must also be cleared only after recovery and after confirming no run is active.

Automatic-end processing caps new rows at the site's EddyPro end and retains NA
observations/QC for absent LI-710 logger intervals. Existing cumulative rows beyond
that end are not deleted; cumulative Level 4 figures continue to show the full
saved history. The period audit's pending coverage describes actual stage rows,
whereas its figure period label describes the shared run envelope.

## Period audit and quality statistics

The `*.csv.period.csv` file remains a single-row audit. It includes total output
coverage, pending coverage/counts, calculation start/duration, reprocessing mode,
and the figure period label. All quality statistics describe **pending rows**,
excluding historical MDS context and the rest of the cumulative dataset.

For LE, H, FC, LE_710, and H_710, available columns generate these metrics:

| Metric suffix | Numerator | Denominator |
| --- | --- | --- |
| `raw_observed` | Finite incoming flux values | Pending rows |
| `QC_<test>_code_<value>` | Rows with that exact QC code | Rows with a nonmissing code for that test |
| `QC_<test>_missing` | Missing QC codes | Pending rows |
| `qc_removed` | Observed raw flux that is missing/nonfinite after the actual final QC filter | Observed raw flux rows |
| `gapfilled` | Missing filtered flux with a finite filled result | Pending rows |
| `gapfill_success` | Missing filtered flux with a finite filled result | Missing filtered flux rows |
| `remaining_gap` | Missing/nonfinite filled result | Pending rows |
| `good_gapfilled` | Actual filled gaps with REddyProc `_fqc == 1` | Successfully filled gap rows |

Each metric has `_n`, `_denominator_n`, and `_pct` columns. Percentages range
from 0 to 100; a zero denominator yields NA. For example,
`LE_qc_removed_pct = 20` means 20% of observed incoming LE was actually removed.
`LE_good_gapfilled_pct = 80` means 80% of successfully filled LE gaps have MDS
quality code 1. Predictions at observed timestamps from `FillAll = TRUE` are
excluded from filling statistics.

L3 supplies incoming flux and QC-code metrics; L4 additionally supplies actual
removal and filling metrics. Earlier stages with no matching flux columns record
the period/count metadata only. Missing metrics are not assumed to be zero.
QC tests can overlap, so their code percentages must not be added to estimate
removal. The actual filtered columns determine removal instead.

"Raw" here means the incoming flux retained in the L3/L4 table, including the
pipeline's existing source-selection/fallback behavior. It is not a measure of
EddyPro-only file completeness. No QC thresholds or gap-filling algorithms are
changed by collecting these statistics.
