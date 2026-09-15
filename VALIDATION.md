# Operational layout validation

- Restored the missing runner and execution examples.
- Confirmed all five stage paths resolve under R/.
- Loaded the actual manual-flag and soil metadata from meta_files/ while running
  from another working directory.
- Checked the default six-site selection and all 30 site/stage combinations with
  deterministic test calculations. Every L2 target received all six L1 references.
- Checked no-new-data reruns, partial history merging, explicit replacement,
  locking, and failure recovery using the existing synthetic regression tests.
- Executed the actual Level 3 EC script with 1,488 synthetic half-hour observations,
  including QAQC, diagnostic generation, output publication and a skipped rerun.

The test environment uses R 4.4.2. Package build-version and synthetic regression
warnings may be emitted. Full scientific Level 1 through Level 4 validation on
real observations is not claimed: access to the Z: data root was denied in the
execution environment. No operational network outputs were modified.

The five scientific stage files match the supplied archive except for the two
metadata file-path corrections in Level 1 and Level 2. Shared I/O and the supplied
metadata are unchanged. The runner was restored from the previous delivered
version and adapted to this archive's directory layout and all-site default.
