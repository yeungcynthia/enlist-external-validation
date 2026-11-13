# Source and adaptation record

Prepared 13 November 2026 for external LoT/PFS validation.

## Alignment source

Public upstream: [elenahuseni/EnLiST_LoT_PFS_Pipeline](https://github.com/elenahuseni/EnLiST_LoT_PFS_Pipeline), revision `a7dd647dc100451601472c7c327cbf2c25e2449c`, `transforms-python/src/myproject/nw_align.py`.

The alignment functions, scoring constants, drug dictionaries, input parsers, and selected field comparators in `alignment.py` were extracted from the paper evaluation checkout at revision `6e6d10785fb30a5e30c01c3ab6337d1c3f982b21`. Its alignment/scoring/normalization functions are identical to those in the public revision. The paper checkout adds `U_TOLERANCE`, `U_TOLERANCE_DEFAULT`, and `tol_for_u`, which are also retained verbatim. Unneeded field normalizers and Foundry adapters were omitted. The CSV adapter uses the original regimen/LoT-label episode filter. This file is a focused extraction, not a complete copy of `nw_align.py`.

Source full-file SHA-256 hashes:

- Public `nw_align.py`: `0227447821ae98880e15cd3d27314ccb89252aaeead474fb8256e4a75a3b51dc`.
- Paper checkout `nw_align.py`: `9a408ad8a80590d6c5543122c15d4f44c03b50e6e5044ddab6c48e2c2094cb0b`.

## Metric sources

The standalone wrapper adapts the metric definitions from these files in the same paper checkout:

- `datasets/stage4_eval_vs_gs.py`: per-field comparators on NW-aligned pairs; line sensitivity, PPV, and F1; gold-side U-code date tolerances.
- `datasets/paper_cluster1_accuracy.py`: episode recovery, PFS flag agreement, directional event/censor discordance, patient-clustered bootstrap (2,000 iterations, seed 20260904).
- `datasets/paper_pfs_conditional.py`: conditional progression/censor date subsets and patient-clustered percentile intervals.

The paper checkout source was identified from the existing analysis workspace and the analysis folder supplied by the project owner. This is a snapshot, not a claim that all Foundry branches have been synchronized to it.

## Deliberate adapter differences

The wrapper uses CSV instead of Foundry/Spark. It includes only regimen, LoT label, start/stop, progression/censor dates, and event/censor flag. It has no extraction prompts, figures, manuscript material, internal dataset identifiers, development-patient exclusions, clinical records, credentials, or model calls.

Following the owner's requested correction, `normalize_mrn`, `parse_date`, and `norm_flag` are retained verbatim from the public source. The adapter uses their original behavior, compares all predefined key fields even when columns are absent, and skips rows with neither regimen nor LoT label, exactly as the source episode filter does. It keeps CSV column mapping and structural checks. Fixed date tolerances are now the default for the command line and Python alignment interface, matching the public pipeline: start/progression ±14 days and stop/censor ±30 days. Uncertainty mode remains an explicit opt-in. Default cohort scope includes every gold patient so absent predictions count as misses; `--patient-scope common` restores the original intersection rule. F1 is zero for a nonempty one-sided timeline and unavailable only if both episode totals are zero. Additional gold-date-present metrics and known-flag agreement separate missingness from populated-value agreement. No lung-specific synonyms were added. Raw confusion counts do not receive bootstrap intervals. The cohort, including reviewed zero-episode patients, supplies the bootstrap sampling frame; source helpers sampled patients represented in detail rows. These choices are described in the README.

The MIT notice from the public upstream is retained in LICENSE.

Event-status sensitivity, specificity, PPV, NPV, and F1 are additional standard calculations from the source-aligned binary confusion counts, conditional on aligned episodes. The wrapper also reports excluded unknown/nonbinary status pairs. These formulas and their bootstrap intervals were added for the collaborator aggregate-results handoff; they are not claimed as verbatim source functions.
