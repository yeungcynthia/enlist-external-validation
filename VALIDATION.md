# Verification

Verified 13 November 2026 with synthetic records only.

- 16 unit tests passed. They cover uncertainty boundaries, fixed mode, null semantics, chronological alignment, missed/extra episodes, cohort scope, original ID/date/flag parsing, the original regimen/label episode filter, absent-column scoring, CSV mapping and structure, and reproducible bootstrap sampling.
- 400 synthetic sequence alignments matched the public source and another 400 matched the paper source exactly.
- All extracted function bodies matched their source ASTs: the shared public functions and the paper uncertainty function `tol_for_u`. The restored input functions are also checked against the public source.
- Ten core line-recovery, field-agreement, and conditional PFS point metrics matched the original paper helper functions on the synthetic shared-patient cohort.
- The CLI example completed with the default 2,000 bootstrap iterations and produced the documented episode counts.

These checks verify the extraction and adapter. Performance on Stanford lung annotations has not yet been measured. Adapt the local schema and review regimen normalization before the clinical validation run.
