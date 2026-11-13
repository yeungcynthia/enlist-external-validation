# Verification

Verified 13 November 2026 with synthetic records only.

- 13 unit tests passed, covering uncertainty boundaries, fixed mode, null semantics, chronological alignment with an inserted episode, low-similarity correspondence, missed patients, extra episodes, cohort scope, date-only annotations, ID preservation, column mapping, input errors, and reproducible bootstrap sampling.
- 400 synthetic sequence alignments matched the public source and another 400 matched the paper source exactly.
- All extracted function bodies matched their source ASTs: eight shared public functions and nine paper functions including `tol_for_u`.
- Ten core line-recovery, field-agreement, and conditional PFS point metrics matched the original paper helper functions on the synthetic shared-patient cohort.
- The CLI example completed with the default 2,000 bootstrap iterations and produced the documented episode counts.

These checks verify the extraction and adapter. Performance on Stanford lung annotations has not yet been measured. Adapt the local schema and review regimen normalization before the clinical validation run.
