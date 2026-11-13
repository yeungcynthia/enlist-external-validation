# EnLiST external validation

Align manually annotated and extracted treatment episodes, then benchmark LoT start/stop dates, PFS event/censor status, progression dates, and censor dates. This is a small, local Python package for collaborators adapting the [EnLiST pipeline](https://github.com/elenahuseni/EnLiST_LoT_PFS_Pipeline). No Foundry or model API access is required.

Default date scoring uses **the gold annotation's uncertainty code**, with a 30-day fallback. The alignment algorithm and regimen normalization are unchanged from the source implementation. See [PROVENANCE.md](PROVENANCE.md) for the distinction between the public pipeline and paper evaluation code.

## Run the synthetic example

Requires Python 3.10 or later.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
python evaluate.py \
  --gold example_gold.csv \
  --predicted example_predicted.csv \
  --output results/example
python -m unittest -v
```

The output directory must be new. This protects prior results from accidental replacement.

Expected example counts: 3 evaluated patients, 3 aligned pairs, 1 missed gold episode, no extra predicted episodes. Line sensitivity = 3/4, PPV = 3/3, F1 = 6/7. Start-date agreement = 2/3 under uncertainty tolerances. The example patients, dates, and regimens are synthetic.

## Adapt your annotations

Supply one CSV row per treatment episode on each side. `patient_id` must identify the same person across files. IDs are kept as strings, including leading zeros. Restrict the gold file to patients and episodes that were fully reviewed for the evaluation scope.

| Column | Meaning |
| --- | --- |
| `patient_id` | Required, local study identifier |
| `systemic_anticancer_agents` | Regimen, preferably generic drug names joined with `+` |
| `enlist_lot_label` | Optional LoT label; ordinary labels such as `1`, `2` are accepted |
| `start_date`, `stop_date` | LoT start and stop, `YYYY-MM-DD` |
| `pfs_censor_flag` | `1` = PFS event, `0` = censored, blank = unavailable |
| `first_pd_date` | Progression date, `YYYY-MM-DD` |
| `censor_date` | Censor date, `YYYY-MM-DD` |
| `start_date_U`, `stop_date_U`, `pd_date_U`, `censor_date_U` | Optional gold uncertainty code for the corresponding date |

`first_pd_date_U` is accepted if `pd_date_U` is blank or absent. Only gold uncertainty codes affect scoring. Use U0–U4 or leave blank. Do not invent codes for annotations that did not capture uncertainty.

Only fields with a column on **both** sides are scored. Omit a column when that field was not collected. An empty cell in a present column is treated as a missing/absent value under the source comparator, including agreement when both values are absent. For dates, separate gold-date-present and both-date-present metrics reveal how much agreement depends on absent dates. If blanks mix unreviewed and confirmed absent values within a column, prepare a consistently reviewed subset before scoring.

A patient row with all evaluable fields blank represents zero episodes. This allows recording a completed review with no episodes. A row containing only a date is retained as an episode. Stable chronological sorting uses start date, with missing starts last and input order breaking ties. When chronology is unavailable, put rows in the intended temporal order.

The default evaluation cohort is all patients in the gold file. Gold episodes for a patient with no prediction rows count as misses. Prediction-only patients are excluded because their manual review status is unknown, and their number is reported in `run.json`. Include a zero-episode gold row to evaluate extra predictions for a reviewed patient with no gold episodes. `--patient-scope common` restores the source's intersection-only patient scope when reproducing that behavior.

For different column names, copy and edit [column_mapping.example.json](column_mapping.example.json). Each mapping is **canonical column name → your CSV column name**, separately for gold and predicted. Remove mappings for columns you do not have. The mapping renames columns; it does not invert status coding or derive clinical dates.

```bash
python evaluate.py \
  --gold /path/to/manual.csv \
  --predicted /path/to/extracted.csv \
  --mapping /path/to/mapping.json \
  --output results/stanford-validation
```

Dates must be complete ISO calendar dates. Invalid dates, nonbinary nonblank event flags, and unknown nonblank U-codes raise errors instead of silently becoming nulls. Translate Stanford event/censor coding explicitly. `first_pd_date` means progression date: if the local PFS event definition includes death without progression, agree on how that endpoint is represented before interpreting this date metric.

## Alignment

For each patient, Needleman–Wunsch aligns the two episode sequences using:

```text
similarity = 0.5 × regimen drug-set Jaccard
           + 0.3 × max(0, 1 − |start difference| / 60 days)
           + 0.2 × max(0, 1 − |stop difference| / 90 days)

pair score = similarity              if similarity >= 0.5
           = 0.5 × similarity − 0.1  if similarity >= 0.2
           = −0.4                    otherwise

gap penalty = −0.5
```

Missing alignment components contribute zero; weights are not renormalized. Traceback prefers a diagonal pair, then a gold-only gap, then a prediction-only gap on ties. Progression dates, censor dates, event flags, and uncertainty codes do not influence the alignment. Field tolerances are applied after alignment.

An aligned pair means a sequence correspondence, not a correct clinical match. A low-similarity pair can score −0.4, which is less costly than two gaps at −1.0. Equal episode counts can therefore produce high line recovery despite poor field agreement. Review `alignment.csv` alongside field metrics. Start and stop dates contribute to the alignment and are subsequently scored, as in the original method.

Regimens are expanded and normalized using the original dictionaries, with leucovorin removed; regimen field agreement uses Jaccard ≥ 0.5. The dictionaries originated in GI cohorts. For lung cases, use the same canonical generic drug names on both sides. Brand names, multiword names, and regimen abbreviations need review in `alignment.py`; this package does not claim a validated lung-specific synonym dictionary. Changing normalization changes the alignment and should be recorded for the validation run.

## Date tolerances

The default `--date-mode uncertainty` uses each **gold date's** code:

| Gold code | Allowed absolute difference, days |
| --- | ---: |
| U0 | 7 |
| U1 | 14 |
| U2 | 31 |
| U3 | 93 |
| U4 | 183 |
| Blank or absent | 30 |

The boundary is inclusive. The public pipeline's fixed comparator is available with `--date-mode fixed`: start/progression ±14 days, stop/censor ±30 days. Always report which mode was used. Date proximity windows used for alignment are separate from these field agreement tolerances.

## Outputs and denominators

- `metrics.csv`: values on a 0–1 scale, numerators, denominators, and patient-clustered 95% percentile bootstrap intervals. Event confusion entries are counts in `numerator`, with value and denominator blank.
- `alignment.csv`: one row per correspondence, missed gold episode, or extra prediction, including source CSV row numbers, compared values, applied tolerances, and match flags. Keep this file within the local clinical environment when using real data.
- `run.json`: cohort counts, fields scored/unavailable, configuration, column mapping, and input SHA-256 hashes.

Line sensitivity = matched / (matched + missed), PPV = matched / (matched + extra), F1 = 2 × matched / (2 × matched + missed + extra). Empty denominators are unavailable. A one-sided nonempty timeline yields zero F1 here; the source paper helper returned unavailable when precision or sensitivity was undefined.

Field metrics use aligned pairs only. Gaps affect line recovery, not field agreement. The main `*_agreement` uses the source's null semantics: absent/absent dates, flags, and labels agree; one absent value disagrees; regimen comparison is unavailable if either drug set is absent. Entirely unavailable columns produce no field metric.

Each date also has `*_acc_gold_date_present` and `*_acc_both_date_present`. These are descriptive subsets of matched pairs, not episode-recovery-inclusive measures. Five conditional PFS date metrics follow the paper's subset definitions: progression when both flags are event, censor date when both flags are censored, progression/censor when both dates are present, and progression when the gold flag is event. The last still uses aligned pairs only. The both-flag subsets retain null/null agreement, exactly as the source helper does; read them with the both-date-present metrics.

Event confusion counts use known binary flags on aligned pairs. `pfs_flag_agreement_both_known` excludes unknown flags, while `pfs_flag_agreement` retains the source's null semantics. False negative means gold event/predicted censored; false positive means gold censored/predicted event. These counts are distinct from missed or extra episodes.

Bootstrap defaults: 2,000 samples, seed 20260904, sampling patients with replacement and retaining all episodes per sampled patient. The wrapper samples the configured cohort, including patients with zero alignment rows. Replicates with undefined denominators are excluded separately for each metric; `n_valid_bootstrap` reports the surviving count. Use `--bootstrap 0` for point estimates only. No confidence intervals are computed for raw confusion counts. Run each manual reviewer separately; there is no implicit pooling or reviewer averaging.

Only code and synthetic examples are included. Clinical inputs and results should remain local, outside version control. Contact: Cynthia Yeung, ckyeung@mdanderson.org.
