"""CSV adapter and metrics for the unchanged EnLiST alignment algorithm."""
import argparse
import csv
from datetime import datetime
import hashlib
import json
from pathlib import Path

import numpy as np

import alignment as nw

DATES = {'start_date': 'start', 'stop_date': 'stop',
         'first_pd_date': 'first_pd', 'censor_date': 'censor_date'}
FIELDS = {'systemic_anticancer_agents': 'regimen', 'enlist_lot_label': 'lot_label',
          **DATES, 'pfs_censor_flag': 'pfs_flag'}
FIXED = {'start': 14, 'stop': 30, 'first_pd': 14, 'censor_date': 30}
U_FIELDS = {'start': 'start_date_U', 'stop': 'stop_date_U',
            'first_pd': 'pd_date_U', 'censor_date': 'censor_date_U'}


def load_csv(path, mapping=None):
    """Adapt CSV columns using Elena's original ID, date, flag, and episode rules."""
    mapping = mapping or {}
    if len(set(mapping.values())) != len(mapping):
        raise ValueError('Column mapping must use distinct source columns')
    with Path(path).open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        headers = reader.fieldnames or []
        if len(headers) != len(set(headers)):
            raise ValueError('Duplicate CSV column headers')
        if any(v not in headers for v in mapping.values()):
            raise ValueError('A mapped source column is absent from the CSV')
        columns = set(headers) | set(mapping)
        if mapping.get('patient_id', 'patient_id') not in headers:
            raise ValueError('patient_id column required (or provide a mapping)')
        rows = []
        for number, raw in enumerate(reader, 2):
            if None in raw:
                raise ValueError(f'CSV row {number}: too many cells')
            row = {k: (raw.get(mapping.get(k, k)) or '').strip() for k in columns}
            patient_id = nw.normalize_mrn(row['patient_id'] or None)
            if not patient_id:
                continue
            episode = {'patient_id': patient_id, 'source_row': number}
            for column, field in FIELDS.items():
                value = row.get(column) or None
                if column in DATES:
                    value = nw.parse_date(value)
                elif column == 'pfs_censor_flag':
                    value = nw.norm_flag(value)
                episode[field] = value
            for field, column in U_FIELDS.items():
                value = row.get(column, '')
                if field == 'first_pd' and not value:
                    value = row.get('first_pd_date_U', '')
                if value and value.upper() not in nw.U_TOLERANCE:
                    raise ValueError(f'CSV row {number}: unknown uncertainty code for {field}')
                episode[field + '_U'] = value.upper() or None
            episode['drugs'] = nw.regimen_to_drug_set(episode['regimen'])
            # Source extract_lines skips rows with neither a regimen nor a LoT label.
            def _blank(v):
                return v is None or str(v).strip() == ''
            episode['placeholder'] = (_blank(episode['regimen']) and
                                      _blank(episode['lot_label']))
            rows.append(episode)
    return rows, columns


def compare(field, gold, pred, mode):
    a, b = gold[field], pred[field]
    if field in FIXED:
        tol = nw.tol_for_u(gold.get(field + '_U')) if mode == 'uncertainty' else FIXED[field]
        return nw.dates_match(a, b, tol), tol
    if field == 'regimen':
        return (nw.jaccard(gold['drugs'], pred['drugs']) >= 0.5
                if gold['drugs'] and pred['drugs'] else None), None
    if field == 'lot_label':
        return nw.norm_lot_label(a) == nw.norm_lot_label(b), None
    return a == b, None


def align(gold, predicted, fields, mode='uncertainty', scope='gold'):
    gold_ids = {r['patient_id'] for r in gold}
    pred_ids = {r['patient_id'] for r in predicted}
    cohort = gold_ids if scope == 'gold' else gold_ids & pred_ids
    detail = []
    for patient in sorted(cohort):
        def episodes(rows):
            return sorted((r for r in rows if r['patient_id'] == patient and not r['placeholder']),
                          key=lambda r: (r['start'] is None, r['start'] or datetime.max))
        gs, ps = episodes(gold), episodes(predicted)
        for gi, pi in nw.needleman_wunsch(gs, ps):
            g, p = gs[gi] if gi is not None else None, ps[pi] if pi is not None else None
            row = {'patient_id': patient,
                   'alignment_type': 'matched' if g and p else 'gs_missed' if g else 'agent_extra',
                   'gold_source_row': g['source_row'] if g else None,
                   'predicted_source_row': p['source_row'] if p else None}
            for f in fields:
                row['gold_' + f] = g[f] if g else None
                row['predicted_' + f] = p[f] if p else None
                row['match_' + f], row['tolerance_' + f] = compare(f, g, p, mode) if g and p else (None, None)
            detail.append(row)
    audit = {'patient_scope': scope, 'n_gold_patients': len(gold_ids),
             'n_predicted_patients': len(pred_ids), 'n_evaluated_patients': len(cohort),
             'n_gold_patients_without_prediction_rows': len(gold_ids - pred_ids),
             'n_prediction_only_patients_excluded': len(pred_ids - gold_ids)}
    return detail, sorted(cohort), audit


def metric_counts(rows, fields):
    """Return numerator, denominator, value for each explicitly named metric."""
    matched = [r for r in rows if r['alignment_type'] == 'matched']
    missed = sum(r['alignment_type'] == 'gs_missed' for r in rows)
    extra = sum(r['alignment_type'] == 'agent_extra' for r in rows)
    m = len(matched)
    stats = {}
    def add(name, numerator, denominator):
        stats[name] = (numerator, denominator, numerator / denominator if denominator else None)
    add('line_sensitivity', m, m + missed)
    add('line_ppv', m, m + extra)
    add('line_f1', 2 * m, 2 * m + missed + extra)
    for field in fields:
        eligible = [r for r in matched if r['match_' + field] is not None]
        add(field + '_agreement', sum(r['match_' + field] for r in eligible), len(eligible))
        if field in FIXED:
            for suffix, subset in [
                ('gold_date_present', [r for r in matched if r['gold_' + field] is not None]),
                ('both_date_present', [r for r in matched if r['gold_' + field] is not None and r['predicted_' + field] is not None])]:
                add(field + '_acc_' + suffix, sum(r['match_' + field] is True for r in subset), len(subset))
    if 'pfs_flag' in fields:
        for name, g, p in [('event_true_positive', 1, 1), ('event_false_negative', 1, 0),
                           ('event_false_positive', 0, 1), ('event_true_negative', 0, 0)]:
            stats[name] = (sum(r['gold_pfs_flag'] == g and r['predicted_pfs_flag'] == p for r in matched), None, None)
        known = [r for r in matched if r['gold_pfs_flag'] in (0, 1) and r['predicted_pfs_flag'] in (0, 1)]
        add('pfs_flag_agreement_both_known', sum(r['gold_pfs_flag'] == r['predicted_pfs_flag'] for r in known), len(known))
        for name, field, flag in [('first_pd_acc_both_event_flag', 'first_pd', 1),
                                  ('censor_date_acc_both_censor_flag', 'censor_date', 0),
                                  ('first_pd_acc_gold_event', 'first_pd', 1)]:
            if field not in fields:
                continue
            subset = [r for r in matched if r['gold_pfs_flag'] == flag and
                      (name.endswith('gold_event') or r['predicted_pfs_flag'] == flag)]
            add(name, sum(r['match_' + field] is True for r in subset), len(subset))
    return stats


def summarize(detail, patients, fields, bootstrap=2000, seed=20260904):
    point = metric_counts(detail, fields)
    values = {name: [] for name in point}
    groups = {p: [r for r in detail if r['patient_id'] == p] for p in patients}
    rng = np.random.default_rng(seed)
    if patients:
        for _ in range(bootstrap):
            sample = rng.choice(patients, size=len(patients), replace=True)
            stats = metric_counts([r for p in sample for r in groups[p]], fields)
            for name, (_, _, value) in stats.items():
                if value is not None:
                    values[name].append(value)
    result = []
    for name, (num, den, value) in point.items():
        vals = values[name]
        lo, hi = np.percentile(vals, [2.5, 97.5]) if vals else (None, None)
        result.append({'metric': name, 'numerator': num, 'denominator': den, 'value': value,
                       'ci_low': lo, 'ci_high': hi, 'n_evaluated_patients': len(patients),
                       'n_valid_bootstrap': len(vals)})
    return result


def write_csv(path, rows, columns=None):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns or list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: v.date().isoformat() if isinstance(v, datetime) else v for k, v in row.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gold', type=Path, required=True)
    parser.add_argument('--predicted', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='New directory, must not exist')
    parser.add_argument('--mapping', type=Path, help='JSON with gold and predicted canonical-to-source maps')
    parser.add_argument('--date-mode', choices=['uncertainty', 'fixed'], default='uncertainty')
    parser.add_argument('--patient-scope', choices=['gold', 'common'], default='gold')
    parser.add_argument('--bootstrap', type=int, default=2000)
    parser.add_argument('--seed', type=int, default=20260904)
    args = parser.parse_args()
    if args.bootstrap < 0:
        parser.error('--bootstrap must be nonnegative')
    mapping = json.loads(args.mapping.read_text()) if args.mapping else {}
    gold, gc = load_csv(args.gold, mapping.get('gold'))
    predicted, pc = load_csv(args.predicted, mapping.get('predicted'))
    # Source stage4 evaluates every predefined field, including absent columns.
    fields = list(FIELDS.values())
    detail, patients, audit = align(gold, predicted, fields, args.date_mode, args.patient_scope)
    summary = summarize(detail, patients, fields, args.bootstrap, args.seed)
    audit.update({'date_mode': args.date_mode, 'fields_scored': fields,
                  'fields_missing_input_columns': [f for c, f in FIELDS.items() if c not in gc or c not in pc],
                  'n_matched': sum(r['alignment_type'] == 'matched' for r in detail),
                  'n_gold_missed': sum(r['alignment_type'] == 'gs_missed' for r in detail),
                  'n_prediction_extra': sum(r['alignment_type'] == 'agent_extra' for r in detail),
                  'bootstrap_iterations': args.bootstrap, 'seed': args.seed,
                  'gold_sha256': hashlib.sha256(args.gold.read_bytes()).hexdigest(),
                  'predicted_sha256': hashlib.sha256(args.predicted.read_bytes()).hexdigest(),
                  'column_mapping': mapping})
    args.output.mkdir(parents=True, exist_ok=False)
    write_csv(args.output / 'metrics.csv', summary)
    write_csv(args.output / 'alignment.csv', detail, list(detail[0]) if detail else ['patient_id', 'alignment_type'])
    (args.output / 'run.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(f"Evaluated {len(patients)} patients; {audit['n_matched']} aligned episode pairs. Results: {args.output}")


if __name__ == '__main__':
    main()
