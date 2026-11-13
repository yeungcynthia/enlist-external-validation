import csv
from datetime import datetime, timedelta
from pathlib import Path
import tempfile
import unittest

import alignment as nw
from evaluate import load_csv, align, metric_counts, summarize, FIELDS

ROOT = Path(__file__).resolve().parent


class ValidationTests(unittest.TestCase):
    def fixture(self):
        gold, gc = load_csv(ROOT / 'example_gold.csv')
        pred, pc = load_csv(ROOT / 'example_predicted.csv')
        fields = list(FIELDS.values())
        return gold, pred, fields

    def test_example_counts_and_uncertainty(self):
        g, p, fs = self.fixture()
        rows, patients, audit = align(g, p, fs)
        stats = metric_counts(rows, fs)
        self.assertEqual(len(patients), 3)
        self.assertEqual(stats['line_sensitivity'], (3, 4, 0.75))
        self.assertEqual(stats['line_ppv'], (3, 3, 1))
        self.assertAlmostEqual(stats['line_f1'][2], 6 / 7)
        self.assertEqual(stats['start_agreement'], (2, 3, 2 / 3))
        self.assertEqual(stats['first_pd_acc_both_event_flag'], (1, 1, 1))
        self.assertEqual(stats['first_pd_acc_gold_event'], (1, 2, 0.5))
        self.assertEqual(stats['event_false_negative'][0], 1)
        self.assertEqual(audit['n_gold_patients_without_prediction_rows'], 1)

    def test_fixed_mode_is_distinct(self):
        g, p, fs = self.fixture()
        rows, _, _ = align(g, p, fs, mode='fixed')
        self.assertEqual(metric_counts(rows, fs)['start_agreement'][2], 1)

    def test_uncertainty_boundaries_and_fallback(self):
        dt = datetime(2020, 1, 1)
        for code, tolerance in [('U0', 7), ('U1', 14), ('U2', 31), ('U3', 93), ('U4', 183), (None, 30)]:
            self.assertEqual(nw.tol_for_u(code), tolerance)
            self.assertTrue(nw.dates_match(dt, dt + timedelta(days=tolerance), tolerance))
            self.assertFalse(nw.dates_match(dt, dt + timedelta(days=tolerance + 1), tolerance))

    def test_null_semantics(self):
        self.assertTrue(nw.dates_match(None, None, 7))
        self.assertFalse(nw.dates_match(None, datetime(2020, 1, 1), 7))

    def test_empty_and_one_sided_sequences(self):
        g, _, fs = self.fixture()
        self.assertEqual(nw.needleman_wunsch([], []), [])
        self.assertEqual(nw.needleman_wunsch(g[:1], []), [(0, None)])
        rows, _, _ = align(g, [], fs)
        stats = metric_counts(rows, fs)
        self.assertEqual(stats['line_sensitivity'][2], 0)
        self.assertIsNone(stats['line_ppv'][2])
        self.assertEqual(stats['line_f1'][2], 0)
        self.assertIsNone(metric_counts([], fs)['line_f1'][2])

    def test_extra_middle_episode_aligns_chronologically(self):
        def ep(drug, day):
            return {'drugs': frozenset([drug]), 'start': datetime(2020, 1, 1) + timedelta(days=day), 'stop': None}
        gold = [ep('a', 0), ep('b', 200)]
        pred = [ep('a', 0), ep('x', 100), ep('b', 200)]
        self.assertEqual(nw.needleman_wunsch(gold, pred), [(0, 0), (None, 1), (1, 2)])

    def test_low_similarity_still_pairs_as_original(self):
        a = {'drugs': frozenset(['a']), 'start': None, 'stop': None}
        b = {'drugs': frozenset(['b']), 'start': None, 'stop': None}
        self.assertEqual(nw.needleman_wunsch([a], [b]), [(0, 0)])

    def test_intersection_scope_is_explicit(self):
        g, p, fs = self.fixture()
        rows, patients, _ = align(g, p, fs, scope='common')
        self.assertEqual(len(patients), 2)
        self.assertEqual(metric_counts(rows, fs)['line_sensitivity'][2], 1)

    def test_bootstrap_reproducible_and_zero_iterations(self):
        g, p, fs = self.fixture()
        rows, patients, _ = align(g, p, fs)
        a = summarize(rows, patients, fs, bootstrap=50)
        self.assertEqual(a, summarize(rows, patients, fs, bootstrap=50))
        self.assertTrue(all(r['ci_low'] is None for r in summarize(rows, patients, fs, bootstrap=0)))
        self.assertTrue(all(r['ci_low'] is None for r in summarize([], [], fs, bootstrap=50)))

    def csv_file(self, text):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        path = Path(temporary.name) / 'input.csv'
        path.write_text(text)
        return path

    def test_original_id_normalization_and_episode_filter(self):
        path = self.csv_file('sid,begin\n001,2020-01-01\n')
        rows, columns = load_csv(path, {'patient_id': 'sid', 'start_date': 'begin'})
        self.assertEqual(rows[0]['patient_id'], '1')
        self.assertTrue(rows[0]['placeholder'])
        self.assertNotIn('pfs_censor_flag', columns)
        detail, _, _ = align(rows, rows, list(FIELDS.values()))
        self.assertEqual(detail, [])

    def test_reviewed_zero_episode_patient_counts_extra_predictions(self):
        g = self.csv_file('patient_id,enlist_lot_label,start_date\nSYN001,,\n')
        p = self.csv_file('patient_id,enlist_lot_label,start_date\nSYN001,1,2020-01-01\nSYN999,1,2021-01-01\n')
        rows, _, audit = align(load_csv(g)[0], load_csv(p)[0], list(FIELDS.values()))
        self.assertEqual([r['alignment_type'] for r in rows], ['agent_extra'])
        self.assertEqual(audit['n_prediction_only_patients_excluded'], 1)
        self.assertEqual(metric_counts(rows, ['start'])['line_f1'][2], 0)

    def test_original_permissive_dates_and_flags(self):
        path = self.csv_file('patient_id,enlist_lot_label,start_date,stop_date,pfs_censor_flag\n001,1,01/02/2020,ongoing,yes\n001,2,not-a-date,,1.0\n001,3,2020-02-30,,unavailable\n001,4,,,2\n')
        rows, _ = load_csv(path)
        self.assertEqual(rows[0]['start'], datetime(2020, 1, 2))
        self.assertIsNone(rows[0]['stop'])
        self.assertEqual(rows[0]['pfs_flag'], 1)
        self.assertIsNone(rows[1]['start'])
        self.assertEqual(rows[1]['pfs_flag'], 1)
        self.assertIsNone(rows[2]['start'])
        self.assertIsNone(rows[2]['pfs_flag'])
        self.assertEqual(rows[3]['pfs_flag'], 2)

    def test_absent_columns_are_compared_as_missing(self):
        path = self.csv_file('patient_id,enlist_lot_label\n001,1\n')
        rows, _ = load_csv(path)
        detail, _, _ = align(rows, rows, list(FIELDS.values()))
        stats = metric_counts(detail, list(FIELDS.values()))
        self.assertEqual(stats['start_agreement'], (1, 1, 1))
        self.assertEqual(stats['pfs_flag_agreement'], (1, 1, 1))
        self.assertEqual(stats['regimen_agreement'], (0, 0, None))

    def test_blank_ids_are_dropped(self):
        path = self.csv_file('patient_id,enlist_lot_label\n,1\n001,1\n')
        rows, _ = load_csv(path)
        self.assertEqual([r['patient_id'] for r in rows], ['1'])

    def test_malformed_csv_structure_is_rejected(self):
        path = self.csv_file('patient_id,start_date,start_date\nSYN001,,\n')
        with self.assertRaises(ValueError):
            load_csv(path)

    def test_gold_date_present_excludes_null_agreement(self):
        path = self.csv_file('patient_id,enlist_lot_label,start_date,first_pd_date\nSYN001,1,2020-01-01,\n')
        rows = load_csv(path)[0]
        detail, _, _ = align(rows, rows, ['start', 'first_pd'])
        stats = metric_counts(detail, ['start', 'first_pd'])
        self.assertEqual(stats['first_pd_agreement'], (1, 1, 1))
        self.assertEqual(stats['first_pd_acc_gold_date_present'], (0, 0, None))


if __name__ == '__main__':
    unittest.main()
