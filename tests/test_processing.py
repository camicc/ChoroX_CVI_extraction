"""Small deterministic regression tests; no real patient data."""

import unittest

import numpy as np

from choroid_cvi.cli import load_parameters
from choroid_cvi.processing import calculate_cvi, calculate_pachy_area


class CalculationTests(unittest.TestCase):
    def test_cvi_median_of_per_scan_column_statistics(self):
        mask = np.ones((2, 3), dtype=bool)
        scan_a = np.array([[0, 1, 1], [0, 0, 1]], dtype=np.uint8)
        scan_b = np.ones((2, 3), dtype=np.uint8)
        mean, sd = calculate_cvi([scan_a, scan_b], [mask, mask])
        self.assertAlmostEqual(mean, (0.5 + 1.0) / 2)
        self.assertAlmostEqual(sd, np.std([0, 0.5, 1]) / 2)

    def test_zero_mask_columns_are_currently_included_as_zero(self):
        mask = np.array([[1, 0], [1, 0]], dtype=bool)
        vessels = np.array([[1, 0], [1, 0]], dtype=np.uint8)
        mean, sd = calculate_cvi([vessels], [mask])
        self.assertEqual(mean, 0.5)
        self.assertEqual(sd, 0.5)

    def test_largest_component_area(self):
        binary = np.zeros((8, 8), dtype=np.uint8)
        binary[1:3, 1:4] = 1  # 6 pixels
        binary[5, 5] = 1
        area = calculate_pachy_area([binary], spacings=[0, 0, 0.01, 0.02])
        self.assertAlmostEqual(area, 6 * 0.01 * 0.02)

    def test_reject_even_threshold_window(self):
        import json
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "params.json"
            path.write_text(json.dumps({"block_size": 40}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "odd integer"):
                load_parameters(path)


if __name__ == "__main__":
    unittest.main()
