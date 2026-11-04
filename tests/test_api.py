"""Unit tests for the user-facing API independent of image segmentation QC."""

import unittest
from unittest.mock import patch

import numpy as np

from choroid_cvi import extract_cvi


class PublicApiTests(unittest.TestCase):
    def test_one_volume_matches_existing_reductions(self):
        bscans = [np.zeros((8, 8), dtype=np.uint8) for _ in range(5)]
        svgs = [object() for _ in range(5)]
        vessels = [np.ones((8, 8), dtype=np.uint8) for _ in range(5)]
        masks = [np.ones((8, 8), dtype=bool) for _ in range(5)]
        with patch("choroid_cvi.api.segment_oct", return_value=(vessels, masks)) as mock:
            result = extract_cvi(bscans, svgs, spacings=[1, 1, 0.01, 0.01])
        self.assertEqual(mock.call_count, 1)
        self.assertAlmostEqual(result.cvi_mean, 1.0)
        self.assertAlmostEqual(result.cvi_spatial_sd, 0.0)
        self.assertAlmostEqual(result.pachy_area_mm2, 64 * 0.01 * 0.01)

    def test_reject_misaligned_input_counts(self):
        with self.assertRaisesRegex(ValueError, "Five aligned"):
            extract_cvi([np.zeros((8, 8))] * 5, [object()] * 4)


if __name__ == "__main__":
    unittest.main()
