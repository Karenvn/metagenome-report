from __future__ import annotations

import unittest

from metagenome_report.tree_figure import _rectangular_y_limits


class RectangularTreeLayoutTests(unittest.TestCase):
    def test_outer_rows_use_half_row_padding(self) -> None:
        self.assertEqual(
            _rectangular_y_limits([0.0, 1.0, 2.0, 3.0, 4.0]),
            (-0.5, 4.5),
        )

    def test_empty_layout_retains_symmetric_limits(self) -> None:
        self.assertEqual(_rectangular_y_limits([]), (-0.5, 0.5))

    def test_single_row_is_padded_around_its_position(self) -> None:
        self.assertEqual(_rectangular_y_limits([3.0]), (2.5, 3.5))


if __name__ == "__main__":
    unittest.main()
