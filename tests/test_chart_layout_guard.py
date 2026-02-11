import unittest
from pathlib import Path


class ChartLayoutGuardTests(unittest.TestCase):
    def test_dynamic_legend_margin_guard_present(self):
        app_path = Path("p53cad/app/main.py")
        source = app_path.read_text()
        self.assertIn("has_visible_legend", source)
        self.assertIn("min_bottom", source)
        self.assertIn("style_plotly_figure", source)


if __name__ == "__main__":
    unittest.main()
