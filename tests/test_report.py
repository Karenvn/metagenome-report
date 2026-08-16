from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

from metagenome_report.report import write_table


class ReportTableTests(unittest.TestCase):
    def test_table_formats_size_and_omits_taxid(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "bin_data.csv"
            source.write_text(
                "bin_id,size,taxon_id,Completeness,Contamination\n"
                "bin-01,5044489,1891229,92.3,1.2\n",
                encoding="utf-8",
            )

            md_path, csv_path = write_table(source, root / "output")

            markdown = md_path.read_text(encoding="utf-8")
            summary = pd.read_csv(csv_path, dtype=str, keep_default_na=False)

        self.assertIn("5,044,489", markdown)
        self.assertEqual(summary.loc[0, "size"], "5,044,489")
        self.assertNotIn("taxon_id", summary.columns)
        self.assertNotIn("taxon_id", markdown)


if __name__ == "__main__":
    unittest.main()
