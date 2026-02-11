import os
import tempfile
import unittest
from pathlib import Path

from p53cad.core.runtime import bootstrap_runtime
from p53cad.engine.drug_generator import DrugGeneratorEngine

bootstrap_runtime()


class DrugRuntimeTests(unittest.TestCase):
    def test_mode_capabilities_template_ready(self):
        engine = DrugGeneratorEngine()
        caps = engine.get_mode_capabilities("Y220C_cavity", "template")
        self.assertTrue(caps["ready"])

    def test_mode_capabilities_docking_reports_missing_runtime(self):
        with tempfile.TemporaryDirectory() as tmp:
            old_cwd = Path.cwd()
            try:
                os.chdir(tmp)
                engine = DrugGeneratorEngine()
                caps = engine.get_mode_capabilities(
                    pocket_key="Y220C_cavity",
                    method="docking",
                    allow_wt_receptor_fallback=False,
                )
                self.assertIn("ready", caps)
                self.assertFalse(caps["ready"])
                self.assertTrue(caps["reason"])
            finally:
                os.chdir(old_cwd)

    def test_receptor_manifest_detects_reuse(self):
        with tempfile.TemporaryDirectory() as tmp:
            old_cwd = Path.cwd()
            try:
                os.chdir(tmp)
                receptors = Path("data/raw/receptors")
                receptors.mkdir(parents=True, exist_ok=True)
                shared = receptors / "shared.pdbqt"
                shared.write_text("RECEPTOR\n")
                (receptors / "y220c_cavity.pdbqt").symlink_to(shared.resolve())
                (receptors / "l1_loop.pdbqt").symlink_to(shared.resolve())

                engine = DrugGeneratorEngine()
                manifest = engine.get_receptor_manifest(allow_wt_receptor_fallback=False)
                self.assertIn("duplicate_receptor_paths", manifest)
                self.assertGreaterEqual(manifest["duplicate_count"], 1)
            finally:
                os.chdir(old_cwd)

    def test_meeko_atom_params_resolver_returns_json_token(self):
        engine = DrugGeneratorEngine()
        try:
            token = engine._resolve_meeko_atom_params()
        except RuntimeError:
            self.skipTest("Meeko atom params unavailable in this environment")
            return
        self.assertTrue(token.endswith(".json") or "/" in token or "\\" in token)


if __name__ == "__main__":
    unittest.main()
