import os
import unittest

from p53cad.core.runtime import bootstrap_runtime, get_runtime_capabilities

bootstrap_runtime()


class RuntimeBootstrapTests(unittest.TestCase):
    def test_bootstrap_sets_required_env_flags(self):
        old_kmp = os.environ.get("KMP_DUPLICATE_LIB_OK")
        old_mps = os.environ.get("PYTORCH_ENABLE_MPS_FALLBACK")
        old_tok = os.environ.get("TOKENIZERS_PARALLELISM")
        try:
            os.environ.pop("KMP_DUPLICATE_LIB_OK", None)
            os.environ.pop("PYTORCH_ENABLE_MPS_FALLBACK", None)
            os.environ.pop("TOKENIZERS_PARALLELISM", None)
            bootstrap_runtime(seed=42)
            self.assertEqual(os.environ.get("KMP_DUPLICATE_LIB_OK"), "TRUE")
            self.assertEqual(os.environ.get("PYTORCH_ENABLE_MPS_FALLBACK"), "1")
            self.assertEqual(os.environ.get("TOKENIZERS_PARALLELISM"), "false")
            self.assertEqual(os.environ.get("PYTHONHASHSEED"), "42")
        finally:
            if old_kmp is not None:
                os.environ["KMP_DUPLICATE_LIB_OK"] = old_kmp
            if old_mps is not None:
                os.environ["PYTORCH_ENABLE_MPS_FALLBACK"] = old_mps
            if old_tok is not None:
                os.environ["TOKENIZERS_PARALLELISM"] = old_tok

    def test_runtime_capabilities_probe_returns_core_keys(self):
        caps = get_runtime_capabilities()
        for key in [
            "python",
            "platform",
            "torch_installed",
            "vina_cli_available",
            "openmm_installed",
        ]:
            self.assertIn(key, caps)


if __name__ == "__main__":
    unittest.main()
