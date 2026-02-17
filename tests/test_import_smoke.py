import unittest

from p53cad.core.runtime import bootstrap_runtime

bootstrap_runtime()


class ImportSmokeTests(unittest.TestCase):
    def test_core_imports(self):
        import p53cad.core.runtime  # noqa: F401
        import p53cad.data.dms  # noqa: F401
        import p53cad.engine.oracle  # noqa: F401
        import p53cad.engine.drug_generator  # noqa: F401
        import p53cad.engine.explainability  # noqa: F401


if __name__ == "__main__":
    unittest.main()
