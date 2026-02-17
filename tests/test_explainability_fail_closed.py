import unittest

from p53cad.core.runtime import bootstrap_runtime
from p53cad.engine.explainability import ExplainabilityEngine

bootstrap_runtime()


class ExplainabilityFailClosedTests(unittest.TestCase):
    def test_probe_reports_unavailable_without_model_stack(self):
        engine = ExplainabilityEngine(model=None, tokenizer=None, oracle=None, embedder=None)
        probe = engine.probe_attention_backend("MEEPQSDPSV")
        self.assertEqual(probe.get("ready"), "false")
        self.assertTrue(probe.get("reason"))


if __name__ == "__main__":
    unittest.main()
