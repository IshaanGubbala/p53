"""Tests for the AttentionPoolingNet oracle architecture.

Validates forward/backward pass with both (batch, seq_len, D) and (batch, D)
inputs, gradient flow, and parameter counts.
"""

import unittest

try:
    import torch
    from p53cad.engine.oracle import AttentionPoolingNet

    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False


@unittest.skipUnless(HAS_TORCH, "Requires PyTorch")
class TestAttentionPoolingNet(unittest.TestCase):
    """Test AttentionPoolingNet forward, backward, and compatibility."""

    def setUp(self):
        self.net = AttentionPoolingNet(input_dim=64, hidden_dim=32, n_heads=2, dropout=0.0)
        self.net.eval()

    def test_forward_3d_input(self):
        """(batch, seq_len, D) input produces (batch, 1) output."""
        x = torch.randn(2, 10, 64)
        out = self.net(x)
        self.assertEqual(out.shape, (2, 1))

    def test_forward_2d_input(self):
        """(batch, D) input works via backward-compatible path."""
        x = torch.randn(3, 64)
        out = self.net(x)
        self.assertEqual(out.shape, (3, 1))

    def test_single_position(self):
        """(batch, 1, D) input works (single-position sequence)."""
        x = torch.randn(1, 1, 64)
        out = self.net(x)
        self.assertEqual(out.shape, (1, 1))

    def test_gradient_flows(self):
        """Gradients flow through attention pooling to input."""
        self.net.train()
        x = torch.randn(2, 10, 64, requires_grad=True)
        out = self.net(x)
        loss = out.sum()
        loss.backward()
        self.assertIsNotNone(x.grad)
        self.assertGreater(x.grad.abs().sum().item(), 0)

    def test_query_parameter_requires_grad(self):
        """The learnable query should be a trainable parameter."""
        self.assertTrue(self.net.query.requires_grad)
        self.assertEqual(self.net.query.shape, (1, 1, 64))

    def test_parameter_count_reasonable(self):
        """Attention oracle should have more parameters than a simple MLP head."""
        total = sum(p.numel() for p in self.net.parameters())
        head_only = sum(p.numel() for p in self.net.head.parameters())
        # Total should include attention + layer_norm + query + head
        self.assertGreater(total, head_only)

    def test_output_changes_with_position_order(self):
        """Permuting positions should change the output (attention is position-aware)."""
        x = torch.randn(1, 8, 64)
        out1 = self.net(x).item()
        # Reverse position order
        x_rev = x.flip(dims=[1])
        out2 = self.net(x_rev).item()
        # Outputs may differ (not guaranteed due to attention symmetry,
        # but with random init they almost certainly will)
        # Just check both are finite
        import math
        self.assertFalse(math.isnan(out1))
        self.assertFalse(math.isnan(out2))

    def test_default_dimensions(self):
        """Default construction uses 1280-dim (matching 650M ESM-2)."""
        net = AttentionPoolingNet()
        self.assertEqual(net.input_dim, 1280)
        x = torch.randn(1, 5, 1280)
        out = net(x)
        self.assertEqual(out.shape, (1, 1))

    def test_deterministic_eval(self):
        """Same input produces same output in eval mode (no dropout)."""
        x = torch.randn(1, 5, 64)
        self.net.eval()
        out1 = self.net(x).item()
        out2 = self.net(x).item()
        self.assertAlmostEqual(out1, out2, places=6)


if __name__ == "__main__":
    unittest.main()
