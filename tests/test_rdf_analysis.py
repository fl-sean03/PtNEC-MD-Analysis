#!/usr/bin/env python
"""
Tests for the RDF analysis module.
"""

import unittest
import numpy as np
from src.rdf_analysis import compute_rdf_statistics, combine_rdfs

class TestRDFAnalysis(unittest.TestCase):
    """Tests for RDF analysis functions."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a simple RDF curve with a clear peak and minimum
        self.bins = np.linspace(0, 10, 100)
        # Gaussian peak at r=3 with width 0.5
        self.rdf_values = 1.0 + 2.0 * np.exp(-((self.bins - 3.0) / 0.5)**2)
        # Add a second smaller peak at r=7
        self.rdf_values += 1.0 * np.exp(-((self.bins - 7.0) / 0.3)**2)
    
    def test_compute_rdf_statistics(self):
        """Test computing RDF statistics."""
        stats = compute_rdf_statistics(self.rdf_values, self.bins)
        
        # Check that all expected keys are present
        expected_keys = ["peak_height", "peak_position", "first_min_position", "coordination_number"]
        for key in expected_keys:
            self.assertIn(key, stats)
        
        # Check specific values
        self.assertAlmostEqual(stats["peak_height"], 3.0, places=1)  # Max height should be ~3.0
        self.assertAlmostEqual(stats["peak_position"], 3.0, places=1)  # Peak should be at r~3.0
        
        # Coordination number should be positive
        self.assertGreater(stats["coordination_number"], 0)
    
    def test_combine_rdfs(self):
        """Test combining RDFs from different systems."""
        # Create a second RDF curve with a peak at a different position
        rdf_values2 = 1.0 + 2.0 * np.exp(-((self.bins - 5.0) / 0.5)**2)
        
        # Combine the RDFs
        system_rdf_data = {
            "system1": (self.bins, self.rdf_values),
            "system2": (self.bins, rdf_values2)
        }
        
        combined = combine_rdfs(system_rdf_data)
        
        # Check that the combined data has the expected structure
        self.assertIn("r", combined)
        self.assertIn("system1", combined)
        self.assertIn("system2", combined)
        
        # Check that the arrays have the correct shape
        self.assertEqual(len(combined["r"]), len(self.bins))
        self.assertEqual(len(combined["system1"]), len(self.rdf_values))
        self.assertEqual(len(combined["system2"]), len(rdf_values2))
        
        # Check that the values are correct
        np.testing.assert_array_equal(combined["r"], self.bins)
        np.testing.assert_array_equal(combined["system1"], self.rdf_values)
        np.testing.assert_array_equal(combined["system2"], rdf_values2)

if __name__ == '__main__':
    unittest.main()