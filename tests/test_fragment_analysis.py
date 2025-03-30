#!/usr/bin/env python
"""
Tests for the fragment analysis module.
"""

import unittest
import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.topology.core import Topology
from MDAnalysis.core.groups import AtomGroup
from src.fragment_analysis import (
    compute_fragment_metrics,
    compute_fragment_normal,
    compute_orientation_angle,
    compute_contact_fraction
)

class TestFragmentAnalysis(unittest.TestCase):
    """Tests for fragment analysis functions."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a simple test universe with a fragment and Pt atoms
        # Fragment: a triangle in the xy-plane
        frag_positions = np.array([
            [0.0, 0.0, 0.0],  # Origin
            [1.0, 0.0, 0.0],  # 1 unit along x
            [0.0, 1.0, 0.0]   # 1 unit along y
        ])
        
        # Pt atoms: a small cluster
        pt_positions = np.array([
            [5.0, 0.0, 0.0],  # 5 units along x
            [6.0, 0.0, 0.0],  # 6 units along x
            [5.5, 1.0, 0.0]   # Between the two, 1 unit up in y
        ])
        
        # Create the universe
        u = mda.Universe.empty(6, n_residues=2, n_segments=2, trajectory=True)
        
        # Add atoms
        u.add_TopologyAttr('name', ['C1', 'C2', 'C3', 'PT1', 'PT2', 'PT3'])
        u.add_TopologyAttr('type', ['C', 'C', 'C', 'PT', 'PT', 'PT'])
        u.add_TopologyAttr('resname', ['FRAG', 'FRAG', 'FRAG', 'PT', 'PT', 'PT'])
        u.add_TopologyAttr('resid', [1, 1, 1, 2, 2, 2])
        u.add_TopologyAttr('segid', ['FRAG', 'FRAG', 'FRAG', 'PT', 'PT', 'PT'])
        
        # Set positions
        all_positions = np.vstack((frag_positions, pt_positions))
        u.trajectory = MemoryReader(np.array([all_positions]), order='afc')
        
        # Create atom groups
        self.frag = u.atoms[:3]
        self.pt_atoms = u.atoms[3:]
        self.pt_com = np.array([5.5, 1/3, 0.0])  # Center of mass of Pt atoms
        
    def test_compute_fragment_normal(self):
        """Test computing fragment normal vector."""
        normal = compute_fragment_normal(self.frag)
        # Expected normal vector for the triangle in xy-plane is along z-axis
        expected = np.array([0.0, 0.0, 1.0])
        self.assertTrue(np.allclose(normal, expected) or np.allclose(normal, -expected))
    
    def test_compute_orientation_angle(self):
        """Test computing orientation angle."""
        angle = compute_orientation_angle(self.frag, self.pt_com)
        # The angle between z-axis and COM-COM vector should be 90 degrees
        self.assertAlmostEqual(angle, 90.0, places=5)
    
    def test_compute_contact_fraction(self):
        """Test computing contact fraction."""
        # No atoms should be in contact with the threshold of 2.5
        fraction = compute_contact_fraction(self.frag, self.pt_atoms, threshold=2.5)
        self.assertEqual(fraction, 0.0)
        
        # All atoms should be in contact with a large threshold of 10
        fraction = compute_contact_fraction(self.frag, self.pt_atoms, threshold=10.0)
        self.assertEqual(fraction, 1.0)
    
    def test_compute_fragment_metrics(self):
        """Test computing all fragment metrics."""
        metrics = compute_fragment_metrics(self.frag, self.pt_atoms, self.pt_com, close_threshold=2.5)
        
        # Check that all expected keys are present
        expected_keys = ["min_dist", "avg_dist", "med_dist", "std_dist", "com_dist", 
                          "num_close_contacts", "contact_fraction", "orientation_angle"]
        for key in expected_keys:
            self.assertIn(key, metrics)
        
        # Check some specific values
        self.assertGreater(metrics["min_dist"], 0)  # Distance should be positive
        self.assertEqual(metrics["num_close_contacts"], 0)  # No atoms in contact with threshold 2.5
        self.assertEqual(metrics["contact_fraction"], 0.0)  # No atoms in contact
        self.assertAlmostEqual(metrics["orientation_angle"], 90.0, places=5)  # Angle should be 90 degrees

if __name__ == '__main__':
    unittest.main()