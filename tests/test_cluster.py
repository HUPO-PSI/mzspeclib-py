import os
import unittest

from mzspeclib.backends import TextSpectralLibrary
from mzspeclib.cluster import SpectrumCluster

from .common import datafile


class TestSpectrumCluster(unittest.TestCase):

    def get_library(self):
        test_file = datafile("clusters_example.mzspeclib.txt")
        return TextSpectralLibrary(test_file)

    def test_text_cluster_parsing(self):
        lib = self.get_library()
        cluster: SpectrumCluster = lib.get_cluster(1)

        assert cluster.key == 1
        assert cluster.size == 6