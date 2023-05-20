# Test cases from https://github.com/maximilianh/crisporWebsite/blob/8063ba87b08b5f8af8f880c661c0c06e11a526df/doenchScore.py

import math
from labplotlib.crispy.doenchScore import calcDoenchScore

def test1():
    assert math.isclose(calcDoenchScore("TATAGCTGCGATCTGAGGTAGGGAGGGACC"), 0.713089368437, abs_tol=1e-7)

def test2():
    assert math.isclose(calcDoenchScore("TCCGCACCTGTCACGGTCGGGGCTTGGCGC"), 0.0189838463593, abs_tol=1e-7)
