# Translation of Mark Tapley's Mathematica notebook, OrbitTools.tm, section L, "delta V and rocket equation"

# AMZ note:  Wonder if we sxhould protect against non-double imputs, though should be a problem as there are hard-coded decimals in each equation.
# written/completed 9/22/2017

import numpy as np
from orbittools import *

def main():
    print("test of massRatioRocket(3.05112,295.0), answer should be 2.86999")
    print(massRatioRocket(3.05112,295.0))

    print("test of dVRocket(215.0,1.25), should be 0.470643")
    print(dVRocket(215.0,1.25))

if __name__ == "__main__":
    main()

