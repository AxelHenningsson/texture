import numpy as np
import matplotlib.pyplot as plt
from xfab.symmetry import Umis

def kernel_average_misorientation(orientation_map, footprint, crystal_system, misorientation_threshold=None):
    """Apply a KAM filter to the input orientation map.
    """
    m = footprint.shape[0] // 2
    n = footprint.shape[1] // 2
    N = np.sum(footprint)
    kam_map = np.zeros((orientation_map.shape[0], orientation_map.shape[1]))
    for i in range(orientation_map.shape[0]):
        for j in range(orientation_map.shape[0]):
            U = orientation_map[i,j]
            mis = 0
            for k in range(footprint.shape[0]):
                for l in range(footprint.shape[1]):
                    if footprint[k,l]>0:
                        row = i - m + k
                        col = j - n + l
                        u = orientation_map[row, col]
                        mis += Umis( U, u, crystal_system )[:,1].min()
            kam_map[i,j] += (mis/N)
    return kam_map

if __name__ == "__main__":
    pass