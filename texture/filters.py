import numpy as np
import matplotlib.pyplot as plt
from xfab.symmetry import Umis

def kernel_average_misorientation(orientation_map,
                                  footprint,
                                  crystal_system,
                                  misorientation_threshold=np.inf,
                                  mask=None,
                                  fill_value=np.nan):
    """Apply a kernel average misorientation filter to the input orientation map.

    Args:
        orientation_map (:obj:`numpy array`): the pixelated orientation matrix field, shape=(M, N, 3, 3).
        footprint (:obj:`numpy array`): boolean array defining the kenrel neighbourhood, shape=(m, n).
        crystal_system (:obj:int): crystal_system number must be one of 1: Triclinic, 2: Monoclinic,
            3: Orthorhombic, 4: Tetragonal, 5: Trigonal, 6: Hexagonal, 7: Cubic
        misorientation_threshold (:obj:`numpy array`): Reject misorientations above this threshold to supress noise in the kam
            map. Defaults to np.inf.
        mask (:obj:`numpy array`): Boolean array, where to skipp pixels, shape=(M, N). Defaults to None.
        fill_value (:obj:`numpy array`): To put where mask is false. Defaults to np.nan.

    Returns:
        :obj:`numpy array`: the scalar kernel average misorientation map, shape=(M, N).
    """
    m = footprint.shape[0] // 2
    n = footprint.shape[1] // 2
    kam_map = np.zeros((orientation_map.shape[0], orientation_map.shape[1]))
    for i in range(m, orientation_map.shape[0]-m ):
        for j in range(n, orientation_map.shape[0]-n):
            if mask is not None and mask[i,j]:
                U = orientation_map[i, j]
                local_average_misorientation = 0
                number_of_pixels = 0
                for k in range(footprint.shape[0]):
                    for l in range(footprint.shape[1]):
                        if footprint[k, l]>0:
                            row = i - m + k
                            col = j - n + l
                            if mask is not None and mask[row, col]:
                                u = orientation_map[row, col]
                                misorientation = Umis( U, u, crystal_system )[:,1].min()
                                if misorientation < misorientation_threshold:
                                    local_average_misorientation += misorientation
                                    number_of_pixels += 1
                if number_of_pixels>0:
                    kam_map[i, j] += (local_average_misorientation / number_of_pixels)
            else:
                kam_map[i, j] = fill_value
    return kam_map

if __name__ == "__main__":
    pass