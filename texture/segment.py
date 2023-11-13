import numpy as np
import matplotlib.pyplot as plt
from xfab.symmetry import Umis
from collections import deque
from scipy.ndimage import binary_fill_holes
from scipy.spatial.transform import Rotation

def flood_fill( orientation_map, 
                seed_point, 
                footprint, 
                crystal_system, 
                local_disorientation_tolerance,
                global_disorientation_tolerance, 
                mask=None, 
                background_value=np.nan,
                fill_holes=False,
                max_grains = 99,
                min_grain_size = 0,
                verbose=False
                ):
    assert footprint.shape[0]%2!=0
    assert footprint.shape[1]%2!=0

    M, N = orientation_map.shape[0], orientation_map.shape[1]
    segmentation = np.zeros((M,N))
    skipps = np.ones((M,N), dtype=bool)
    label = 1
    segmentation[seed_point[0], seed_point[1]:seed_point[1] + 3] = label
    done = False
    iteration = 0
    while( not done and iteration < max_grains ):
        rows, cols = np.where( mask*(segmentation==0)*skipps )
        if len(rows)>0:
            n = np.random.randint( 0, len(rows) )
            seed_point = ( rows[n], cols[n] )
            grain_mask = _flood(orientation_map, 
                                seed_point, 
                                footprint, 
                                crystal_system, 
                                local_disorientation_tolerance,
                                global_disorientation_tolerance,
                                mask)
            if fill_holes: grain_mask = binary_fill_holes(grain_mask)
            if np.sum(grain_mask) > min_grain_size:
                segmentation[grain_mask] = label
                label += 1
            else:
                skipps[grain_mask] = False
            iteration += 1
            if verbose:
                print('Iteration ', iteration, ', found grain with : ', np.sum(grain_mask), ' voxels @, seed_point, ', np.sum(grain_mask))
        else:
            done = True

    segmentation[segmentation==0] = background_value
    return segmentation

def _flood( orientation_map, 
            seed_point, 
            footprint, 
            crystal_system,
            local_disorientation_tolerance,
            global_disorientation_tolerance,
            mask ):
    m = footprint.shape[0] // 2
    n = footprint.shape[1] // 2
    flood_mask = np.zeros((orientation_map.shape[0], orientation_map.shape[1]), dtype=bool)

    i, j = seed_point
    flood_mask[i, j] = True
    if mask is not None and mask[i, j]==0: raise ValueError('Seed point not in mask')
    
    unchartered_indices = deque([ seed_point ])
    while(  len(unchartered_indices) > 0 ):
        i, j = unchartered_indices.pop()
        U = orientation_map[i, j]
        
        if np.sum(flood_mask) < 20:
            global_U = Rotation.from_matrix( orientation_map[flood_mask,:,:] ).mean().as_matrix()

        for k in range(footprint.shape[0]):
            for l in range(footprint.shape[1]):
                if footprint[k, l]>0:
                    row = i - m + k
                    col = j - n + l
                    if not flood_mask[row, col] and mask is not None and mask[row, col]:
                        u = orientation_map[row, col]
                        misorientation = Umis( U, u, crystal_system )[:,1].min()
                        
                        if misorientation < local_disorientation_tolerance:
                            glob_misorientation = Umis( global_U, u, crystal_system )[:,1].min()
                            if glob_misorientation < global_disorientation_tolerance:
                                flood_mask[row, col] = True
                                unchartered_indices.appendleft(  ( row, col) )
    return flood_mask

if __name__ == "__main__":
    pass