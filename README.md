# texture
A collection of utility functions for analyzing crystallographic texture maps. Curently just an implementation of the [Kernel Average Misorientation filter.](https://mtex-toolbox.github.io/EBSDKAM.html). To achive something like the below, please see the Kernel Average Misorientation Example section

![image](https://github.com/AxelHenningsson/texture/assets/31615210/fa3a7dd2-7b94-437f-a1f5-01fda3f1977c)

# Kernel Average Misorientation Example
Given that `umap` is your `shape=(M,N,3,3)` orientation matrix numpy array map you may do something like this
````python
import texture.filters
kam_map = texture.filters.kernel_average_misorientation(umap,
                                                        footprint = np.ones((5, 5)),
                                                        crystal_system = 7,
                                                        mask = sample_mask)
````

The `footprint` controls the size of the kernel while the `sample_mask` tells the filter what pixels to skip in `umap`. In this case the crystal system is `crystal_system=7` cubic. Please see the docstrings for more details.

# Install
You may install using pip from the source repo.
````
git clone https://github.com/AxelHenningsson/texture.git && cd texture
pip install -e .
````
