# qqconf 1.3.0
* Added function `get_qq_band()` to obtain testing band without making a plot.
* Added argument `prob_pts_method` to control x-coordinates of PP and QQ plots.

# qqconf 1.2.3

* Replaced mm_malloc with fftw_malloc in cpp code for cross-compatibility.

# qqconf 1.2.2

* Changed upper bound in binary search for local level for one_sided and two sided cases.

* Added more accurate calculations for a few grid points for alpha = .01.

* Added fftw3 to SystemRequirements.

* Added bug fix for architectures without mm_malloc.h.

* Slightly changed boundaries at the endpoints of qq and pp plots to enhance visualization.

# qqconf 1.2.1

* Added github repo to URL section in description

# qqconf 1.2.0

* Implemented FFT method for one and two sided bounds to increase speed in `get_level_from_bounds_one_sided()` and `get_level_from_bounds_two_sided()`

# qqconf 1.1.1

* Updated vignette to include table of contents

# qqconf 1.1.0

* Added C code call to `get_level_from_bounds_one_sided()` to improve speed

* Added vignette
