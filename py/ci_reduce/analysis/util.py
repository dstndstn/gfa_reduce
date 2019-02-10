import numpy as np
import ci_reduce.common as common

# should probably also have something available for the case of
# upbinning, although not clear when I might actually ever use that
def ci_downbinned_shape(binfac):
    # assume integer rebinning until I come across a case where
    # arbitrary rebinning would be valuable
    # the native CI image dimensions are very favorable for integer rebinning

    # assume same rebinning factor in both dimensions for now, until
    # I come across a case where I would want otherwise

    assert((type(binfac).__name__ == 'int') or binfac.is_integer())

    par = common.ci_misc_params()

    width_native = par['width_pix_native']
    height_native = par['height_pix_native']

    width_downbinned = float(width_native)/float(binfac)
    height_downbinned = float(height_native)/float(binfac)

    assert(width_downbinned.is_integer())
    assert(height_downbinned.is_integer())

    # note Python convention for (height, width)
    return int(height_downbinned), int(width_downbinned)
