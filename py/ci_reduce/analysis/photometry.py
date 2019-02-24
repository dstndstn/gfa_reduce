import ci_reduce.analysis.segment as segment
import numpy as np
from photutils import source_properties
from photutils import deblend_sources

def get_segmentation_map(image, bitmask, extname, get_kernel=True):
    # here image is just a 2D numpy array, not a CI_image object

    filler_value = np.median(image)

    # revisit this for saturation case, where I probably want to handle
    # this differently -- also should do something like djs_maskinterp
    # instead anyway

    image[bitmask != 0] = filler_value

    segmap, kernel = segment.segmentation_map(image, extname, get_kernel=True)

    if not get_kernel:
        return segmap
    else:
        return segmap, kernel

def do_deblend_sources(segmap, image, kernel):

    npixels = 5
    segmap_deblend = deblend_sources(image, segmap, npixels,
                                     filter_kernel=kernel, nlevels=32,
                                     contrast=0.001)
    return segmap_deblend

def get_source_list(image, bitmask, extname):
    segmap, kernel = get_segmentation_map(image, bitmask, extname, 
                                          get_kernel=True)

    segmap_deblend = do_deblend_sources(segmap, image, kernel)

    cat = source_properties(image-np.median(image), segmap_deblend)

    cat = cat.to_table()

    return cat
