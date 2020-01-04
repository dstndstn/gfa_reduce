import ci_reduce.common as common
import ci_reduce.imred.load_calibs as load_calibs
import ci_reduce.dark_current as dark_current
import astropy.io.fits as fits
import numpy as np
import ci_reduce.analysis.util as util

class CI_exposure:
    """Object encapsulating the contents of a single CI exposure"""

    def __init__(self, image_list, exp_header=None, bintables=None):
        # images is a dictionary of CI_image objects

        par = common.ci_misc_params()
        self.images = dict(zip(common.valid_image_extname_list(), 
                               par['n_cameras']*[None]))

        self.assign_image_list(image_list)

        # exposure-level header
        self.exp_header = exp_header
        self.pixels_calibrated = None
        self.bintables = bintables
        
    def assign_one_image(self, image):
        extname = (image.header)['EXTNAME']
        self.images[extname] = image

    def assign_image_list(self, image_list):
        for image in image_list:
            self.assign_one_image(image)

    def subtract_bias(self):
        print('Attempting to subtract bias...')
        for extname in self.images.keys():
            if self.images[extname] is not None:
                self.images[extname].image = self.images[extname].image - \
                    load_calibs.read_bias_image(extname)

                amps = common.valid_amps_list()
                for amp in amps:
                    bdy = common.amp_bdy_coords(amp)
                    self.images[extname].image[bdy['y_l']:bdy['y_u'],
                                               bdy['x_l']:bdy['x_u']] -= \
                        self.images[extname].overscan.overscan_medians[amp]
                
                self.images[extname].bias_subtracted = True
                self.images[extname].calc_variance_e_squared()

    def subtract_dark_current(self):
        print('Attempting to subtract dark current...')
        for extname in self.images.keys():
            if self.images[extname] is not None:
                
                acttime = self.images[extname].try_retrieve_meta_keyword('EXPTIME')
                if acttime is None:
                    print('trying REQTIME instead of EXPTIME')
                    acttime = self.images[extname].try_retrieve_meta_keyword('REQTIME')
                if acttime is None:
                    print('could not find an exposure time !!!!')
                    assert(False) # die

                t_c = self.images[extname].try_retrieve_meta_keyword('GCCDTEMP')
                if t_c is None:
                    print('trying GCOLDTEC instead of GCCDTEMP')
                    t_c = self.images[extname].try_retreive_meta_keyword('GCOLDTEC')
                if t_c is None:
                    print('could not find a CCD temperature !!!!')
                    t_c = 11.0 # HACK !!!!
                
                dark_image = dark_current.total_dark_image_adu(extname,
                                                               acttime, t_c)
                self.images[extname].image = self.images[extname].image - \
                    dark_image
                self.images[extname].dark_subtracted = True
                self.images[extname].create_dq_mask(dark_image)

    def apply_flatfield(self):
        print('Attempting to apply flat field...')
        for extname in self.images.keys():
            if self.images[extname] is not None:
                flatfield = load_calibs.read_flat_image(extname)
                self.images[extname].image = self.images[extname].image / \
                    flatfield
                self.images[extname].flatfielded = True
                self.images[extname].calc_variance_adu(flatfield)

    def calibrate_pixels(self):
        self.subtract_bias()
        self.subtract_dark_current()
        self.apply_flatfield()
        self.pixels_calibrated = True

    def num_images_populated(self):
        return sum( im != None for im in self.images.values() )

    def populated_extnames(self):
        return [k for k,v in self.images.items() if v is not None]

    def to_hdulist(self, flavor=''):
        extname_list = common.valid_image_extname_list()
        hdulist = []
        for extname in extname_list:
            if self.images[extname] is not None:
                hdu = self.images[extname].to_hdu(primary=(len(hdulist) == 0), 
                                                  flavor=flavor)
                hdulist.append(hdu)

        return fits.HDUList(hdulist)

    def estimate_all_sky_mags(self, careful_sky=False):
        for im in self.images.values():
            if im is not None:
                im.estimate_sky_mag(careful_sky=careful_sky)

    def estimate_all_sky_sigmas(self, careful_sky=False):
        for im in self.images.values():
            if im is not None:
                im.set_empirical_bg_sigma(careful_sky=careful_sky)
                print('empirical ' + im.header['EXTNAME'] + 
                      ' background sigma = ' + 
                      '{:.2f}'.format(im.empirical_bg_sigma) + ' ADU')

    def all_source_catalogs(self):
        tables = dict(zip(self.images.keys(), len(self.images.keys())*[None]))

        for extname, im in self.images.items():
            if im is not None:
                tables[extname] = im.catalog_sources()

        # do I also want to store the tables as an attribute belonging to
        # this exposure object?
        return tables

    def update_wcs(self, astr):

        for a in astr:
            extname = a['extname']
            self.images[extname].update_wcs(a)

    def recompute_catalog_radec(self, cat):

        extnames = np.unique(cat['camera'])

        for extname in extnames:
            cat[cat['camera'] == extname] = self.images[extname].catalog_add_radec(cat[cat['camera'] == extname])

    def set_bintable_rows(self):
        for image in self.images.values():
            extname = image.header['EXTNAME'].strip()
            if image.cube_index is None:
                image.bintable_row = None
            elif image.cube_index == -1:
                image.bintable_row = util.average_bintable_metadata(self.bintables[extname])
            else:
                image.bintable_row = self.bintables[extname][image.cube_index]
