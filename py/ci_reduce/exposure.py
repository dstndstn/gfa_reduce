import ci_reduce.common as common
import ci_reduce.imred.load_calibs as load_calibs
import ci_reduce.dark_current as dark_current
import astropy.io.fits as fits

class CI_exposure:
    """Object encapsulating the contents of a single CI exposure"""

    def __init__(self, image_list, dummy_fz_header=None):
        # images is a dictionary of CI_image objects

        par = common.ci_misc_params()
        self.images = dict(zip(common.valid_image_extname_list(), 
                               par['n_cameras']*[None]))

        self.assign_image_list(image_list)
        self.dummy_fz_header = dummy_fz_header
        self.pixels_calibrated = None

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
                self.images[extname].bias_subtracted = True
                self.images[extname].calc_variance_e_squared()

    def subtract_dark_current(self):
        print('Attempting to subtract dark current...')
        for extname in self.images.keys():
            acttime = self.images[extname].header['EXPTIME']
            t_c = self.images[extname].header['CCDTEMP']
            self.images[extname].image = self.images[extname].image - \
                dark_current.total_dark_current_adu(extname, acttime, t_c)
            self.images[extname].dark_subtracted = True

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

    def create_all_bitmasks(self):
        for image in self.images.values():
            if image is None:
                continue
            print('Attempting to create image quality bitmask, ' + 
                  'extension name : ' + image.header['EXTNAME'])
            image.create_dq_mask()

    def to_hdulist(self, flavor=''):
        # I don't plan on writing fpack'ed outputs
        fz = False

        extnum_list = common.valid_image_extnum_list(fz=fz)

        hdulist = []
        for extnum in extnum_list:
            extname = common.ci_extnum_to_extname(extnum, fz=fz)
            hdu = self.images[extname].to_hdu(primary=(extnum == 0), 
                                              flavor=flavor)
            hdulist.append(hdu)

        return fits.HDUList(hdulist)

    def estimate_all_sky_mags(self, careful_sky=False):
        for im in self.images.values():
            im.estimate_sky_mag(careful_sky=careful_sky)

    def estimate_all_sky_sigmas(self, careful_sky=False):
        for im in self.images.values():
            im.set_empirical_bg_sigma(careful_sky=careful_sky)
            print('empirical ' + im.header['EXTNAME'] + 
                  ' background sigma = ' + 
                  '{:.2f}'.format(im.empirical_bg_sigma) + ' ADU')

    def all_source_catalogs(self):
        tables = dict(zip(self.images.keys(), len(self.images.keys())*[None]))

        for extname, im in self.images.items():
            tables[extname] = im.catalog_sources()

        # do I also want to store the tables as an attribute belonging to
        # this exposure object?
        return tables
