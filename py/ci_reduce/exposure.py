import ci_reduce.common as common

class CI_exposure:
    """Object encapsulating the contents of a single CI exposure"""

    def __init__(self, image_list):
        # images is a dictionary of CI_image objects

        par = common.ci_misc_params()
        self.images = dict(zip(common.valid_image_extname_list(), 
                               par['n_cameras']*[None]))

        self.process_image_list(image_list)

    def assign_one_image(self, image):
        extname = (image.header)['EXTNAME']
        self.images[extname] = image

    def process_image_list(self, image_list):
        for image in image_list:
            self.assign_one_image(image)
