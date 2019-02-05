import ci_reduce.common as common

class CI_exposure:
    """Object encapsulating the contents of a single CI exposure"""

    def __init__(self, image_list):
        self.images = {1: None, 2: None, 3: None, 4: None, 5: None}

        self.process_image_list(image_list)

    def assign_one_image(self, image):
        extname = (image.header)['EXTNAME']
        self.images[common.ci_extname_to_ci_number(extname)] = image

    def process_image_list(self, image_list):
        for image in image_list:
            self.assign_one_image(image)
