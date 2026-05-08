from cellpose import models

class NucleiSegmenter:
    def __init__(self):
        self.model = models.Cellpose()

    def segment(self, image):
        masks, *_ = self.model.eval(image, channels=[0, 0])
        return masks
