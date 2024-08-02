import numpy as np
from scipy.ndimage import distance_transform_edt

def distance_transform(image):
    distance_f = np.copy(image)
    distance_f = distance_transform_edt(distance_f)
    distance_b = np.copy(1-image)
    distance_b = distance_transform_edt(distance_b)
    signed_distance = distance_f - distance_b

    return signed_distance