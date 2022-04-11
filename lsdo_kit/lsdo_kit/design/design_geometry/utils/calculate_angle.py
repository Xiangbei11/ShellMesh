import numpy as np

def calculate_angle(v1, v2):

    angle = np.arccos( (np.dot(v1,v2)) /   (np.linalg.norm(v1) + np.linalg.norm(v2)))

    return angle
