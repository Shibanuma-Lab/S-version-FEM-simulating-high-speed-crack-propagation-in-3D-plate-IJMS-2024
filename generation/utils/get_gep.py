import numpy as np

def get_gep(first_coordinate, first_element_size, common_ratio, num_elem):
    """
    docstr
    """
    if common_ratio == 1 or common_ratio == 1.:        
        return first_element_size * np.arange(0, num_elem+1) + first_coordinate
    else:
        gep = [first_coordinate]
        for i in range(2, num_elem+2):
            gep.append(first_coordinate + (first_element_size * (1 - common_ratio**(i - 1))) / (1 - common_ratio))
        return gep